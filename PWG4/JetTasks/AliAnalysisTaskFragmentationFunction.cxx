// *************************************************************************
// *                                                                       *
// * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
// *                                                                       *
// *************************************************************************


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

/* $Id: */

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom3.h"

#include "AliAODInputHandler.h" 
#include "AliAODHandler.h" 
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliVParticle.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskFragmentationFunction.h"

ClassImp(AliAnalysisTaskFragmentationFunction)

//____________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction()
   : AliAnalysisTaskSE()
   ,fESD(0)
   ,fAOD(0)
   ,fBranchRecJets("jets")
   ,fBranchRecBackJets("")
   ,fBranchRecBckgClusters("")
   ,fBranchGenJets("")
   ,fTrackTypeGen(0)
   ,fJetTypeGen(0)
   ,fJetTypeRecEff(0)
   ,fFilterMask(0)
   ,fUsePhysicsSelection(kTRUE)
   ,fEventClass(0)
   ,fMaxVertexZ(10)
   ,fTrackPtCut(0)
   ,fTrackEtaMin(0)
   ,fTrackEtaMax(0)
   ,fTrackPhiMin(0)
   ,fTrackPhiMax(0)
   ,fJetPtCut(0)
   ,fJetEtaMin(0)
   ,fJetEtaMax(0)
   ,fJetPhiMin(0)
   ,fJetPhiMax(0)
   ,fDiJetCut(0)
   ,fDiJetDeltaPhiCut(0)
   ,fDiJetPtFractionCut(0)
   ,fDiJetCDFCut(0)
   ,fDiJetKindBins(0)
   ,fFFRadius(0)
   ,fFFBckgRadius(0)
   ,fBckgMode(0)
   ,fIJMode(0)
   ,fQAMode(0)
   ,fFFMode(0)
   ,fDJMode(0)
   ,fEffMode(0)
   ,fPhiCorrMode(0)
   ,fUseRecEffRecJetPtBins(0)
   ,fUseResponseRecJetPtBins(1)
   ,fAvgTrials(0)
   ,fTracksRec(0)
   ,fTracksRecCuts(0)
   ,fTracksGen(0)
   ,fTracksAODMCCharged(0)
   ,fTracksRecQualityCuts(0)
   ,fJetsRec(0)
   ,fJetsRecCuts(0)
   ,fJetsGen(0)
   ,fJetsRecEff(0)
   ,fBckgJetsRec(0)
   ,fBckgJetsRecCuts(0)
   ,fBckgJetsGen(0)
   ,fQATrackHistosRec(0)
   ,fQATrackHistosRecCuts(0)
   ,fQATrackHistosGen(0)
   ,fQAJetHistosRec(0)
   ,fQAJetHistosRecCuts(0)
   ,fQAJetHistosRecCutsLeading(0)
   ,fQAJetHistosGen(0)
   ,fQAJetHistosGenLeading(0)
   ,fQAJetHistosRecEffLeading(0)
   ,fFFHistosRecCuts(0)
   ,fFFHistosRecLeading(0)
   ,fFFHistosRecLeadingTrack(0)
   ,fFFHistosGen(0)
   ,fFFHistosGenLeading(0)
   ,fFFHistosGenLeadingTrack(0)
   ,fIJHistosRecCuts(0)
   ,fIJHistosRecLeading(0)
   ,fIJHistosRecLeadingTrack(0)
   ,fIJHistosGen(0)
   ,fIJHistosGenLeading(0)
   ,fIJHistosGenLeadingTrack(0)
   ,fFFDiJetHistosRecCuts(0)
   ,fFFDiJetHistosRecLeading(0)
   ,fFFDiJetHistosRecLeadingTrack(0)
   ,fFFDiJetHistosGen(0)
   ,fFFDiJetHistosGenLeading(0)
   ,fFFDiJetHistosGenLeadingTrack(0)
   ,fQADiJetHistosRecCuts(0)
   ,fQADiJetHistosGen(0)
   ,fPhiCorrHistosJetArea(0)
   ,fPhiCorrHistosTransverseArea(0)
   ,fPhiCorrHistosAwayArea(0)
   ,fQATrackHighPtThreshold(0)
   ,fFFNBinsJetPt(0)    
   ,fFFJetPtMin(0) 
   ,fFFJetPtMax(0)
   ,fFFNBinsPt(0)      
   ,fFFPtMin(0)        
   ,fFFPtMax(0)        
   ,fFFNBinsXi(0)      
   ,fFFXiMin(0)        
   ,fFFXiMax(0)        
   ,fFFNBinsZ(0)       
   ,fFFZMin(0)         
   ,fFFZMax(0)         
   ,fQAJetNBinsPt(0)   
   ,fQAJetPtMin(0)     
   ,fQAJetPtMax(0)     
   ,fQAJetNBinsEta(0)  
   ,fQAJetEtaMin(0)    
   ,fQAJetEtaMax(0)    
   ,fQAJetNBinsPhi(0)  
   ,fQAJetPhiMin(0)    
   ,fQAJetPhiMax(0)    
   ,fQATrackNBinsPt(0) 
   ,fQATrackPtMin(0)   
   ,fQATrackPtMax(0)   
   ,fQATrackNBinsEta(0)
   ,fQATrackEtaMin(0)  
   ,fQATrackEtaMax(0)  
   ,fQATrackNBinsPhi(0)
   ,fQATrackPhiMin(0)  
   ,fQATrackPhiMax(0)
   ,fIJNBinsJetPt(0)
   ,fIJJetPtMin(0)
   ,fIJJetPtMax(0)
   ,fIJNBinsPt(0)
   ,fIJPtMin(0)
   ,fIJPtMax(0)
   ,fIJNBinsZ(0)
   ,fIJZMin(0)
   ,fIJZMax(0)
   ,fIJNBinsCosTheta(0)
   ,fIJCosThetaMin(0)
   ,fIJCosThetaMax(0)
   ,fIJNBinsTheta(0)
   ,fIJThetaMin(0)
   ,fIJThetaMax(0)
   ,fIJNBinsJt(0)
   ,fIJJtMin(0)
   ,fIJJtMax(0)
   ,fDiJetNBinsJetInvMass(0)
   ,fDiJetJetInvMassMin(0)
   ,fDiJetJetInvMassMax(0)
   ,fDiJetNBinsJetPt(0)
   ,fDiJetJetPtMin(0)
   ,fDiJetJetPtMax(0)
   ,fDiJetNBinsPt(0)
   ,fDiJetPtMin(0)
   ,fDiJetPtMax(0)
   ,fDiJetNBinsXi(0)
   ,fDiJetXiMin(0)
   ,fDiJetXiMax(0)
   ,fDiJetNBinsZ(0)
   ,fDiJetZMin(0)
   ,fDiJetZMax(0)
   ,fQADiJetNBinsInvMass(0)
   ,fQADiJetInvMassMin(0)
   ,fQADiJetInvMassMax(0)
   ,fQADiJetNBinsJetPt(0)
   ,fQADiJetJetPtMin(0)
   ,fQADiJetJetPtMax(0)
   ,fQADiJetNBinsDeltaPhi(0)
   ,fQADiJetDeltaPhiMin(0)
   ,fQADiJetDeltaPhiMax(0)
   ,fQADiJetNBinsDeltaEta(0)
   ,fQADiJetDeltaEtaMin(0)
   ,fQADiJetDeltaEtaMax(0)
   ,fQADiJetNBinsDeltaPt(0)
   ,fQADiJetDeltaPtMin(0)
   ,fQADiJetDeltaPtMax(0)
   ,fQADiJetNBinsInBal(0)  
   ,fQADiJetInBalMin(0)    
   ,fQADiJetInBalMax(0)    
   ,fPhiCorrNBinsPt(0)
   ,fPhiCorrPtMin(0)
   ,fPhiCorrPtMax(0)
   ,fPhiCorrNBinsEta(0)
   ,fPhiCorrEtaMin(0)
   ,fPhiCorrEtaMax(0)
   ,fPhiCorrNBinsPhi(0)
   ,fPhiCorrPhiMin(0)
   ,fPhiCorrPhiMax(0)
   ,fCommonHistList(0)
   ,fh1EvtSelection(0)
   ,fh1VertexNContributors(0)
   ,fh1VertexZ(0)
   ,fh1EvtMult(0)
   ,fh1EvtCent(0)
   ,fh1Xsec(0)
   ,fh1Trials(0)
   ,fh1PtHard(0)
   ,fh1PtHardTrials(0)
   ,fh1nRecJetsCuts(0)
   ,fh1nGenJets(0)
   ,fh1nRecEffJets(0)
   ,fh1nRecBckgJetsCuts(0)
   ,fh1nGenBckgJets(0)
   ,fh2PtRecVsGenPrim(0)
   ,fQATrackHistosRecEffGen(0)  
   ,fQATrackHistosRecEffRec(0)  
   ,fFFHistosRecEffGen(0)    
   ,fFFHistosRecEffRec(0)
   ,fhnResponseSinglePt(0)  
   ,fhnResponseJetTrackPt(0)  
   ,fhnResponseJetZ(0)        
   ,fhnResponseJetXi(0)       
   // Background
   ,fh1OutLeadingMult(0)
   ,fh1OutLeadingStatMult(0)
   ,fh1PerpMult(0)
   ,fh1ASideMult(0)
   ,fh1ASideWindowMult(0)
   ,fh1PerpWindowMult(0)
   ,fh1Out2JetsMult(0)
   ,fh1Out3JetsMult(0)
   ,fh1MedianClustersMult(0)
   ,fh1OutClustersMult(0)
   ,fQABckgHisto0RecCuts(0)  
   ,fQABckgHisto0Gen(0)      
   ,fQABckgHisto1RecCuts(0)  
   ,fQABckgHisto1Gen(0)      
   ,fQABckgHisto2RecCuts(0)  
   ,fQABckgHisto2Gen(0)
   ,fQABckgHisto3RecCuts(0)
   ,fQABckgHisto3Gen(0)
   ,fFFBckgHisto0RecCuts(0)
   ,fFFBckgHisto0RecLeading(0)
   ,fFFBckgHisto0Gen(0)       
   ,fFFBckgHisto0GenLeading(0)
   ,fFFBckgHisto1RecCuts(0)
   ,fFFBckgHisto1RecLeading(0)
   ,fFFBckgHisto1Gen(0)       
   ,fFFBckgHisto1GenLeading(0)
   ,fFFBckgHisto2RecCuts(0)
   ,fFFBckgHisto2RecLeading(0)
   ,fFFBckgHisto2Gen(0)       
   ,fFFBckgHisto2GenLeading(0)
   ,fFFBckgHisto3RecCuts(0)
   ,fFFBckgHisto3RecLeading(0)
   ,fFFBckgHisto3Gen(0)       
   ,fFFBckgHisto3GenLeading(0)
   ,fIJBckgHisto0RecCuts(0)   
   ,fIJBckgHisto0RecLeading(0)
   ,fIJBckgHisto0Gen(0)       
   ,fIJBckgHisto0GenLeading(0)
   ,fIJBckgHisto1RecCuts(0)   
   ,fIJBckgHisto1RecLeading(0)
   ,fIJBckgHisto1Gen(0)       
   ,fIJBckgHisto1GenLeading(0)
   ,fIJBckgHisto2RecCuts(0)   
   ,fIJBckgHisto2RecLeading(0)
   ,fIJBckgHisto2Gen(0)       
   ,fIJBckgHisto2GenLeading(0)
   ,fIJBckgHisto3RecCuts(0)   
   ,fIJBckgHisto3RecLeading(0)
   ,fIJBckgHisto3Gen(0)       
   ,fIJBckgHisto3GenLeading(0)
   ,fRandom(0)
   ,fBckgSubMethod(0)
{
   // default constructor
  fBckgType[0] = 0;
  fBckgType[1] = 0;
  fBckgType[2] = 0;
  fBckgType[3] = 0;
}

//__________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fBranchRecJets("jets")
  ,fBranchRecBackJets("")
  ,fBranchRecBckgClusters("")
  ,fBranchGenJets("")
  ,fTrackTypeGen(0)
  ,fJetTypeGen(0)
  ,fJetTypeRecEff(0)
  ,fFilterMask(0)
  ,fUsePhysicsSelection(kTRUE)
  ,fEventClass(0)
  ,fMaxVertexZ(10)
  ,fTrackPtCut(0)
  ,fTrackEtaMin(0)
  ,fTrackEtaMax(0)
  ,fTrackPhiMin(0)
  ,fTrackPhiMax(0)
  ,fJetPtCut(0)
  ,fJetEtaMin(0)
  ,fJetEtaMax(0)
  ,fJetPhiMin(0)
  ,fJetPhiMax(0)
  ,fDiJetCut(0)
  ,fDiJetDeltaPhiCut(0)
  ,fDiJetPtFractionCut(0)
  ,fDiJetCDFCut(0)
  ,fDiJetKindBins(0)
  ,fFFRadius(0)
  ,fFFBckgRadius(0)
  ,fBckgMode(0)
  ,fIJMode(0)
  ,fQAMode(0)
  ,fFFMode(0)
  ,fDJMode(0)
  ,fEffMode(0)
  ,fPhiCorrMode(0)
  ,fUseRecEffRecJetPtBins(0)
  ,fUseResponseRecJetPtBins(1)
  ,fAvgTrials(0)
  ,fTracksRec(0)
  ,fTracksRecCuts(0)
  ,fTracksGen(0)
  ,fTracksAODMCCharged(0)
  ,fTracksRecQualityCuts(0)
  ,fJetsRec(0)
  ,fJetsRecCuts(0)
  ,fJetsGen(0)
  ,fJetsRecEff(0)
  ,fBckgJetsRec(0)
  ,fBckgJetsRecCuts(0)
  ,fBckgJetsGen(0)
  ,fQATrackHistosRec(0)
  ,fQATrackHistosRecCuts(0)
  ,fQATrackHistosGen(0)
  ,fQAJetHistosRec(0)
  ,fQAJetHistosRecCuts(0)
  ,fQAJetHistosRecCutsLeading(0)
  ,fQAJetHistosGen(0)
  ,fQAJetHistosGenLeading(0)
  ,fQAJetHistosRecEffLeading(0)
  ,fFFHistosRecCuts(0)
  ,fFFHistosRecLeading(0)
  ,fFFHistosRecLeadingTrack(0)
  ,fFFHistosGen(0)
  ,fFFHistosGenLeading(0)
  ,fFFHistosGenLeadingTrack(0)
  ,fIJHistosRecCuts(0)
  ,fIJHistosRecLeading(0)
  ,fIJHistosRecLeadingTrack(0)
  ,fIJHistosGen(0)
  ,fIJHistosGenLeading(0)
  ,fIJHistosGenLeadingTrack(0)
  ,fFFDiJetHistosRecCuts(0)
  ,fFFDiJetHistosRecLeading(0)
  ,fFFDiJetHistosRecLeadingTrack(0)
  ,fFFDiJetHistosGen(0)
  ,fFFDiJetHistosGenLeading(0)
  ,fFFDiJetHistosGenLeadingTrack(0)
  ,fQADiJetHistosRecCuts(0)
  ,fQADiJetHistosGen(0)
  ,fPhiCorrHistosJetArea(0)
  ,fPhiCorrHistosTransverseArea(0)
  ,fPhiCorrHistosAwayArea(0)
  ,fQATrackHighPtThreshold(0) 
  ,fFFNBinsJetPt(0)    
  ,fFFJetPtMin(0) 
  ,fFFJetPtMax(0)
  ,fFFNBinsPt(0)      
  ,fFFPtMin(0)        
  ,fFFPtMax(0)        
  ,fFFNBinsXi(0)      
  ,fFFXiMin(0)        
  ,fFFXiMax(0)        
  ,fFFNBinsZ(0)       
  ,fFFZMin(0)         
  ,fFFZMax(0)         
  ,fQAJetNBinsPt(0)   
  ,fQAJetPtMin(0)     
  ,fQAJetPtMax(0)     
  ,fQAJetNBinsEta(0)  
  ,fQAJetEtaMin(0)    
  ,fQAJetEtaMax(0)    
  ,fQAJetNBinsPhi(0)  
  ,fQAJetPhiMin(0)    
  ,fQAJetPhiMax(0)    
  ,fQATrackNBinsPt(0) 
  ,fQATrackPtMin(0)   
  ,fQATrackPtMax(0)   
  ,fQATrackNBinsEta(0)
  ,fQATrackEtaMin(0)  
  ,fQATrackEtaMax(0)  
  ,fQATrackNBinsPhi(0)
  ,fQATrackPhiMin(0)  
  ,fQATrackPhiMax(0)  
  ,fIJNBinsJetPt(0)
  ,fIJJetPtMin(0)
  ,fIJJetPtMax(0)
  ,fIJNBinsPt(0)
  ,fIJPtMin(0)
  ,fIJPtMax(0)
  ,fIJNBinsZ(0)
  ,fIJZMin(0)
  ,fIJZMax(0)
  ,fIJNBinsCosTheta(0)
  ,fIJCosThetaMin(0)
  ,fIJCosThetaMax(0)
  ,fIJNBinsTheta(0)
  ,fIJThetaMin(0)
  ,fIJThetaMax(0)
  ,fIJNBinsJt(0)
  ,fIJJtMin(0)
  ,fIJJtMax(0)
  ,fDiJetNBinsJetInvMass(0)
  ,fDiJetJetInvMassMin(0)
  ,fDiJetJetInvMassMax(0)
  ,fDiJetNBinsJetPt(0)
  ,fDiJetJetPtMin(0)
  ,fDiJetJetPtMax(0)
  ,fDiJetNBinsPt(0)
  ,fDiJetPtMin(0)
  ,fDiJetPtMax(0)
  ,fDiJetNBinsXi(0)
  ,fDiJetXiMin(0)
  ,fDiJetXiMax(0)
  ,fDiJetNBinsZ(0)
  ,fDiJetZMin(0)
  ,fDiJetZMax(0)
  ,fQADiJetNBinsInvMass(0)
  ,fQADiJetInvMassMin(0)
  ,fQADiJetInvMassMax(0)
  ,fQADiJetNBinsJetPt(0)
  ,fQADiJetJetPtMin(0)
  ,fQADiJetJetPtMax(0)
  ,fQADiJetNBinsDeltaPhi(0)
  ,fQADiJetDeltaPhiMin(0)
  ,fQADiJetDeltaPhiMax(0)
  ,fQADiJetNBinsDeltaEta(0)
  ,fQADiJetDeltaEtaMin(0)
  ,fQADiJetDeltaEtaMax(0)
  ,fQADiJetNBinsDeltaPt(0)
  ,fQADiJetDeltaPtMin(0)
  ,fQADiJetDeltaPtMax(0)
  ,fQADiJetNBinsInBal(0)  
  ,fQADiJetInBalMin(0)    
  ,fQADiJetInBalMax(0)    
  ,fPhiCorrNBinsPt(0)
  ,fPhiCorrPtMin(0)
  ,fPhiCorrPtMax(0)
  ,fPhiCorrNBinsEta(0)
  ,fPhiCorrEtaMin(0)
  ,fPhiCorrEtaMax(0)
  ,fPhiCorrNBinsPhi(0)
  ,fPhiCorrPhiMin(0)
  ,fPhiCorrPhiMax(0)
  ,fCommonHistList(0)
  ,fh1EvtSelection(0)
  ,fh1VertexNContributors(0)
  ,fh1VertexZ(0)
  ,fh1EvtMult(0)
  ,fh1EvtCent(0)
  ,fh1Xsec(0)
  ,fh1Trials(0)
  ,fh1PtHard(0)
  ,fh1PtHardTrials(0)
  ,fh1nRecJetsCuts(0)
  ,fh1nGenJets(0)
  ,fh1nRecEffJets(0)
  ,fh1nRecBckgJetsCuts(0)
  ,fh1nGenBckgJets(0)
  ,fh2PtRecVsGenPrim(0)
  ,fQATrackHistosRecEffGen(0)  
  ,fQATrackHistosRecEffRec(0)  
  ,fFFHistosRecEffGen(0)    
  ,fFFHistosRecEffRec(0)
  ,fhnResponseSinglePt(0)  
  ,fhnResponseJetTrackPt(0)  
  ,fhnResponseJetZ(0)        
  ,fhnResponseJetXi(0)       
  // Background
  ,fh1OutLeadingMult(0)
  ,fh1OutLeadingStatMult(0)
  ,fh1PerpMult(0)
  ,fh1ASideMult(0)
  ,fh1ASideWindowMult(0)
  ,fh1PerpWindowMult(0)
  ,fh1Out2JetsMult(0)
  ,fh1Out3JetsMult(0)
  ,fh1MedianClustersMult(0)
  ,fh1OutClustersMult(0)
  ,fQABckgHisto0RecCuts(0)  
  ,fQABckgHisto0Gen(0)      
  ,fQABckgHisto1RecCuts(0)  
  ,fQABckgHisto1Gen(0)      
  ,fQABckgHisto2RecCuts(0)  
  ,fQABckgHisto2Gen(0) 
  ,fQABckgHisto3RecCuts(0)  
  ,fQABckgHisto3Gen(0)
  ,fFFBckgHisto0RecCuts(0)
  ,fFFBckgHisto0RecLeading(0)
  ,fFFBckgHisto0Gen(0)       
  ,fFFBckgHisto0GenLeading(0)
  ,fFFBckgHisto1RecCuts(0)
  ,fFFBckgHisto1RecLeading(0)
  ,fFFBckgHisto1Gen(0)       
  ,fFFBckgHisto1GenLeading(0)
  ,fFFBckgHisto2RecCuts(0)
  ,fFFBckgHisto2RecLeading(0)
  ,fFFBckgHisto2Gen(0)       
  ,fFFBckgHisto2GenLeading(0)
  ,fFFBckgHisto3RecCuts(0)
  ,fFFBckgHisto3RecLeading(0)
  ,fFFBckgHisto3Gen(0)       
  ,fFFBckgHisto3GenLeading(0)
  ,fIJBckgHisto0RecCuts(0)   
  ,fIJBckgHisto0RecLeading(0)
  ,fIJBckgHisto0Gen(0)       
  ,fIJBckgHisto0GenLeading(0)
  ,fIJBckgHisto1RecCuts(0)   
  ,fIJBckgHisto1RecLeading(0)
  ,fIJBckgHisto1Gen(0)       
  ,fIJBckgHisto1GenLeading(0)
  ,fIJBckgHisto2RecCuts(0)   
  ,fIJBckgHisto2RecLeading(0)
  ,fIJBckgHisto2Gen(0)       
  ,fIJBckgHisto2GenLeading(0)
  ,fIJBckgHisto3RecCuts(0)   
  ,fIJBckgHisto3RecLeading(0)
  ,fIJBckgHisto3Gen(0)       
  ,fIJBckgHisto3GenLeading(0)
  ,fRandom(0)
  ,fBckgSubMethod(0)
{
  // constructor
  fBckgType[0] = 0;
  fBckgType[1] = 0;
  fBckgType[2] = 0;
  fBckgType[3] = 0;

  DefineOutput(1,TList::Class());
  

}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const  AliAnalysisTaskFragmentationFunction &copy)
  : AliAnalysisTaskSE()
  ,fESD(copy.fESD)
  ,fAOD(copy.fAOD)
  ,fBranchRecJets(copy.fBranchRecJets)
  ,fBranchRecBackJets(copy.fBranchRecBackJets)
  ,fBranchRecBckgClusters(copy.fBranchRecBckgClusters)
  ,fBranchGenJets(copy.fBranchGenJets)
  ,fTrackTypeGen(copy.fTrackTypeGen)
  ,fJetTypeGen(copy.fJetTypeGen)
  ,fJetTypeRecEff(copy.fJetTypeRecEff)
  ,fFilterMask(copy.fFilterMask)
  ,fUsePhysicsSelection(copy.fUsePhysicsSelection)
  ,fEventClass(copy.fEventClass)
  ,fMaxVertexZ(copy.fMaxVertexZ)
  ,fTrackPtCut(copy.fTrackPtCut)
  ,fTrackEtaMin(copy.fTrackEtaMin)
  ,fTrackEtaMax(copy.fTrackEtaMax)
  ,fTrackPhiMin(copy.fTrackPhiMin)
  ,fTrackPhiMax(copy.fTrackPhiMax)
  ,fJetPtCut(copy.fJetPtCut)
  ,fJetEtaMin(copy.fJetEtaMin)
  ,fJetEtaMax(copy.fJetEtaMax)
  ,fJetPhiMin(copy.fJetPhiMin)
  ,fJetPhiMax(copy.fJetPhiMax)
  ,fDiJetCut(copy.fDiJetCut)
  ,fDiJetDeltaPhiCut(copy.fDiJetDeltaPhiCut)
  ,fDiJetPtFractionCut(copy.fDiJetPtFractionCut)
  ,fDiJetCDFCut(copy.fDiJetCDFCut)
  ,fDiJetKindBins(copy.fDiJetKindBins)
  ,fFFRadius(copy.fFFRadius)
  ,fFFBckgRadius(copy.fFFBckgRadius)
  ,fBckgMode(copy.fBckgMode)
  ,fIJMode(copy.fIJMode)
  ,fQAMode(copy.fQAMode)
  ,fFFMode(copy.fFFMode)
  ,fDJMode(copy.fDJMode)
  ,fEffMode(copy.fEffMode)
  ,fPhiCorrMode(copy.fPhiCorrMode)
  ,fUseRecEffRecJetPtBins(copy.fUseRecEffRecJetPtBins)
  ,fUseResponseRecJetPtBins(copy.fUseResponseRecJetPtBins)
  ,fAvgTrials(copy.fAvgTrials)
  ,fTracksRec(copy.fTracksRec)
  ,fTracksRecCuts(copy.fTracksRecCuts)
  ,fTracksGen(copy.fTracksGen)
  ,fTracksAODMCCharged(copy.fTracksAODMCCharged)
  ,fTracksRecQualityCuts(copy.fTracksRecQualityCuts)
  ,fJetsRec(copy.fJetsRec)
  ,fJetsRecCuts(copy.fJetsRecCuts)
  ,fJetsGen(copy.fJetsGen)
  ,fJetsRecEff(copy.fJetsRecEff)
  ,fBckgJetsRec(copy.fBckgJetsRec)
  ,fBckgJetsRecCuts(copy.fBckgJetsRecCuts)
  ,fBckgJetsGen(copy.fBckgJetsGen)
  ,fQATrackHistosRec(copy.fQATrackHistosRec)
  ,fQATrackHistosRecCuts(copy.fQATrackHistosRecCuts)
  ,fQATrackHistosGen(copy.fQATrackHistosGen)
  ,fQAJetHistosRec(copy.fQAJetHistosRec)
  ,fQAJetHistosRecCuts(copy.fQAJetHistosRecCuts)
  ,fQAJetHistosRecCutsLeading(copy.fQAJetHistosRecCutsLeading)
  ,fQAJetHistosGen(copy.fQAJetHistosGen)
  ,fQAJetHistosGenLeading(copy.fQAJetHistosGenLeading)
  ,fQAJetHistosRecEffLeading(copy.fQAJetHistosRecEffLeading)
  ,fFFHistosRecCuts(copy.fFFHistosRecCuts)
  ,fFFHistosRecLeading(copy.fFFHistosRecLeading)
  ,fFFHistosRecLeadingTrack(copy.fFFHistosRecLeadingTrack)
  ,fFFHistosGen(copy.fFFHistosGen)
  ,fFFHistosGenLeading(copy.fFFHistosGenLeading)
  ,fFFHistosGenLeadingTrack(copy.fFFHistosGenLeadingTrack)
  ,fIJHistosRecCuts(copy.fIJHistosRecCuts)
  ,fIJHistosRecLeading(copy.fIJHistosRecLeading)
  ,fIJHistosRecLeadingTrack(copy.fIJHistosRecLeadingTrack)
  ,fIJHistosGen(copy.fIJHistosGen)
  ,fIJHistosGenLeading(copy.fIJHistosGenLeading)
  ,fIJHistosGenLeadingTrack(copy.fIJHistosGenLeadingTrack)
  ,fFFDiJetHistosRecCuts(copy.fFFDiJetHistosRecCuts)
  ,fFFDiJetHistosRecLeading(copy.fFFDiJetHistosRecLeading)
  ,fFFDiJetHistosRecLeadingTrack(copy.fFFDiJetHistosRecLeadingTrack)
  ,fFFDiJetHistosGen(copy.fFFDiJetHistosGen)
  ,fFFDiJetHistosGenLeading(copy.fFFDiJetHistosGenLeading)
  ,fFFDiJetHistosGenLeadingTrack(copy.fFFDiJetHistosGenLeadingTrack)
  ,fQADiJetHistosRecCuts(copy.fQADiJetHistosRecCuts)
  ,fQADiJetHistosGen(copy.fQADiJetHistosGen)
  ,fPhiCorrHistosJetArea(copy.fPhiCorrHistosJetArea)
  ,fPhiCorrHistosTransverseArea(copy.fPhiCorrHistosTransverseArea)
  ,fPhiCorrHistosAwayArea(copy.fPhiCorrHistosAwayArea)
  ,fQATrackHighPtThreshold(copy.fQATrackHighPtThreshold) 
  ,fFFNBinsJetPt(copy.fFFNBinsJetPt)    
  ,fFFJetPtMin(copy.fFFJetPtMin) 
  ,fFFJetPtMax(copy.fFFJetPtMax)
  ,fFFNBinsPt(copy.fFFNBinsPt)      
  ,fFFPtMin(copy.fFFPtMin)        
  ,fFFPtMax(copy.fFFPtMax)        
  ,fFFNBinsXi(copy.fFFNBinsXi)      
  ,fFFXiMin(copy.fFFXiMin)        
  ,fFFXiMax(copy.fFFXiMax)        
  ,fFFNBinsZ(copy.fFFNBinsZ)       
  ,fFFZMin(copy.fFFZMin)         
  ,fFFZMax(copy.fFFZMax)         
  ,fQAJetNBinsPt(copy.fQAJetNBinsPt)   
  ,fQAJetPtMin(copy.fQAJetPtMin)     
  ,fQAJetPtMax(copy.fQAJetPtMax)     
  ,fQAJetNBinsEta(copy.fQAJetNBinsEta)  
  ,fQAJetEtaMin(copy.fQAJetEtaMin)    
  ,fQAJetEtaMax(copy.fQAJetEtaMax)    
  ,fQAJetNBinsPhi(copy.fQAJetNBinsPhi)  
  ,fQAJetPhiMin(copy.fQAJetPhiMin)    
  ,fQAJetPhiMax(copy.fQAJetPhiMax)    
  ,fQATrackNBinsPt(copy.fQATrackNBinsPt) 
  ,fQATrackPtMin(copy.fQATrackPtMin)   
  ,fQATrackPtMax(copy.fQATrackPtMax)   
  ,fQATrackNBinsEta(copy.fQATrackNBinsEta)
  ,fQATrackEtaMin(copy.fQATrackEtaMin)  
  ,fQATrackEtaMax(copy.fQATrackEtaMax)  
  ,fQATrackNBinsPhi(copy.fQATrackNBinsPhi)
  ,fQATrackPhiMin(copy.fQATrackPhiMin)  
  ,fQATrackPhiMax(copy.fQATrackPhiMax)
  ,fIJNBinsJetPt(copy.fIJNBinsJetPt)
  ,fIJJetPtMin(copy.fIJJetPtMin)
  ,fIJJetPtMax(copy.fIJJetPtMax)
  ,fIJNBinsPt(copy.fIJNBinsPt)
  ,fIJPtMin(copy.fIJPtMin)
  ,fIJPtMax(copy.fIJPtMax)
  ,fIJNBinsZ(copy.fIJNBinsZ)
  ,fIJZMin(copy.fIJZMin)
  ,fIJZMax(copy.fIJZMax)
  ,fIJNBinsCosTheta(copy.fIJNBinsCosTheta)
  ,fIJCosThetaMin(copy.fIJCosThetaMin)
  ,fIJCosThetaMax(copy.fIJCosThetaMax)
  ,fIJNBinsTheta(copy.fIJNBinsTheta)
  ,fIJThetaMin(copy.fIJThetaMin)
  ,fIJThetaMax(copy.fIJThetaMax)
  ,fIJNBinsJt(copy.fIJNBinsJt)
  ,fIJJtMin(copy.fIJJtMin)
  ,fIJJtMax(copy.fIJJtMax)
  ,fDiJetNBinsJetInvMass(copy.fDiJetNBinsJetInvMass)
  ,fDiJetJetInvMassMin(copy.fDiJetJetInvMassMin)
  ,fDiJetJetInvMassMax(copy.fDiJetJetInvMassMax)
  ,fDiJetNBinsJetPt(copy.fDiJetNBinsJetPt)
  ,fDiJetJetPtMin(copy.fDiJetJetPtMin)
  ,fDiJetJetPtMax(copy.fDiJetJetPtMax)
  ,fDiJetNBinsPt(copy.fDiJetNBinsPt)
  ,fDiJetPtMin(copy.fDiJetPtMin)
  ,fDiJetPtMax(copy.fDiJetPtMax)
  ,fDiJetNBinsXi(copy.fDiJetNBinsXi)
  ,fDiJetXiMin(copy.fDiJetXiMin)
  ,fDiJetXiMax(copy.fDiJetXiMax)
  ,fDiJetNBinsZ(copy.fDiJetNBinsZ)
  ,fDiJetZMin(copy.fDiJetZMin)
  ,fDiJetZMax(copy.fDiJetZMax)
  ,fQADiJetNBinsInvMass(copy.fQADiJetNBinsInvMass)
  ,fQADiJetInvMassMin(copy.fQADiJetInvMassMin)
  ,fQADiJetInvMassMax(copy.fQADiJetInvMassMax)
  ,fQADiJetNBinsJetPt(copy.fQADiJetNBinsJetPt)
  ,fQADiJetJetPtMin(copy.fQADiJetJetPtMin)
  ,fQADiJetJetPtMax(copy.fQADiJetJetPtMax)
  ,fQADiJetNBinsDeltaPhi(copy.fQADiJetNBinsDeltaPhi)
  ,fQADiJetDeltaPhiMin(copy.fQADiJetDeltaPhiMin)
  ,fQADiJetDeltaPhiMax(copy.fQADiJetDeltaPhiMax)
  ,fQADiJetNBinsDeltaEta(copy.fQADiJetNBinsDeltaEta)
  ,fQADiJetDeltaEtaMin(copy.fQADiJetDeltaEtaMin)
  ,fQADiJetDeltaEtaMax(copy.fQADiJetDeltaEtaMax)
  ,fQADiJetNBinsDeltaPt(copy.fQADiJetNBinsDeltaPt)
  ,fQADiJetDeltaPtMin(copy.fQADiJetDeltaPtMin)
  ,fQADiJetDeltaPtMax(copy.fQADiJetDeltaPtMax)
  ,fQADiJetNBinsInBal(copy.fQADiJetNBinsInBal)
  ,fQADiJetInBalMin(copy.fQADiJetInBalMin)
  ,fQADiJetInBalMax(copy.fQADiJetInBalMax)
  ,fPhiCorrNBinsPt(copy.fPhiCorrNBinsPt)
  ,fPhiCorrPtMin(copy.fPhiCorrPtMin)
  ,fPhiCorrPtMax(copy.fPhiCorrPtMax)
  ,fPhiCorrNBinsEta(copy.fPhiCorrNBinsEta)
  ,fPhiCorrEtaMin(copy.fPhiCorrEtaMin)
  ,fPhiCorrEtaMax(copy.fPhiCorrEtaMax)
  ,fPhiCorrNBinsPhi(copy.fPhiCorrNBinsPhi)
  ,fPhiCorrPhiMin(copy.fPhiCorrPhiMin)
  ,fPhiCorrPhiMax(copy.fPhiCorrPhiMax)
  ,fCommonHistList(copy.fCommonHistList)
  ,fh1EvtSelection(copy.fh1EvtSelection)
  ,fh1VertexNContributors(copy.fh1VertexNContributors)
  ,fh1VertexZ(copy.fh1VertexZ)
  ,fh1EvtMult(copy.fh1EvtMult)
  ,fh1EvtCent(copy.fh1EvtCent)
  ,fh1Xsec(copy.fh1Xsec)
  ,fh1Trials(copy.fh1Trials)
  ,fh1PtHard(copy.fh1PtHard)  
  ,fh1PtHardTrials(copy.fh1PtHardTrials)  
  ,fh1nRecJetsCuts(copy.fh1nRecJetsCuts)
  ,fh1nGenJets(copy.fh1nGenJets)
  ,fh1nRecEffJets(copy.fh1nRecEffJets)
  ,fh1nRecBckgJetsCuts(copy.fh1nRecBckgJetsCuts)
  ,fh1nGenBckgJets(copy.fh1nGenBckgJets)
  ,fh2PtRecVsGenPrim(copy.fh2PtRecVsGenPrim)
  ,fQATrackHistosRecEffGen(copy.fQATrackHistosRecEffGen)  
  ,fQATrackHistosRecEffRec(copy.fQATrackHistosRecEffRec)  
  ,fFFHistosRecEffGen(copy.fFFHistosRecEffGen)    
  ,fFFHistosRecEffRec(copy.fFFHistosRecEffRec)   
  ,fhnResponseSinglePt(copy.fhnResponseSinglePt)
  ,fhnResponseJetTrackPt(copy.fhnResponseJetTrackPt)
  ,fhnResponseJetZ(copy.fhnResponseJetZ)
  ,fhnResponseJetXi(copy.fhnResponseJetXi)
  // Background
  ,fh1OutLeadingMult(copy.fh1OutLeadingMult)
  ,fh1OutLeadingStatMult(copy.fh1OutLeadingStatMult)
  ,fh1PerpMult(copy.fh1PerpMult)
  ,fh1ASideMult(copy.fh1ASideMult)
  ,fh1ASideWindowMult(copy.fh1ASideWindowMult)
  ,fh1PerpWindowMult(copy.fh1PerpWindowMult)
  ,fh1Out2JetsMult(copy.fh1Out2JetsMult)
  ,fh1Out3JetsMult(copy.fh1Out3JetsMult)
  ,fh1MedianClustersMult(copy.fh1MedianClustersMult)
  ,fh1OutClustersMult(copy.fh1OutClustersMult)
  ,fQABckgHisto0RecCuts(copy.fQABckgHisto0RecCuts)  
  ,fQABckgHisto0Gen(copy.fQABckgHisto0Gen)      
  ,fQABckgHisto1RecCuts(copy.fQABckgHisto1RecCuts)  
  ,fQABckgHisto1Gen(copy.fQABckgHisto1Gen)      
  ,fQABckgHisto2RecCuts(copy.fQABckgHisto2RecCuts)  
  ,fQABckgHisto2Gen(copy.fQABckgHisto2Gen)
  ,fQABckgHisto3RecCuts(copy.fQABckgHisto3RecCuts)  
  ,fQABckgHisto3Gen(copy.fQABckgHisto3Gen)     
  ,fFFBckgHisto0RecCuts(copy.fFFBckgHisto0RecCuts)
  ,fFFBckgHisto0RecLeading(copy.fFFBckgHisto0RecLeading)
  ,fFFBckgHisto0Gen(copy.fFFBckgHisto0Gen)       
  ,fFFBckgHisto0GenLeading(copy.fFFBckgHisto0GenLeading)
  ,fFFBckgHisto1RecCuts(copy.fFFBckgHisto1RecCuts)
  ,fFFBckgHisto1RecLeading(copy.fFFBckgHisto1RecLeading)
  ,fFFBckgHisto1Gen(copy.fFFBckgHisto1Gen)       
  ,fFFBckgHisto1GenLeading(copy.fFFBckgHisto1GenLeading)
  ,fFFBckgHisto2RecCuts(copy.fFFBckgHisto2RecCuts)
  ,fFFBckgHisto2RecLeading(copy.fFFBckgHisto2RecLeading)
  ,fFFBckgHisto2Gen(copy.fFFBckgHisto2Gen)       
  ,fFFBckgHisto2GenLeading(copy.fFFBckgHisto2GenLeading)
  ,fFFBckgHisto3RecCuts(copy.fFFBckgHisto3RecCuts)
  ,fFFBckgHisto3RecLeading(copy.fFFBckgHisto3RecLeading)
  ,fFFBckgHisto3Gen(copy.fFFBckgHisto3Gen)       
  ,fFFBckgHisto3GenLeading(copy.fFFBckgHisto3GenLeading)
  ,fIJBckgHisto0RecCuts(copy.fIJBckgHisto0RecCuts)   
  ,fIJBckgHisto0RecLeading(copy.fIJBckgHisto0RecLeading)
  ,fIJBckgHisto0Gen(copy.fIJBckgHisto0Gen)       
  ,fIJBckgHisto0GenLeading(copy.fIJBckgHisto0GenLeading)
  ,fIJBckgHisto1RecCuts(copy.fIJBckgHisto1RecCuts)   
  ,fIJBckgHisto1RecLeading(copy.fIJBckgHisto1RecLeading)
  ,fIJBckgHisto1Gen(copy.fIJBckgHisto1Gen)       
  ,fIJBckgHisto1GenLeading(copy.fIJBckgHisto1GenLeading)
  ,fIJBckgHisto2RecCuts(copy.fIJBckgHisto2RecCuts)   
  ,fIJBckgHisto2RecLeading(copy.fIJBckgHisto2RecLeading)
  ,fIJBckgHisto2Gen(copy.fIJBckgHisto2Gen)       
  ,fIJBckgHisto2GenLeading(copy.fIJBckgHisto2GenLeading)
  ,fIJBckgHisto3RecCuts(copy.fIJBckgHisto3RecCuts)   
  ,fIJBckgHisto3RecLeading(copy.fIJBckgHisto3RecLeading)
  ,fIJBckgHisto3Gen(copy.fIJBckgHisto3Gen)       
  ,fIJBckgHisto3GenLeading(copy.fIJBckgHisto3GenLeading)
  ,fRandom(copy.fRandom)
  ,fBckgSubMethod(copy.fBckgSubMethod)
{
  // copy constructor
  fBckgType[0] = copy.fBckgType[0];
  fBckgType[1] = copy.fBckgType[1];
  fBckgType[2] = copy.fBckgType[2];
  fBckgType[3] = copy.fBckgType[3];
}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction& AliAnalysisTaskFragmentationFunction::operator=(const AliAnalysisTaskFragmentationFunction& o)
{
  // assignment
  
  if(this!=&o){

    AliAnalysisTaskSE::operator=(o);
    fESD                          = o.fESD;
    fAOD                          = o.fAOD;
    fBranchRecJets                = o.fBranchRecJets;
    fBranchRecBackJets            = o.fBranchRecBackJets;
    fBranchRecBckgClusters        = o.fBranchRecBckgClusters;
    fBranchGenJets                = o.fBranchGenJets;
    fTrackTypeGen                 = o.fTrackTypeGen;
    fJetTypeGen                   = o.fJetTypeGen;
    fJetTypeRecEff                = o.fJetTypeRecEff;
    fFilterMask                   = o.fFilterMask;
    fUsePhysicsSelection          = o.fUsePhysicsSelection;
    fEventClass                   = o.fEventClass;
    fMaxVertexZ                   = o.fMaxVertexZ;
    fTrackPtCut                   = o.fTrackPtCut;
    fTrackEtaMin                  = o.fTrackEtaMin;
    fTrackEtaMax                  = o.fTrackEtaMax;
    fTrackPhiMin                  = o.fTrackPhiMin;
    fTrackPhiMax                  = o.fTrackPhiMax;
    fJetPtCut                     = o.fJetPtCut;
    fJetEtaMin                    = o.fJetEtaMin;
    fJetEtaMax                    = o.fJetEtaMax;
    fJetPhiMin                    = o.fJetPhiMin;
    fJetPhiMax                    = o.fJetPhiMin;
    fDiJetCut                     = o.fDiJetCut;
    fDiJetDeltaPhiCut             = o.fDiJetDeltaPhiCut;
    fDiJetPtFractionCut           = o.fDiJetPtFractionCut;
    fDiJetCDFCut                  = o.fDiJetCDFCut;
    fDiJetKindBins                = o.fDiJetKindBins;
    fFFRadius                     = o.fFFRadius;
    fFFBckgRadius                 = o.fFFBckgRadius;
    fBckgMode                     = o.fBckgMode;
    fIJMode                       = o.fIJMode;
    fQAMode                       = o.fQAMode;
    fFFMode                       = o.fFFMode;
    fDJMode                       = o.fDJMode;
    fEffMode                      = o.fEffMode;
    fPhiCorrMode                  = o.fPhiCorrMode;
    fBckgType[0]                  = o.fBckgType[0];
    fBckgType[1]                  = o.fBckgType[1];
    fBckgType[2]                  = o.fBckgType[2];
    fBckgType[3]                  = o.fBckgType[3];
    fUseRecEffRecJetPtBins        = o.fUseRecEffRecJetPtBins;
    fUseResponseRecJetPtBins      = o.fUseResponseRecJetPtBins;
    fAvgTrials                    = o.fAvgTrials;
    fTracksRec                    = o.fTracksRec;
    fTracksRecCuts                = o.fTracksRecCuts;
    fTracksGen                    = o.fTracksGen;
    fTracksAODMCCharged           = o.fTracksAODMCCharged;
    fTracksRecQualityCuts         = o.fTracksRecQualityCuts;
    fJetsRec                      = o.fJetsRec;
    fJetsRecCuts                  = o.fJetsRecCuts;
    fJetsGen                      = o.fJetsGen;
    fJetsRecEff                   = o.fJetsRecEff;
    fBckgJetsRec                  = o.fBckgJetsRec;
    fBckgJetsRecCuts              = o.fBckgJetsRecCuts;
    fBckgJetsGen                  = o.fBckgJetsGen;
    fQATrackHistosRec             = o.fQATrackHistosRec;
    fQATrackHistosRecCuts         = o.fQATrackHistosRecCuts;
    fQATrackHistosGen             = o.fQATrackHistosGen;
    fQAJetHistosRec               = o.fQAJetHistosRec;
    fQAJetHistosRecCuts           = o.fQAJetHistosRecCuts;
    fQAJetHistosRecCutsLeading    = o.fQAJetHistosRecCutsLeading;
    fQAJetHistosGen               = o.fQAJetHistosGen;
    fQAJetHistosGenLeading        = o.fQAJetHistosGenLeading;
    fQAJetHistosRecEffLeading     = o.fQAJetHistosRecEffLeading;
    fFFHistosRecCuts              = o.fFFHistosRecCuts;
    fFFHistosRecLeading           = o.fFFHistosRecLeading;
    fFFHistosRecLeadingTrack      = o.fFFHistosRecLeadingTrack;
    fFFHistosGen                  = o.fFFHistosGen;
    fFFHistosGenLeading           = o.fFFHistosGenLeading;
    fFFHistosGenLeadingTrack      = o.fFFHistosGenLeadingTrack;
    fIJHistosRecCuts              = o.fIJHistosRecCuts;
    fIJHistosRecLeading           = o.fIJHistosRecLeading;
    fIJHistosRecLeadingTrack      = o.fIJHistosRecLeadingTrack;
    fIJHistosGen                  = o.fIJHistosGen;
    fIJHistosGenLeading           = o.fIJHistosGenLeading;
    fIJHistosGenLeadingTrack      = o.fIJHistosGenLeadingTrack;
    fFFDiJetHistosRecCuts         = o.fFFDiJetHistosRecCuts;
    fFFDiJetHistosRecLeading      = o.fFFDiJetHistosRecLeading;
    fFFDiJetHistosRecLeadingTrack = o.fFFDiJetHistosRecLeadingTrack;
    fFFDiJetHistosGen             = o.fFFDiJetHistosGen;
    fFFDiJetHistosGenLeading      = o.fFFDiJetHistosGenLeading;
    fFFDiJetHistosGenLeadingTrack = o.fFFDiJetHistosGenLeadingTrack;
    fQADiJetHistosRecCuts         = o.fQADiJetHistosRecCuts;
    fQADiJetHistosGen             = o.fQADiJetHistosGen;
    fPhiCorrHistosJetArea         = o.fPhiCorrHistosJetArea;
    fPhiCorrHistosTransverseArea  = o.fPhiCorrHistosTransverseArea;
    fPhiCorrHistosAwayArea        = o.fPhiCorrHistosAwayArea;
    fQATrackHighPtThreshold       = o.fQATrackHighPtThreshold; 
    fFFNBinsJetPt                 = o.fFFNBinsJetPt;    
    fFFJetPtMin                   = o.fFFJetPtMin; 
    fFFJetPtMax                   = o.fFFJetPtMax;
    fFFNBinsPt                    = o.fFFNBinsPt;      
    fFFPtMin                      = o.fFFPtMin;        
    fFFPtMax                      = o.fFFPtMax;        
    fFFNBinsXi                    = o.fFFNBinsXi;      
    fFFXiMin                      = o.fFFXiMin;        
    fFFXiMax                      = o.fFFXiMax;        
    fFFNBinsZ                     = o.fFFNBinsZ;       
    fFFZMin                       = o.fFFZMin;         
    fFFZMax                       = o.fFFZMax;         
    fQAJetNBinsPt                 = o.fQAJetNBinsPt;   
    fQAJetPtMin                   = o.fQAJetPtMin;     
    fQAJetPtMax                   = o.fQAJetPtMax;     
    fQAJetNBinsEta                = o.fQAJetNBinsEta;  
    fQAJetEtaMin                  = o.fQAJetEtaMin;    
    fQAJetEtaMax                  = o.fQAJetEtaMax;    
    fQAJetNBinsPhi                = o.fQAJetNBinsPhi;  
    fQAJetPhiMin                  = o.fQAJetPhiMin;    
    fQAJetPhiMax                  = o.fQAJetPhiMax;    
    fQATrackNBinsPt               = o.fQATrackNBinsPt; 
    fQATrackPtMin                 = o.fQATrackPtMin;   
    fQATrackPtMax                 = o.fQATrackPtMax;   
    fQATrackNBinsEta              = o.fQATrackNBinsEta;
    fQATrackEtaMin                = o.fQATrackEtaMin;  
    fQATrackEtaMax                = o.fQATrackEtaMax;  
    fQATrackNBinsPhi              = o.fQATrackNBinsPhi;
    fQATrackPhiMin                = o.fQATrackPhiMin;  
    fQATrackPhiMax                = o.fQATrackPhiMax;  
    fIJNBinsJetPt                 = o.fIJNBinsJetPt;
    fIJJetPtMin                   = o.fIJJetPtMin;
    fIJJetPtMax                   = o.fIJJetPtMax;
    fIJNBinsPt                    = o.fIJNBinsPt;
    fIJPtMin                      = o.fIJPtMin;
    fIJPtMax                      = o.fIJPtMax;
    fIJNBinsZ                     = o.fIJNBinsZ;
    fIJZMin                       = o.fIJZMin;
    fIJZMax                       = o.fIJZMax;
    fIJNBinsCosTheta              = o.fIJNBinsCosTheta;
    fIJCosThetaMin                = o.fIJCosThetaMin;
    fIJCosThetaMax                = o.fIJCosThetaMax;
    fIJNBinsTheta                 = o.fIJNBinsTheta;
    fIJThetaMin                   = o.fIJThetaMin;
    fIJThetaMax                   = o.fIJThetaMax;
    fIJNBinsJt                    = o.fIJNBinsJt;
    fIJJtMin                      = o.fIJJtMin;
    fIJJtMax                      = o.fIJJtMax;
    fDiJetNBinsJetInvMass         = o.fDiJetNBinsJetInvMass;
    fDiJetJetInvMassMin           = o.fDiJetJetInvMassMin;
    fDiJetJetInvMassMax           = o.fDiJetJetInvMassMax;
    fDiJetNBinsJetPt              = o.fDiJetNBinsJetPt;
    fDiJetJetPtMin                = o.fDiJetJetPtMin;
    fDiJetJetPtMax                = o.fDiJetJetPtMax;
    fDiJetNBinsPt                 = o.fDiJetNBinsPt;
    fDiJetPtMin                   = o.fDiJetPtMin;
    fDiJetPtMax                   = o.fDiJetPtMax;
    fDiJetNBinsXi                 = o.fDiJetNBinsXi;
    fDiJetXiMin                   = o.fDiJetXiMin;
    fDiJetXiMax                   = o.fDiJetXiMax;
    fDiJetNBinsZ                  = o.fDiJetNBinsZ;
    fDiJetZMin                    = o.fDiJetZMin;
    fDiJetZMax                    = o.fDiJetZMax;
    fQADiJetNBinsInvMass          = o.fQADiJetNBinsInvMass;
    fQADiJetInvMassMin            = o.fQADiJetInvMassMin;
    fQADiJetInvMassMax            = o.fQADiJetInvMassMax;
    fQADiJetNBinsJetPt            = o.fQADiJetNBinsJetPt;
    fQADiJetJetPtMin              = o.fQADiJetJetPtMin;
    fQADiJetJetPtMax              = o.fQADiJetJetPtMax;
    fQADiJetNBinsDeltaPhi         = o.fQADiJetNBinsDeltaPhi;
    fQADiJetDeltaPhiMin           = o.fQADiJetDeltaPhiMin;
    fQADiJetDeltaPhiMax           = o.fQADiJetDeltaPhiMax;
    fQADiJetNBinsDeltaEta         = o.fQADiJetNBinsDeltaEta;
    fQADiJetDeltaEtaMin           = o.fQADiJetDeltaEtaMin;
    fQADiJetDeltaEtaMax           = o.fQADiJetDeltaEtaMax;
    fQADiJetNBinsDeltaPt          = o.fQADiJetNBinsDeltaPt;
    fQADiJetDeltaPtMin            = o.fQADiJetDeltaPtMin;
    fQADiJetDeltaPtMax            = o.fQADiJetDeltaPtMax;
    fQADiJetNBinsInBal            = o.fQADiJetNBinsInBal;
    fQADiJetInBalMin              = o.fQADiJetInBalMin;
    fQADiJetInBalMax              = o.fQADiJetInBalMax;
    fPhiCorrNBinsPt               = o.fPhiCorrNBinsPt;
    fPhiCorrPtMin                 = o.fPhiCorrPtMin;
    fPhiCorrPtMax                 = o.fPhiCorrPtMax;
    fPhiCorrNBinsEta              = o.fPhiCorrNBinsEta;
    fPhiCorrEtaMin                = o.fPhiCorrEtaMin;
    fPhiCorrEtaMax                = o.fPhiCorrEtaMax;
    fPhiCorrNBinsPhi              = o.fPhiCorrNBinsPhi;
    fPhiCorrPhiMin                = o.fPhiCorrPhiMin;
    fPhiCorrPhiMax                = o.fPhiCorrPhiMax;
    fCommonHistList               = o.fCommonHistList;
    fh1EvtSelection               = o.fh1EvtSelection;
    fh1VertexNContributors        = o.fh1VertexNContributors;
    fh1VertexZ                    = o.fh1VertexZ;
    fh1EvtMult                    = o.fh1EvtMult;
    fh1EvtCent                    = o.fh1EvtCent;
    fh1Xsec                       = o.fh1Xsec;
    fh1Trials                     = o.fh1Trials;
    fh1PtHard                     = o.fh1PtHard;
    fh1PtHardTrials               = o.fh1PtHardTrials;
    fh1nRecJetsCuts               = o.fh1nRecJetsCuts;
    fh1nGenJets                   = o.fh1nGenJets; 
    fh1nRecEffJets                = o.fh1nRecEffJets;
    fh2PtRecVsGenPrim             = o.fh2PtRecVsGenPrim;
    fQATrackHistosRecEffGen       = o.fQATrackHistosRecEffGen;  
    fQATrackHistosRecEffRec       = o.fQATrackHistosRecEffRec;  
    fFFHistosRecEffGen            = o.fFFHistosRecEffGen;    
    fFFHistosRecEffRec            = o.fFFHistosRecEffRec;   
    fhnResponseSinglePt           = o.fhnResponseSinglePt;
    fhnResponseJetTrackPt         = o.fhnResponseJetTrackPt;
    fhnResponseJetZ               = o.fhnResponseJetZ;
    fhnResponseJetXi              = o.fhnResponseJetXi;
    // Background
    fh1OutLeadingMult             = o.fh1OutLeadingMult;
    fh1OutLeadingStatMult         = o.fh1OutLeadingStatMult;
    fh1PerpMult                   = o.fh1PerpMult;
    fh1ASideMult                  = o.fh1ASideMult;
    fh1ASideWindowMult            = o.fh1ASideWindowMult;
    fh1PerpWindowMult             = o.fh1PerpWindowMult;
    fh1Out2JetsMult               = o.fh1Out2JetsMult;
    fh1Out3JetsMult               = o.fh1Out3JetsMult;
    fh1MedianClustersMult         = o.fh1MedianClustersMult;
    fh1OutClustersMult            = o.fh1OutClustersMult;
    fQABckgHisto0RecCuts          = o.fQABckgHisto0RecCuts;  
    fQABckgHisto0Gen              = o.fQABckgHisto0Gen;      
    fQABckgHisto1RecCuts          = o.fQABckgHisto1RecCuts;  
    fQABckgHisto1Gen              = o.fQABckgHisto1Gen;      
    fQABckgHisto2RecCuts          = o.fQABckgHisto2RecCuts;  
    fQABckgHisto2Gen              = o.fQABckgHisto2Gen;  
    fQABckgHisto3RecCuts          = o.fQABckgHisto3RecCuts;  
    fQABckgHisto3Gen              = o.fQABckgHisto3Gen;  
    fFFBckgHisto0RecCuts          = o.fFFBckgHisto0RecCuts;
    fFFBckgHisto0RecLeading       = o.fFFBckgHisto0RecLeading;
    fFFBckgHisto0Gen              = o.fFFBckgHisto0Gen;       
    fFFBckgHisto0GenLeading       = o.fFFBckgHisto0GenLeading;
    fFFBckgHisto1RecCuts          = o.fFFBckgHisto1RecCuts;
    fFFBckgHisto1RecLeading       = o.fFFBckgHisto1RecLeading;
    fFFBckgHisto1Gen              = o.fFFBckgHisto1Gen;       
    fFFBckgHisto1GenLeading       = o.fFFBckgHisto1GenLeading;
    fFFBckgHisto2RecCuts          = o.fFFBckgHisto2RecCuts;
    fFFBckgHisto2RecLeading       = o.fFFBckgHisto2RecLeading;
    fFFBckgHisto2Gen              = o.fFFBckgHisto2Gen;       
    fFFBckgHisto2GenLeading       = o.fFFBckgHisto2GenLeading;
    fFFBckgHisto3RecCuts          = o.fFFBckgHisto3RecCuts;
    fFFBckgHisto3RecLeading       = o.fFFBckgHisto3RecLeading;
    fFFBckgHisto3Gen              = o.fFFBckgHisto3Gen;       
    fFFBckgHisto3GenLeading       = o.fFFBckgHisto3GenLeading;
    fIJBckgHisto0RecCuts          = o.fIJBckgHisto0RecCuts;   
    fIJBckgHisto0RecLeading       = o.fIJBckgHisto0RecLeading;
    fIJBckgHisto0Gen              = o.fIJBckgHisto0Gen;       
    fIJBckgHisto0GenLeading       = o.fIJBckgHisto0GenLeading;
    fIJBckgHisto1RecCuts          = o.fIJBckgHisto1RecCuts;   
    fIJBckgHisto1RecLeading       = o.fIJBckgHisto1RecLeading;
    fIJBckgHisto1Gen              = o.fIJBckgHisto1Gen;       
    fIJBckgHisto1GenLeading       = o.fIJBckgHisto1GenLeading;
    fIJBckgHisto2RecCuts          = o.fIJBckgHisto2RecCuts;   
    fIJBckgHisto2RecLeading       = o.fIJBckgHisto2RecLeading;
    fIJBckgHisto2Gen              = o.fIJBckgHisto2Gen;       
    fIJBckgHisto2GenLeading       = o.fIJBckgHisto2GenLeading;
    fIJBckgHisto3Gen              = o.fIJBckgHisto3Gen;       
    fIJBckgHisto3GenLeading       = o.fIJBckgHisto3GenLeading;
    fRandom                       = o.fRandom;
    fBckgSubMethod                = o.fBckgSubMethod;
  }
    
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::~AliAnalysisTaskFragmentationFunction()
{
  // destructor
  
  if(fTracksRec)            delete fTracksRec;
  if(fTracksRecCuts)        delete fTracksRecCuts;
  if(fTracksGen)            delete fTracksGen;
  if(fTracksAODMCCharged)   delete fTracksAODMCCharged;  
  if(fTracksRecQualityCuts) delete fTracksRecQualityCuts; 
  if(fJetsRec)              delete fJetsRec;
  if(fJetsRecCuts)          delete fJetsRecCuts;
  if(fJetsGen)              delete fJetsGen;
  if(fJetsRecEff)           delete fJetsRecEff;
  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || 
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
      fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)){

    if(fBckgJetsRec)          delete fBckgJetsRec;
    if(fBckgJetsRecCuts)      delete fBckgJetsRecCuts;
    if(fBckgJetsGen)          delete fBckgJetsGen;
  }
  if(fRandom)               delete fRandom;
}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AliFragFuncHistos(const char* name, 
							 Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
							 Int_t nPt, Float_t ptMin, Float_t ptMax,
							 Int_t nXi, Float_t xiMin, Float_t xiMax,
							 Int_t nZ , Float_t zMin , Float_t zMax )
  : TObject()
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsXi(nXi) 
  ,fXiMin(xiMin)   
  ,fXiMax(xiMax)   
  ,fNBinsZ(nZ)  
  ,fZMin(zMin)    
  ,fZMax(zMax)    
  ,fh2TrackPt(0)
  ,fh2Xi(0)
  ,fh2Z(0)
  ,fh1JetPt(0)
  ,fNameFF(name)
{
  // default constructor

}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AliFragFuncHistos(const AliFragFuncHistos& copy)
  : TObject()
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   
  ,fPtMax(copy.fPtMax)   
  ,fNBinsXi(copy.fNBinsXi) 
  ,fXiMin(copy.fXiMin)   
  ,fXiMax(copy.fXiMax)   
  ,fNBinsZ(copy.fNBinsZ)  
  ,fZMin(copy.fZMin)    
  ,fZMax(copy.fZMax)    
  ,fh2TrackPt(copy.fh2TrackPt)
  ,fh2Xi(copy.fh2Xi)
  ,fh2Z(copy.fh2Z)
  ,fh1JetPt(copy.fh1JetPt)
  ,fNameFF(copy.fNameFF)
{
  // copy constructor
}

//_______________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsJetPt = o.fNBinsJetPt;
    fJetPtMin   = o.fJetPtMin;
    fJetPtMax   = o.fJetPtMax;
    fNBinsPt    = o.fNBinsPt; 
    fPtMin      = o.fPtMin;   
    fPtMax      = o.fPtMax;   
    fNBinsXi    = o.fNBinsXi; 
    fXiMin      = o.fXiMin;   
    fXiMax      = o.fXiMax;   
    fNBinsZ     = o.fNBinsZ;  
    fZMin       = o.fZMin;    
    fZMax       = o.fZMax;    
    fh2TrackPt  = o.fh2TrackPt;
    fh2Xi       = o.fh2Xi;
    fh2Z        = o.fh2Z;
    fh1JetPt    = o.fh1JetPt;
    fNameFF     = o.fNameFF;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::~AliFragFuncHistos()
{
  // destructor 

  if(fh1JetPt)   delete fh1JetPt;
  if(fh2TrackPt) delete fh2TrackPt;
  if(fh2Xi)      delete fh2Xi;
  if(fh2Z)       delete fh2Z;
}

//_________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::DefineHistos()
{
  // book FF histos

  fh1JetPt   = new TH1F(Form("fh1FFJetPt%s", fNameFF.Data()),"",fNBinsJetPt,fJetPtMin,fJetPtMax);
  fh2TrackPt = new TH2F(Form("fh2FFTrackPt%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsPt, fPtMin, fPtMax);
  fh2Xi      = new TH2F(Form("fh2FFXi%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsXi, fXiMin, fXiMax);
  fh2Z       = new TH2F(Form("fh2FFZ%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsZ, fZMin, fZMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{T} [GeV/c]", "entries"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPt,"jet p_{T} [GeV/c]","p_{T} [GeV/c]","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi,"jet p_{T} [GeV/c]","#xi", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z,"jet p_{T} [GeV/c]","z","entries");
}

//_______________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::FillFF(Float_t trackPt, Float_t jetPt, Bool_t incrementJetPt, Float_t norm)
{
  // fill FF
 
  if(incrementJetPt && norm) fh1JetPt->Fill(jetPt,1/norm);
  else if(incrementJetPt) fh1JetPt->Fill(jetPt);

 // Added for proper normalization of FF background estimation
  // when zero track are found in the background region
  if(trackPt==-1.) return;
 
  if(norm)fh2TrackPt->Fill(jetPt,trackPt,1/norm);
  else fh2TrackPt->Fill(jetPt,trackPt);
  
  Double_t z = 0.;
  if(jetPt>0) z = trackPt / jetPt;
  Double_t xi = 0;
  if(z>0) xi = TMath::Log(1/z);
  
  if(norm){
    fh2Xi->Fill(jetPt,xi,1/norm);
    fh2Z->Fill(jetPt,z,1/norm);
  }
  else {
    fh2Xi->Fill(jetPt,xi);
    fh2Z->Fill(jetPt,z);
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  
  list->Add(fh2TrackPt);
  list->Add(fh2Xi);
  list->Add(fh2Z);
}

//_________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const char* name,
							       Int_t nPt,   Float_t ptMin,   Float_t ptMax,
							       Int_t nEta,  Float_t etaMin,  Float_t etaMax,
							       Int_t nPhi,  Float_t phiMin,  Float_t phiMax)
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsEta(nEta)
  ,fEtaMin(etaMin)
  ,fEtaMax(etaMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fh2EtaPhi(0)
  ,fh1Pt(0)
  ,fNameQAJ(name)
{
  // default constructor
}

//____________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const AliFragFuncQAJetHistos& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsEta(copy.fNBinsEta)
  ,fEtaMin(copy.fEtaMin)
  ,fEtaMax(copy.fEtaMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fh2EtaPhi(copy.fh2EtaPhi)
  ,fh1Pt(copy.fh1Pt)
  ,fNameQAJ(copy.fNameQAJ)
{
  // copy constructor
}

//________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt  = o.fNBinsPt;
    fPtMin    = o.fPtMin;
    fPtMax    = o.fPtMax;
    fNBinsEta = o.fNBinsEta;
    fEtaMin   = o.fEtaMin;
    fEtaMax   = o.fEtaMax;
    fNBinsPhi = o.fNBinsPhi;
    fPhiMin   = o.fPhiMin;
    fPhiMax   = o.fPhiMax;
    fh2EtaPhi = o.fh2EtaPhi;
    fh1Pt     = o.fh1Pt;
    fNameQAJ  = o.fNameQAJ;
  }
  
  return *this;
}

//______________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::~AliFragFuncQAJetHistos()
{
  // destructor 
  
  if(fh2EtaPhi) delete fh2EtaPhi;
  if(fh1Pt)     delete fh1Pt;
}

//____________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::DefineHistos()
{
  // book jet QA histos

  fh2EtaPhi  = new TH2F(Form("fh2JetQAEtaPhi%s", fNameQAJ.Data()), Form("%s: #eta - #phi distribution", fNameQAJ.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt      = new TH1F(Form("fh1JetQAPt%s", fNameQAJ.Data()), Form("%s: p_{T} distribution", fNameQAJ.Data()), fNBinsPt, fPtMin, fPtMax);
	
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
}

//____________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::FillJetQA(Float_t eta, Float_t phi, Float_t pt)
{
  // fill jet QA histos 

  fh2EtaPhi->Fill( eta, phi);
  fh1Pt->Fill( pt );
}

//____________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQAJetHistos::AddToOutput(TList* list) const 
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh1Pt);
}

//___________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const char* name,
								   Int_t nPt, Float_t ptMin, Float_t ptMax,
								   Int_t nEta, Float_t etaMin, Float_t etaMax,
								   Int_t nPhi, Float_t phiMin, Float_t phiMax,
								   Float_t ptThresh) 
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsEta(nEta)
  ,fEtaMin(etaMin)
  ,fEtaMax(etaMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fHighPtThreshold(ptThresh)
  ,fh2EtaPhi(0)
  ,fh1Pt(0)
  ,fh2HighPtEtaPhi(0)
  ,fh2PhiPt(0)
  ,fNameQAT(name)
{
  // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const AliFragFuncQATrackHistos& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsEta(copy.fNBinsEta)
  ,fEtaMin(copy.fEtaMin)
  ,fEtaMax(copy.fEtaMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fHighPtThreshold(copy.fHighPtThreshold)
  ,fh2EtaPhi(copy.fh2EtaPhi)
  ,fh1Pt(copy.fh1Pt)
  ,fh2HighPtEtaPhi(copy.fh2HighPtEtaPhi)
  ,fh2PhiPt(copy.fh2PhiPt)
  ,fNameQAT(copy.fNameQAT)
{
  // copy constructor
}

// _____________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt         = o.fNBinsPt;
    fPtMin           = o.fPtMin;
    fPtMax           = o.fPtMax;
    fNBinsEta        = o.fNBinsEta;
    fEtaMin          = o.fEtaMin;
    fEtaMax          = o.fEtaMax;
    fNBinsPhi        = o.fNBinsPhi;
    fPhiMin          = o.fPhiMin;
    fPhiMax          = o.fPhiMax;
    fHighPtThreshold = o.fHighPtThreshold;
    fh2EtaPhi        = o.fh2EtaPhi;
    fh1Pt            = o.fh1Pt;
    fh2HighPtEtaPhi  = o.fh2HighPtEtaPhi;
    fh2PhiPt         = o.fh2PhiPt;
    fNameQAT         = o.fNameQAT;
  }
  
  return *this;
}

//___________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::~AliFragFuncQATrackHistos()
{
  // destructor 
  
  if(fh2EtaPhi)       delete fh2EtaPhi;
  if(fh2HighPtEtaPhi) delete fh2HighPtEtaPhi;
  if(fh1Pt)           delete fh1Pt;
  if(fh2PhiPt)        delete fh2PhiPt;
}

//______________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::DefineHistos()
{
  // book track QA histos

  fh2EtaPhi       = new TH2F(Form("fh2TrackQAEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh2HighPtEtaPhi = new TH2F(Form("fh2TrackQAHighPtEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution for high-p_{T}", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt           = new TH1F(Form("fh1TrackQAPt%s", fNameQAT.Data()), Form("%s: p_{T} distribution", fNameQAT.Data()), fNBinsPt, fPtMin, fPtMax);
    fh2PhiPt        = new TH2F(Form("fh2TrackQAPhiPt%s", fNameQAT.Data()), Form("%s: #phi - #p_{T} distribution", fNameQAT.Data()), fNBinsPhi, fPhiMin, fPhiMax, fNBinsPt, fPtMin, fPtMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2HighPtEtaPhi, "#eta", "#phi");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2PhiPt, "#phi", "p_{T} [GeV/c]"); 
}

//________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::FillTrackQA(Float_t eta, Float_t phi, Float_t pt, Bool_t weightPt, Float_t norm)
{
  // fill track QA histos
  Float_t weight = 1.;
  if(weightPt) weight = pt;  
  fh2EtaPhi->Fill( eta, phi, weight);
  if(pt > fHighPtThreshold) fh2HighPtEtaPhi->Fill( eta, phi, weight);
  if(norm) fh1Pt->Fill( pt, 1/norm );
  else fh1Pt->Fill( pt );
  fh2PhiPt->Fill(phi, pt);
}

//______________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh2HighPtEtaPhi);
  list->Add(fh1Pt);
  list->Add(fh2PhiPt);
}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::AliFragFuncIntraJetHistos(const char* name, 
							 Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
							 Int_t nPt, Float_t ptMin, Float_t ptMax,
							 Int_t nZ , Float_t zMin , Float_t zMax,
							 Int_t nCosTheta , Float_t costhetaMin , Float_t costhetaMax,
							 Int_t nTheta , Float_t thetaMin , Float_t thetaMax,
							 Int_t nJt , Float_t jtMin , Float_t jtMax)
  : TObject()
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsZ(nZ) 
  ,fZMin(zMin)   
  ,fZMax(zMax)   
  ,fNBinsJt(nJt)
  ,fJtMin(jtMin)
  ,fJtMax(jtMax)
  ,fNBinsTheta(nTheta)
  ,fThetaMin(thetaMin)
  ,fThetaMax(thetaMax)
  ,fNBinsCosTheta(nCosTheta)
  ,fCosThetaMin(costhetaMin)
  ,fCosThetaMax(costhetaMax)
  ,fh2CosTheta(0)
  ,fh2PtZ(0)
  ,fh3ThetaZ(0)
  ,fh3JtTheta(0)
  ,fh3JtZ(0)
  ,fNameIJ(name)
{
  // default constructor

}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::AliFragFuncIntraJetHistos(const AliFragFuncIntraJetHistos& copy)
  : TObject()
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   
  ,fPtMax(copy.fPtMax)   
  ,fNBinsZ(copy.fNBinsZ) 
  ,fZMin(copy.fZMin)   
  ,fZMax(copy.fZMax)   
  ,fNBinsJt(copy.fNBinsJt)
  ,fJtMin(copy.fJtMin)
  ,fJtMax(copy.fJtMax)
  ,fNBinsTheta(copy.fNBinsTheta)
  ,fThetaMin(copy.fThetaMin)
  ,fThetaMax(copy.fThetaMax)
  ,fNBinsCosTheta(copy.fNBinsCosTheta)
  ,fCosThetaMin(copy.fCosThetaMin)
  ,fCosThetaMax(copy.fCosThetaMax)
  ,fh2CosTheta(copy.fh2CosTheta)
  ,fh2PtZ(copy.fh2PtZ)
  ,fh3ThetaZ(copy.fh3ThetaZ)
  ,fh3JtTheta(copy.fh3JtTheta)
  ,fh3JtZ(copy.fh3JtZ)
  ,fNameIJ(copy.fNameIJ)
{
  // copy constructor
}

//_______________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsJetPt       = o.fNBinsJetPt;
    fJetPtMin         = o.fJetPtMin;
    fJetPtMax         = o.fJetPtMax;
    fNBinsPt          = o.fNBinsPt; 
    fPtMin            = o.fPtMin;   
    fPtMax            = o.fPtMax;   
    fNBinsZ           = o.fNBinsZ; 
    fZMin             = o.fZMin;   
    fZMax             = o.fZMax;   
    fNBinsJt          = o.fNBinsJt;
    fJtMin            = o.fJtMin;
    fJtMax            = o.fJtMax;
    fNBinsTheta       = o.fNBinsTheta;
    fThetaMin         = o.fThetaMin;
    fThetaMax         = o.fThetaMax;
    fNBinsCosTheta    = o.fNBinsCosTheta;
    fCosThetaMin      = o.fCosThetaMin;
    fCosThetaMax      = o.fCosThetaMax;
    fh2CosTheta       = o.fh2CosTheta;
    fh2PtZ            = o.fh2PtZ;
    fh3ThetaZ         = o.fh3ThetaZ;
    fh3JtTheta        = o.fh3JtTheta;
    fh3JtZ            = o.fh3JtZ;
    fNameIJ           = o.fNameIJ;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::~AliFragFuncIntraJetHistos()
{
  // destructor 


  if(fh2CosTheta)       delete fh2CosTheta;
  if(fh2PtZ)            delete fh2PtZ;
  if(fh3ThetaZ)         delete fh3ThetaZ;             
  if(fh3JtTheta)        delete fh3JtTheta;
  if(fh3JtZ)            delete fh3JtZ;

}

//_________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::DefineHistos()
{
  // book FF histos

  fh2CosTheta = new TH2F(Form("fh2IJcosTheta%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsCosTheta, fCosThetaMin, fCosThetaMax);
  fh2PtZ      = new TH2F(Form("fh2IJPtZ%s",fNameIJ.Data()),"",fNBinsPt, fPtMin, fPtMax, fNBinsZ, fZMin, fZMax);
  fh3ThetaZ   = new TH3F(Form("fh3IJThetaZ%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsTheta, fThetaMin, fThetaMax, fNBinsZ, fZMin, fZMax);
  fh3JtTheta  = new TH3F(Form("fh3IJJtTheta%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsJt, fJtMin, fJtMax, fNBinsTheta, fThetaMin, fThetaMax);
  fh3JtZ      = new TH3F(Form("fh3IJJtZ%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsJt, fJtMin, fJtMax, fNBinsZ, fZMin, fZMax);
  
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2CosTheta,"jet p_{T} [GeV/c]","cos(#Theta)", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2PtZ,"had p_{T} [GeV/c]","z=p_{T}^{had}/p_{T}^{jet}","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh3ThetaZ,"jet p_{T} [GeV/c]","#Theta [rad]","z=p_{T}^{had}/p_{T}^{jet}");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh3JtTheta,"jet p_{T} [GeV/c]","j_{T} [GeV/c]","#Theta [rad]");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh3JtZ,"jet p_{T} [GeV/c]","j_{T} [GeV/c]","z=p_{T}^{had}/p_{T}^{jet}");

}

//_______________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::FillIntraJet(const TLorentzVector* trackV, const TLorentzVector* jetV, Float_t norm)
{
  // fill IntraJet histos
 
  Float_t cosTheta = 0.; Float_t theta = 0.; 
  Float_t jt = 0.; Float_t z = 0.; 
  // For Theta distribution
  Float_t pxT  = trackV->Px();
  Float_t pyT  = trackV->Py();
  Float_t pzT  = trackV->Pz();
  Float_t ptT  = trackV->Pt();
  Float_t pT   = trackV->P();
  Float_t etaT = trackV->Eta();
  Float_t phiT = trackV->Phi(); // Check the value returned
  Float_t pxJ = jetV->Px();
  Float_t pyJ = jetV->Py();
  Float_t pzJ = jetV->Pz();
  Float_t ptJ = jetV->Pt();
  Float_t pJ  = jetV->P();

  // Compute z
  if(ptJ>0) z = (Float_t)(ptT/ptJ);

  // Compute theta
  cosTheta = (pxT*pxJ+pyT*pyJ+pzT*pzJ)/(pT*pJ);
  theta = TMath::ACos(cosTheta);

  // Compute jt
  TVector3 trackP; TVector3 jetP;
  jetP[0] = pxJ;
  jetP[1] = pyJ;
  jetP[2] = pzJ;
  trackP.SetPtEtaPhi(ptT,etaT,phiT);
  jt = TMath::Sin(trackP.Angle(jetP))*trackP.Mag();

  // Fill histos
  if(norm){
    fh2CosTheta->Fill(ptJ,cosTheta,1/norm);
    fh2PtZ->Fill(ptT,z,1/norm);
    fh3ThetaZ->Fill(ptJ,theta,z,1/norm);
    fh3JtTheta->Fill(ptJ,jt,theta,1/norm);
    fh3JtZ->Fill(ptJ,jt,z,1/norm);
  }
  else {
    fh2CosTheta->Fill(ptJ,cosTheta);
    fh2PtZ->Fill(ptT,z);
    fh3ThetaZ->Fill(ptJ,theta,z);
    fh3JtTheta->Fill(ptJ,jt,theta);
    fh3JtZ->Fill(ptJ,jt,z);
  }

}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::AliFragFuncDiJetHistos(const char* name, Int_t kindSlices,
							 Int_t nJetInvMass, Float_t jetInvMassMin, Float_t jetInvMassMax,  
							 Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
							 Int_t nPt, Float_t ptMin, Float_t ptMax,
							 Int_t nXi, Float_t xiMin, Float_t xiMax,
							 Int_t nZ , Float_t zMin , Float_t zMax)
  : TObject()
  ,fKindSlices(kindSlices)
  ,fNBinsJetInvMass(nJetInvMass)
  ,fJetInvMassMin(jetInvMassMin)
  ,fJetInvMassMax(jetInvMassMax)
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsXi(nXi) 
  ,fXiMin(xiMin)   
  ,fXiMax(xiMax)   
  ,fNBinsZ(nZ)  
  ,fZMin(zMin)    
  ,fZMax(zMax)
  ,fh2TrackPtJet1(0)
  ,fh2TrackPtJet2(0)
  ,fh2TrackPtJet(0)
  ,fh1Jet1Pt(0)
  ,fh1Jet2Pt(0)
  ,fh1JetPt(0)
  ,fh2Xi1(0)
  ,fh2Xi2(0)
  ,fh2Xi(0)
  ,fh2Z1(0)
  ,fh2Z2(0)
  ,fh2Z(0)
  ,fh2Pt1(0)
  ,fh2Pt2(0)
  ,fh2Pt(0)
  ,fNameDJ(name)
{
  // default constructor

}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::AliFragFuncDiJetHistos(const AliFragFuncDiJetHistos& copy)
  : TObject()
  ,fKindSlices(copy.fKindSlices)
  ,fNBinsJetInvMass(copy.fNBinsJetInvMass)
  ,fJetInvMassMin(copy.fJetInvMassMin)
  ,fJetInvMassMax(copy.fJetInvMassMax)
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   
  ,fPtMax(copy.fPtMax)   
  ,fNBinsXi(copy.fNBinsXi) 
  ,fXiMin(copy.fXiMin)   
  ,fXiMax(copy.fXiMax)   
  ,fNBinsZ(copy.fNBinsZ)  
  ,fZMin(copy.fZMin)    
  ,fZMax(copy.fZMax)
  ,fh2TrackPtJet1(copy.fh2TrackPtJet1)
  ,fh2TrackPtJet2(copy.fh2TrackPtJet2)
  ,fh2TrackPtJet(copy.fh2TrackPtJet)
  ,fh1Jet1Pt(copy.fh1Jet1Pt)
  ,fh1Jet2Pt(copy.fh1Jet2Pt)
  ,fh1JetPt(copy.fh1JetPt)
  ,fh2Xi1(copy.fh2Xi1)
  ,fh2Xi2(copy.fh2Xi2)
  ,fh2Xi(copy.fh2Xi2)
  ,fh2Z1(copy.fh2Z1)
  ,fh2Z2(copy.fh2Z2)
  ,fh2Z(copy.fh2Z)
  ,fh2Pt1(copy.fh2Pt1)
  ,fh2Pt2(copy.fh2Pt2)
  ,fh2Pt(copy.fh2Pt)
  ,fNameDJ(copy.fNameDJ)
{
  // default constructor

}

//_______________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fKindSlices      = o.fKindSlices;
    fNBinsJetInvMass = o.fNBinsJetInvMass;
    fJetInvMassMin   = o.fJetInvMassMin;
    fJetInvMassMax   = o.fJetInvMassMax;
    fNBinsJetPt      = o.fNBinsJetPt;
    fJetPtMin        = o.fJetPtMin;
    fJetPtMax        = o.fJetPtMax;
    fNBinsPt         = o.fNBinsPt; 
    fPtMin           = o.fPtMin;   
    fPtMax           = o.fPtMax;   
    fNBinsXi         = o.fNBinsXi; 
    fXiMin           = o.fXiMin;   
    fXiMax           = o.fXiMax;   
    fNBinsZ          = o.fNBinsZ;  
    fZMin            = o.fZMin;    
    fZMax            = o.fZMax;   
    fh2TrackPtJet1   = o.fh2TrackPtJet1;
    fh2TrackPtJet2   = o.fh2TrackPtJet2;
    fh2TrackPtJet    = o.fh2TrackPtJet;
    fh1Jet1Pt        = o.fh1Jet1Pt;
    fh1Jet2Pt        = o.fh1Jet2Pt;
    fh1JetPt         = o.fh1JetPt;
    fh2Xi1           = o.fh2Xi1;
    fh2Xi2           = o.fh2Xi2;
    fh2Xi            = o.fh2Xi;
    fh2Z1            = o.fh2Z1;
    fh2Z2            = o.fh2Z2;
    fh2Z             = o.fh2Z;
    fh2Pt1           = o.fh2Pt1;
    fh2Pt2           = o.fh2Pt2;
    fh2Pt            = o.fh2Pt;
    fNameDJ          = o.fNameDJ;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::~AliFragFuncDiJetHistos()
{
  // destructor 

  if(fh2TrackPtJet1) delete fh2TrackPtJet1;
  if(fh2TrackPtJet2) delete fh2TrackPtJet2;
  if(fh2TrackPtJet ) delete fh2TrackPtJet;
  if(fh1Jet1Pt)      delete fh1Jet1Pt;
  if(fh1Jet2Pt)      delete fh1Jet2Pt;
  if(fh1JetPt)       delete fh1JetPt;
  if(fh2Xi1)         delete fh2Xi1;
  if(fh2Xi2)         delete fh2Xi2;
  if(fh2Xi)          delete fh2Xi;
  if(fh2Z1)          delete fh2Z1;
  if(fh2Z2)          delete fh2Z2;
  if(fh2Z)           delete fh2Z;
  if(fh2Pt1)         delete fh2Pt1;
  if(fh2Pt2)         delete fh2Pt2;
  if(fh2Pt)          delete fh2Pt;
}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::DefineDiJetHistos()
{
  // book DiJet histos
  
  Int_t nBins = 0;
  Double_t min = 0.;
  Double_t max = 0.;
  const char *xaxis = "";
  if(fKindSlices == 1)
    {
      nBins = fNBinsJetInvMass;
      min   = fJetInvMassMin;
      max   = fJetInvMassMax;
      xaxis = "M_{JJ} [GeV]";
    }
  if(fKindSlices == 2 || fKindSlices == 3)
    {
      nBins = fNBinsJetPt;
      min   = fJetPtMin;
      max   = fJetPtMax;
      if(fKindSlices == 2) xaxis = "E_{Tmean} [GeV]";
      if(fKindSlices == 3) xaxis ="leading jet p_{T} [GeV/c]";
    }
  
  fh1Jet1Pt      = new TH1F(Form("fh1DJJet1Pt%s", fNameDJ.Data()), "", fNBinsJetPt, fJetPtMin, fJetPtMax);
  fh1Jet2Pt      = new TH1F(Form("fh1DJJet2Pt%s", fNameDJ.Data()), "", fNBinsJetPt, fJetPtMin, fJetPtMax);
  fh1JetPt       = new TH1F(Form("fh1DJJetPt%s",  fNameDJ.Data()), "", fNBinsJetPt, fJetPtMin, fJetPtMax);
  
  fh2TrackPtJet1 = new TH2F(Form("fh2DJTrackPtJet1%s", fNameDJ.Data()), "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
  fh2TrackPtJet2 = new TH2F(Form("fh2DJTrackPtJet2%s", fNameDJ.Data()), "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
  fh2TrackPtJet  = new TH2F(Form("fh2DJTrackPtJet%s", fNameDJ.Data()),  "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
  
  fh2Xi1         = new TH2F(Form("fh2DJXi1%s", fNameDJ.Data()), "",nBins, min, max, fNBinsXi, fXiMin, fXiMax);
  fh2Xi2         = new TH2F(Form("fh2DJXi2%s", fNameDJ.Data()), "",nBins, min, max, fNBinsXi, fXiMin, fXiMax);
  fh2Xi          = new TH2F(Form("fh2DJXi%s", fNameDJ.Data()),  "",nBins, min, max, fNBinsXi, fXiMin, fXiMax);
  
  fh2Z1          = new TH2F(Form("fh2DJZ1%s", fNameDJ.Data()), "",nBins, min, max, fNBinsZ, fZMin, fZMax);
  fh2Z2          = new TH2F(Form("fh2DJZ2%s", fNameDJ.Data()), "",nBins, min, max, fNBinsZ, fZMin, fZMax);
  fh2Z           = new TH2F(Form("fh2DJZ%s", fNameDJ.Data()),  "",nBins, min, max, fNBinsZ, fZMin, fZMax);
  
  fh2Pt1         = new TH2F(Form("fh2DJPt1%s", fNameDJ.Data()), "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
  fh2Pt2         = new TH2F(Form("fh2DJPt2%s", fNameDJ.Data()), "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
  fh2Pt          = new TH2F(Form("fh2DJPtZ%s", fNameDJ.Data()),  "",nBins, min, max, fNBinsPt, fPtMin, fPtMax);
      
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Jet1Pt, "p_{T} [GeV/c]", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Jet2Pt, "p_{T} [GeV/c]", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{T} [GeV/c]", "entries");

  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPtJet1, xaxis, "p_{T} [GeV/c]", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPtJet2, xaxis, "p_{T} [GeV/c]", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPtJet, xaxis, "p_{T} [GeV/c]", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi1, xaxis, "#xi", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi2, xaxis, "#xi", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi, xaxis, "#xi", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z1, xaxis, "z", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z2, xaxis, "z", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z, xaxis, "z", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Pt1, xaxis, "p_{T} [GeV/c]", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Pt2, xaxis, "p_{T} [GeV/c]", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Pt, xaxis, "p_{T} [GeV/c]", "Entries");
}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::FillDiJetFF(Int_t jetType, Float_t trackPt, Float_t jetPt, Double_t jetBin, Bool_t incrementJetPt)
{
  // fill DiJet FF

  if(jetType == 0)
    {
      if(incrementJetPt) fh1JetPt->Fill(jetPt);  
      
      fh2TrackPtJet->Fill(jetBin, trackPt);
      
      Double_t z = trackPt / jetPt;
      Double_t xi = 0;
      if(z>0) xi = TMath::Log(1/z);
      
      fh2Xi->Fill(jetBin, xi);
      fh2Z->Fill(jetBin, z);
    }
  if(jetType == 1)
    {
      if(incrementJetPt) fh1Jet1Pt->Fill(jetPt);
      
      fh2TrackPtJet1->Fill(jetBin, trackPt);
      
      Double_t z = trackPt / jetPt;
      Double_t xi = 0;
      if(z>0) xi = TMath::Log(1/z);
      
      fh2Xi1->Fill(jetBin, xi);
      fh2Z1->Fill(jetBin, z);
    }
  if(jetType == 2)
    {
      if(incrementJetPt) fh1Jet2Pt->Fill(jetPt);
      
      fh2TrackPtJet2->Fill(jetBin, trackPt);
      
      Double_t z = trackPt / jetPt;
      Double_t xi = 0;
      if(z>0) xi = TMath::Log(1/z);
      
      fh2Xi2->Fill(jetBin, xi);
      fh2Z2->Fill(jetBin, z);
    }


}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncDiJetHistos::AddToOutput(TList* list)const
{
  // add histos to list

  list->Add(fh1Jet1Pt);
  list->Add(fh1Jet2Pt);
  list->Add(fh1JetPt);
  list->Add(fh2TrackPtJet1);
  list->Add(fh2TrackPtJet2);
  list->Add(fh2TrackPtJet);
  list->Add(fh2Xi1);
  list->Add(fh2Xi2);
  list->Add(fh2Xi);
  list->Add(fh2Z1);
  list->Add(fh2Z2);
  list->Add(fh2Z);
}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::AliFragFuncQADiJetHistos(const char* name, Int_t kindSlices,
					     Int_t nInvMass, Float_t invMassMin, Float_t invMassMax, 
					     Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
					     Int_t nDeltaPhi, Float_t deltaPhiMin, Float_t deltaPhiMax, 
					     Int_t nDeltaEta, Float_t deltaEtaMin, Float_t deltaEtaMax, 
					     Int_t nDeltaPt, Float_t deltaPtMin, Float_t deltaPtMax, 
                                             Int_t nInBal, Float_t inBalMin, Float_t inBalMax)
  : TObject()
  ,fKindSlices(kindSlices)
  ,fNBinsJetInvMass(nInvMass)
  ,fJetInvMassMin(invMassMin)
  ,fJetInvMassMax(invMassMax)
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsDeltaPhi(nDeltaPhi)
  ,fDeltaPhiMin(deltaPhiMin)
  ,fDeltaPhiMax(deltaPhiMax)
  ,fNBinsDeltaEta(nDeltaEta)
  ,fDeltaEtaMin(deltaEtaMin)
  ,fDeltaEtaMax(deltaEtaMax)
  ,fNBinsDeltaPt(nDeltaPt)
  ,fDeltaPtMin(deltaPtMin)
  ,fDeltaPtMax(deltaPtMax)
  ,fNBinsInBal(nInBal)
  ,fInBalMin(inBalMin)
  ,fInBalMax(inBalMax)
  ,fh2InvMass(0)
  ,fh2DeltaPhi(0)
  ,fh2DeltaEta(0)
  ,fh2DeltaPt(0)
  ,fh2InBal(0)
  ,fNameQADJ(name)
{
  // default constructor

}

//______________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::AliFragFuncQADiJetHistos(const AliFragFuncQADiJetHistos& copy)
  : TObject()
  ,fKindSlices(copy.fKindSlices)
  ,fNBinsJetInvMass(copy.fNBinsJetInvMass)
  ,fJetInvMassMin(copy.fJetInvMassMin)
  ,fJetInvMassMax(copy.fJetInvMassMax)
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsDeltaPhi(copy.fNBinsDeltaPhi)
  ,fDeltaPhiMin(copy.fDeltaPhiMin)
  ,fDeltaPhiMax(copy.fDeltaPhiMax)
  ,fNBinsDeltaEta(copy.fNBinsDeltaEta)
  ,fDeltaEtaMin(copy.fDeltaEtaMin)
  ,fDeltaEtaMax(copy.fDeltaEtaMax)
  ,fNBinsDeltaPt(copy.fNBinsDeltaPt)
  ,fDeltaPtMin(copy.fDeltaPtMin)
  ,fDeltaPtMax(copy.fDeltaPtMax)
  ,fNBinsInBal(copy.fNBinsInBal)
  ,fInBalMin(copy.fInBalMin)
  ,fInBalMax(copy.fInBalMax)
  ,fh2InvMass(copy.fh2InvMass)
  ,fh2DeltaPhi(copy.fh2DeltaPhi)
  ,fh2DeltaEta(copy.fh2DeltaEta)
  ,fh2DeltaPt(copy.fh2DeltaPt)
  ,fh2InBal(copy.fh2InBal)
  ,fNameQADJ(copy.fNameQADJ)
{
  // default constructor

}

//_______________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos& AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::operator=(const AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fKindSlices       = o.fKindSlices;
    fNBinsJetInvMass  = o.fNBinsJetInvMass;
    fJetInvMassMin    = o.fJetInvMassMin;
    fJetInvMassMax    = o.fJetInvMassMax;
    fNBinsJetPt       = o.fNBinsJetPt;
    fJetPtMin         = o.fJetPtMin;
    fJetPtMax         = o.fJetPtMax;
    fNBinsDeltaPhi    = o.fNBinsDeltaPhi;
    fDeltaPhiMin      = o.fDeltaPhiMin;
    fDeltaPhiMax      = o.fDeltaPhiMax;
    fNBinsDeltaEta    = o.fNBinsDeltaEta;
    fDeltaEtaMin      = o.fDeltaEtaMin;
    fDeltaEtaMax      = o.fDeltaEtaMax;
    fNBinsDeltaPt     = o.fNBinsDeltaPt;
    fDeltaPtMin       = o.fDeltaPtMin;
    fDeltaPtMax       = o.fDeltaPtMax;
    fNBinsInBal       = o.fNBinsInBal;
    fInBalMin         = o.fInBalMin;
    fInBalMax         = o.fInBalMax;
    fh2InvMass        = o.fh2InvMass;
    fh2DeltaPhi       = o.fh2DeltaPhi;
    fh2DeltaEta       = o.fh2DeltaEta;
    fh2DeltaPt        = o.fh2DeltaPt;
    fh2InBal          = o.fh2InBal;
    fNameQADJ         = o.fNameQADJ;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::~AliFragFuncQADiJetHistos()
{
  // destructor 

  if(fh2InvMass)  delete fh2InvMass;
  if(fh2DeltaPhi) delete fh2DeltaPhi;
  if(fh2DeltaEta) delete fh2DeltaEta;
  if(fh2DeltaPt)  delete fh2DeltaPt;
  if(fh2InBal)    delete fh2InBal;
}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::DefineQADiJetHistos()
{
  // define histos
  
  Int_t nBins = 0;
  Double_t min = 0.;
  Double_t max = 0.;
  const char *xaxis = "";
  if(fKindSlices == 1)
    {
      nBins = fNBinsJetInvMass;
      min   = fJetInvMassMin;
      max   = fJetInvMassMax;
      xaxis = "M_{JJ} [GeV]";
    }
  if(fKindSlices == 2 || fKindSlices == 3)
    {
      nBins = fNBinsJetPt;
      min   = fJetPtMin;
      max   = fJetPtMax;
      if(fKindSlices == 2) xaxis = "E_{Tmean} [GeV]";
      if(fKindSlices == 3) xaxis ="leading jet p_{T} [GeV/c]";
    }
  
  
  fh2InvMass  = new TH2F(Form("fh2DJInvMassPositionCut%s",  fNameQADJ.Data()), "",nBins, min, max, fNBinsJetInvMass, fJetInvMassMin, fJetInvMassMax);
  fh2DeltaPhi = new TH2F(Form("fh2DJDeltaPhiPositionCut%s", fNameQADJ.Data()), "",nBins, min, max, fNBinsDeltaPhi, fDeltaPhiMin, fDeltaPhiMax);
  fh2DeltaEta = new TH2F(Form("fh2DJDeltaEtaPositionCut%s", fNameQADJ.Data()), "",nBins, min, max, fNBinsDeltaEta, fDeltaEtaMin, fDeltaEtaMax);
  fh2DeltaPt  = new TH2F(Form("fh2DJDeltaPtPositionCut%s",  fNameQADJ.Data()), "",nBins, min, max, fNBinsDeltaPt, fDeltaPtMin, fDeltaPtMax);
  fh2InBal  = new TH2F(Form("fh2DJInBalPositionCut%s",  fNameQADJ.Data()), "",nBins, min, max, fNBinsInBal, fInBalMin, fInBalMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh2InvMass, xaxis, "Invariant Mass", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaPhi, xaxis, "#Delta #phi", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaEta, xaxis, "#Delta #eta", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaPt, xaxis, "#Delta p_{T}", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2InBal, xaxis, "(p_{T}^{1}-p_{T}^{2})/(p_{T}^{1}+p_{T}^{2})", "Entries");

}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::FillDiJetQA(Double_t invMass, Double_t deltaPhi, Double_t deltaEta,Double_t deltaPt, Double_t inbal, Double_t jetBin)
{
  // fill dijet QA

  fh2InvMass->Fill(jetBin, invMass);
  fh2DeltaPhi->Fill(jetBin, deltaPhi);
  fh2DeltaEta->Fill(jetBin, deltaEta);
  fh2DeltaPt->Fill(jetBin, deltaPt);
  fh2InBal->Fill(jetBin, inbal);
}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::AddToOutput(TList* list)const
{
  // add histos to list

  list->Add(fh2InvMass);
  list->Add(fh2DeltaPhi);
  list->Add(fh2DeltaEta);
  list->Add(fh2DeltaPt);
  list->Add(fh2InBal);
}

//_________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2CosTheta);
  list->Add(fh2PtZ);
  list->Add(fh3ThetaZ);
  list->Add(fh3JtTheta);
  list->Add(fh3JtZ);

}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskFragmentationFunction::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // (taken from AliAnalysisTaskJetSpectrum2)
  // 
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }

  // Set seed for backg study
  fRandom = new TRandom3();
  fRandom->SetSeed(0);

  return kTRUE;
}



//__________________________________________________________________
void AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()
{
  // create output objects

  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()");
 
  // create list of tracks and jets 
  
  fTracksRec = new TList();
  fTracksRec->SetOwner(kFALSE);  

  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE);  

  fTracksGen = new TList();
  fTracksGen->SetOwner(kFALSE);

  fTracksAODMCCharged = new TList();
  fTracksAODMCCharged->SetOwner(kFALSE);
    
  fTracksRecQualityCuts = new TList(); 
  fTracksRecQualityCuts->SetOwner(kFALSE);

  fJetsRec = new TList();
  fJetsRec->SetOwner(kFALSE);
  if(fBranchRecJets.Contains("KT") && fBckgSubMethod) fJetsRec->SetOwner(kTRUE);

  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);
  if(fBranchRecJets.Contains("KT") && fBckgSubMethod) fJetsRecCuts->SetOwner(kTRUE);

  fJetsGen = new TList();
  fJetsGen->SetOwner(kFALSE);

  fJetsRecEff = new TList();
  fJetsRecEff->SetOwner(kFALSE);

  // fJetsKine = new TList();
  // fJetsKine->SetOwner(kTRUE); // delete AOD jets using mom from Kine Tree via TList::Clear()

  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters ||  fBckgType[3]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
      fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)){
    
    fBckgJetsRec = new TList();
    fBckgJetsRec->SetOwner(kFALSE);

    fBckgJetsRecCuts = new TList();
    fBckgJetsRecCuts->SetOwner(kFALSE);

    fBckgJetsGen = new TList();
    fBckgJetsGen->SetOwner(kFALSE);
  }

  //
  // Create histograms / output container
  //

  OpenFile(1);
  fCommonHistList = new TList();
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  
  // Histograms	
  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1EvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fh1EvtSelection->GetXaxis()->SetBinLabel(2,"event selection: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(3,"event class: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(4,"vertex Ncontr: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(5,"vertex z: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(6,"vertex type: rejected");
  
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 2500,-.5, 2499.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1EvtMult 	             = new TH1F("fh1EvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",120,0.,12000.);
  fh1EvtCent 	             = new TH1F("fh1EvtCent","centrality",100,0.,100.);
  fh1Xsec                    = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Trials                  = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1PtHard                  = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fh1PtHardTrials            = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);

  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",10,-0.5,9.5);
  fh1nGenJets                = new TH1F("fh1nGenJets","generated jets per event",10,-0.5,9.5);
  fh1nRecEffJets             = new TH1F("fh1nRecEffJets","reconstruction effiency: jets per event",10,-0.5,9.5);
  fh2PtRecVsGenPrim          = new TH2F("fh2PtRecVsGenPrim","rec vs gen pt",fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax,fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax);

  // Background
  if(fBckgMode) {
    if(fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || 
       fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
       fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading){
      
      fh1nRecBckgJetsCuts        = new TH1F("fh1nRecBckgJetsCuts","reconstructed background jets per event",10,-0.5,9.5);
      fh1nGenBckgJets            = new TH1F("fh1nGenBckgJets","generated background jets per event",10,-0.5,9.5);
    }

    if(fBckgType[0]==kBckgPerp || fBckgType[1]==kBckgPerp || fBckgType[2]==kBckgPerp || fBckgType[3]==kBckgPerp)
      fh1PerpMult                = new TH1F("fh1PerpMult","Background multiplicity - Cone perpendicular to leading jet axis",500,0.,500.);
    if(fBckgType[0]==kBckgASide || fBckgType[1]==kBckgASide || fBckgType[2]==kBckgASide || fBckgType[3]==kBckgASide)
      fh1ASideMult               = new TH1F("fh1ASideMult","Background multiplicity - Cone in the away side of leading jet axis",500,0.,500.);
    if(fBckgType[0]==kBckgASideWindow || fBckgType[1]==kBckgASideWindow || fBckgType[2]==kBckgASideWindow || fBckgType[3]==kBckgASideWindow)
      fh1ASideWindowMult         = new TH1F("fh1ASideWindowMult","Background multiplicity - Cone in the away side of leading jet axis",500,0.,500.);
    if(fBckgType[0]==kBckgPerpWindow || fBckgType[1]==kBckgPerpWindow || fBckgType[2]==kBckgPerpWindow || fBckgType[3]==kBckgPerpWindow)
      fh1PerpWindowMult         = new TH1F("fh1PerpWindowMult","Background multiplicity - Cone in the perp direction of leading jet axis",500,0.,500.);
    if(fBckgType[0]==kBckgOutLJ || fBckgType[1]==kBckgOutLJ || fBckgType[2]==kBckgOutLJ || fBckgType[3]==kBckgOutLJ)
      fh1OutLeadingMult          = new TH1F("fh1OutLeadingMult","Background multiplicity - Cone outside leading jet",500,0,500.);
    if(fBckgType[0]==kBckgOutLJStat || fBckgType[1]==kBckgOutLJStat || fBckgType[2]==kBckgOutLJStat || fBckgType[3]==kBckgOutLJStat)
      fh1OutLeadingStatMult      = new TH1F("fh1OutLeadingStatMult","Background multiplicity - Cone outside leading jet",3000,0,3000.);
    if(fBckgType[0]==kBckgOut2J || fBckgType[1]==kBckgOut2J || fBckgType[2]==kBckgOut2J || fBckgType[3]==kBckgOut2J)
      fh1Out2JetsMult            = new TH1F("fh1Out2JetsMult","Background multiplicity - Cone outside 2 jets",500,0.,500.);
    if(fBckgType[0]==kBckgOut3J || fBckgType[1]==kBckgOut3J || fBckgType[2]==kBckgOut3J || fBckgType[3]==kBckgOut3J)
      fh1Out3JetsMult            = new TH1F("fh1Out3JetsMult","Background multiplicity - Cone outside 3 jets",500,0.,500.);
    if(fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters)
      fh1MedianClustersMult  = new TH1F("fh1MedianClustersMult","Background multiplicity - median cluster",500,0.,500.);
    if(fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)
      fh1OutClustersMult            = new TH1F("fh1OutClustersMult","Background multiplicity - clusters outside leading jet",3000,0.,3000.);
  }
  
  if(fQAMode){
    if(fQAMode&1){ // track QA
      fQATrackHistosRec          = new AliFragFuncQATrackHistos("Rec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								fQATrackHighPtThreshold);
      fQATrackHistosRecCuts      = new AliFragFuncQATrackHistos("RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								fQATrackHighPtThreshold);
      fQATrackHistosGen          = new AliFragFuncQATrackHistos("Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								fQATrackHighPtThreshold);
    }

    if(fQAMode&2){ // jet QA
      fQAJetHistosRec            = new AliFragFuncQAJetHistos("Rec", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							      fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
							      fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
      fQAJetHistosRecCuts        = new AliFragFuncQAJetHistos("RecCuts", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							      fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
							      fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
      fQAJetHistosRecCutsLeading = new AliFragFuncQAJetHistos("RecCutsLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							      fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
							      fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
      fQAJetHistosGen            = new AliFragFuncQAJetHistos("Gen", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							      fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
							      fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
      fQAJetHistosGenLeading     = new AliFragFuncQAJetHistos("GenLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							      fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,
							      fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);  
      if(fEffMode) fQAJetHistosRecEffLeading  = new AliFragFuncQAJetHistos("RecEffLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
									   fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);
    }
  } // end: QA

  if(fFFMode){
    fFFHistosRecCuts   	     = new AliFragFuncHistos("RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFHistosRecLeading        = new AliFragFuncHistos("RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						       fFFNBinsPt, fFFPtMin, fFFPtMax, 
						       fFFNBinsXi, fFFXiMin, fFFXiMax,  
						       fFFNBinsZ , fFFZMin , fFFZMax);
    fFFHistosRecLeadingTrack   = new AliFragFuncHistos("RecLeadingTrack", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						       fFFNBinsPt, fFFPtMin, fFFPtMax, 
						       fFFNBinsXi, fFFXiMin, fFFXiMax,  
						       fFFNBinsZ , fFFZMin , fFFZMax);
    fFFHistosGen   	     = new AliFragFuncHistos("Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFHistosGenLeading        = new AliFragFuncHistos("GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						       fFFNBinsPt, fFFPtMin, fFFPtMax, 
						       fFFNBinsXi, fFFXiMin, fFFXiMax,  
						       fFFNBinsZ , fFFZMin , fFFZMax);
    fFFHistosGenLeadingTrack   = new AliFragFuncHistos("GenLeadingTrack", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						       fFFNBinsPt, fFFPtMin, fFFPtMax, 
						       fFFNBinsXi, fFFXiMin, fFFXiMax,  
						       fFFNBinsZ , fFFZMin , fFFZMax);
  } // end: FF

  if(fIJMode)
    {
      fIJHistosRecCuts   	     = new AliFragFuncIntraJetHistos("RecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJHistosRecLeading        = new AliFragFuncIntraJetHistos("RecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								 fIJNBinsPt, fIJPtMin, fIJPtMax, 
								 fIJNBinsZ, fIJZMin, fIJZMax,  
								 fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								 fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								 fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJHistosRecLeadingTrack   = new AliFragFuncIntraJetHistos("RecLeadingTrack", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								 fIJNBinsPt, fIJPtMin, fIJPtMax, 
								 fIJNBinsZ, fIJZMin, fIJZMax,  
								 fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								 fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								 fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJHistosGen   	     = new AliFragFuncIntraJetHistos("Gen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
							     fIJNBinsPt, fIJPtMin, fIJPtMax, 
							     fIJNBinsZ, fIJZMin, fIJZMax,  
							     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
							     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
							     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJHistosGenLeading        = new AliFragFuncIntraJetHistos("GenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								 fIJNBinsPt, fIJPtMin, fIJPtMax, 
								 fIJNBinsZ, fIJZMin, fIJZMax,  
								 fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								 fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								 fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJHistosGenLeadingTrack   = new AliFragFuncIntraJetHistos("GenLeadingTrack", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								 fIJNBinsPt, fIJPtMin, fIJPtMax, 
								 fIJNBinsZ, fIJZMin, fIJZMax,  
								 fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								 fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								 fIJNBinsJt , fIJJtMin , fIJJtMax);
    } // end: intra-jet
  
  if(fDJMode){
    if(fDJMode&1){
      fFFDiJetHistosRecCuts         = new AliFragFuncDiJetHistos("RecCuts", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax,
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax, 
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax, 
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax);
      fFFDiJetHistosRecLeading      = new AliFragFuncDiJetHistos("RecLeading", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax, 
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax, 
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax, 
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax); 
      fFFDiJetHistosRecLeadingTrack = new AliFragFuncDiJetHistos("RecLeadingTrack", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax, 
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax, 
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax, 
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax);
      
      fFFDiJetHistosGen             = new AliFragFuncDiJetHistos("Gen", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax,
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax, 
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax,
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax);
      fFFDiJetHistosGenLeading      = new AliFragFuncDiJetHistos("GenLeading", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax,
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax,
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax,
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax);
      fFFDiJetHistosGenLeadingTrack = new AliFragFuncDiJetHistos("GenLeadingTrack", fDiJetKindBins, 
								 fDiJetNBinsJetInvMass, fDiJetJetInvMassMin, fDiJetJetInvMassMax,
								 fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax, 
								 fDiJetNBinsPt, fDiJetPtMin, fDiJetPtMax, 
								 fDiJetNBinsXi, fDiJetXiMin, fDiJetXiMax, 
								 fDiJetNBinsZ, fDiJetZMin, fDiJetZMax);
    }
    
    if(fDJMode&2){
      fQADiJetHistosRecCuts = new AliFragFuncQADiJetHistos("RecCuts", fDiJetKindBins, 
							   fQADiJetNBinsInvMass, fQADiJetInvMassMin, fQADiJetInvMassMax,
							   fQADiJetNBinsJetPt, fQADiJetJetPtMin, fQADiJetJetPtMax,
							   fQADiJetNBinsDeltaPhi, fQADiJetDeltaPhiMin, fQADiJetDeltaPhiMax , 
							   fQADiJetNBinsDeltaEta, fQADiJetDeltaEtaMin, fQADiJetDeltaEtaMax , 
							   fQADiJetNBinsDeltaPt, fQADiJetDeltaPtMin, fQADiJetDeltaPtMax,
							   fQADiJetNBinsInBal, fQADiJetInBalMin, fQADiJetInBalMax);
      fQADiJetHistosGen     = new AliFragFuncQADiJetHistos("Gen", fDiJetKindBins, 
							   fQADiJetNBinsInvMass, fQADiJetInvMassMin,  fQADiJetInvMassMax, 
							   fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax,
							   fQADiJetNBinsDeltaPhi, fQADiJetDeltaPhiMin, fQADiJetDeltaPhiMax,
							   fQADiJetNBinsDeltaEta, fQADiJetDeltaEtaMin, fQADiJetDeltaEtaMax,
							   fQADiJetNBinsDeltaPt, fQADiJetDeltaPtMin, fQADiJetDeltaPtMax,
							   fQADiJetNBinsInBal, fQADiJetInBalMin, fQADiJetInBalMax);
    }
  } // end: di-jet

  // efficiency

  if(fEffMode){
    if(fQAMode&1){
      fQATrackHistosRecEffGen = new AliFragFuncQATrackHistos("RecEffGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							     fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							     fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							     fQATrackHighPtThreshold);
      
      fQATrackHistosRecEffRec = new AliFragFuncQATrackHistos("RecEffRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							     fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							     fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							     fQATrackHighPtThreshold);


      Int_t    nBinsResponseSinglePt[2]     = {fFFNBinsPt, fFFNBinsPt};
      Double_t binMinResponseSinglePt[2]    = {fFFPtMin, fFFPtMin};
      Double_t binMaxResponseSinglePt[2]    = {fFFPtMax, fFFPtMax};
      const char* labelsResponseSinglePt[2] = {"rec p_{T} [GeV/c]", "gen p_{T} [GeV/c]"};

      fhnResponseSinglePt  = new THnSparseF("fhnResponseSinglePt","track pt gen : track pt rec",2,
					    nBinsResponseSinglePt,binMinResponseSinglePt,binMaxResponseSinglePt);
     
      AliAnalysisTaskFragmentationFunction::SetProperties(fhnResponseSinglePt,2,labelsResponseSinglePt);
    }
    if(fFFMode){
      fFFHistosRecEffGen      = new AliFragFuncHistos("RecEffGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFHistosRecEffRec      = new AliFragFuncHistos("RecEffRec", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);


      Int_t    nBinsResponseJetTrackPt[3]     = {fFFNBinsJetPt,fFFNBinsPt,fFFNBinsPt};
      Double_t binMinResponseJetTrackPt[3]    = {fFFJetPtMin,fFFPtMin, fFFPtMin};
      Double_t binMaxResponseJetTrackPt[3]    = {fFFJetPtMax,fFFPtMax,fFFPtMax};
      const char* labelsResponseJetTrackPt[3] = { "jet p_{T} [GeV/c]","rec p_{T} [GeV/c]", "gen p_{T} [GeV/c]"};

      fhnResponseJetTrackPt  = new THnSparseF("fhnResponseJetTrackPt","jet pt:track pt rec:track pt gen",3,
					      nBinsResponseJetTrackPt,binMinResponseJetTrackPt,binMaxResponseJetTrackPt);
     
      AliAnalysisTaskFragmentationFunction::SetProperties(fhnResponseJetTrackPt,3,labelsResponseJetTrackPt);

      Int_t    nBinsResponseJetZ[3]     = {fFFNBinsJetPt, fFFNBinsZ,fFFNBinsZ};
      Double_t binMinResponseJetZ[3]    = {fFFJetPtMin, fFFZMin, fFFZMin};
      Double_t binMaxResponseJetZ[3]    = {fFFJetPtMax, fFFZMax, fFFZMax};
      const char* labelsResponseJetZ[3] = { "jet p_{T} [GeV/c]","rec z","gen z"};

      fhnResponseJetZ  = new THnSparseF("fhnResponseJetZ","jet pt:track pt rec:track pt gen",3,
					nBinsResponseJetZ,binMinResponseJetZ,binMaxResponseJetZ);
      
      AliAnalysisTaskFragmentationFunction::SetProperties(fhnResponseJetZ,3,labelsResponseJetZ);
      
      Int_t    nBinsResponseJetXi[3]     = {fFFNBinsJetPt, fFFNBinsXi,fFFNBinsXi};
      Double_t binMinResponseJetXi[3]    = {fFFJetPtMin, fFFXiMin, fFFXiMin};
      Double_t binMaxResponseJetXi[3]    = {fFFJetPtMax, fFFXiMax, fFFXiMax};
      const char* labelsResponseJetXi[3] = { "jet p_{T} [GeV/c]","rec xi","gen xi"};

      fhnResponseJetXi  = new THnSparseF("fhnResponseJetXi","jet pt:track xi rec:track xi gen",3,
					nBinsResponseJetXi,binMinResponseJetXi,binMaxResponseJetXi);
      
      AliAnalysisTaskFragmentationFunction::SetProperties(fhnResponseJetXi,3,labelsResponseJetXi);

    }
  } // end: efficiency

  // Background
  if(fBckgMode){
    // Track QA
    TString title[4];
    for(Int_t i=0; i<4; i++){
      if(fBckgType[i]==kBckgPerp) title[i]="Perp";
      else if(fBckgType[i]==kBckgPerpWindow) title[i]="PerpW";
      else if(fBckgType[i]==kBckgASide) title[i]="ASide";
      else if(fBckgType[i]==kBckgASideWindow) title[i]="ASideW";
      else if(fBckgType[i]==kBckgOutLJ) title[i]="OutLeadingJet";
      else if(fBckgType[i]==kBckgOut2J) title[i]="Out2Jets";
      else if(fBckgType[i]==kBckgOut3J) title[i]="Out3Jets";
      else if(fBckgType[i]==kBckgOutAJ) title[i]="AllJets";
      else if(fBckgType[i]==kBckgOutLJStat) title[i]="OutLeadingJetStat";
      else if(fBckgType[i]==kBckgOut2JStat) title[i]="Out2JetsStat";
      else if(fBckgType[i]==kBckgOut3JStat) title[i]="Out3JetsStat";
      else if(fBckgType[i]==kBckgOutAJStat) title[i]="AllJetsStat";
      else if(fBckgType[i]==kBckgClustersOutLeading) title[i]="OutClusters";
      else if(fBckgType[i]==kBckgClusters) title[i]="MedianClusters";
      else printf("Please chose background method number %d!",i);
    }

    if(fQAMode&1){
      fQABckgHisto0RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[0]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto0Gen          = new AliFragFuncQATrackHistos("Bckg"+title[0]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto1RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[1]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto1Gen          = new AliFragFuncQATrackHistos("Bckg"+title[1]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto2RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[2]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto2Gen          = new AliFragFuncQATrackHistos("Bckg"+title[2]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto3RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[3]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto3Gen          = new AliFragFuncQATrackHistos("Bckg"+title[3]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
    } // end: background QA

    if(fFFMode){
      // outside leading jet or 2 jets or more
      fFFBckgHisto0RecCuts    = new AliFragFuncHistos("Bckg"+title[0]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto0RecLeading = new AliFragFuncHistos("Bckg"+title[0]+"RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto0Gen        = new AliFragFuncHistos("Bckg"+title[0]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto0GenLeading = new AliFragFuncHistos("Bckg"+title[0]+"GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFBckgHisto1RecCuts    = new AliFragFuncHistos("Bckg"+title[1]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto1RecLeading = new AliFragFuncHistos("Bckg"+title[1]+"RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto1Gen        = new AliFragFuncHistos("Bckg"+title[1]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto1GenLeading = new AliFragFuncHistos("Bckg"+title[1]+"GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFBckgHisto2RecCuts    = new AliFragFuncHistos("Bckg"+title[2]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto2RecLeading = new AliFragFuncHistos("Bckg"+title[2]+"RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto2Gen        = new AliFragFuncHistos("Bckg"+title[2]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto2GenLeading = new AliFragFuncHistos("Bckg"+title[2]+"GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto3RecCuts    = new AliFragFuncHistos("Bckg"+title[3]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto3RecLeading = new AliFragFuncHistos("Bckg"+title[3]+"RecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto3Gen        = new AliFragFuncHistos("Bckg"+title[3]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      fFFBckgHisto3GenLeading = new AliFragFuncHistos("Bckg"+title[3]+"GenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
    } // end: background FF

    if(fIJMode){    
      fIJBckgHisto0RecCuts   	     = new AliFragFuncIntraJetHistos("Bckg"+title[0]+"RecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto0RecLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[0]+"RecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto0Gen   	     = new AliFragFuncIntraJetHistos("Bckg"+title[0]+"Gen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto0GenLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[0]+"GenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      
      
      fIJBckgHisto1RecCuts   	     = new AliFragFuncIntraJetHistos("Bckg"+title[1]+"RecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto1RecLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[1]+"RecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto1Gen   	     = new AliFragFuncIntraJetHistos("Bckg"+title[1]+"Gen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto1GenLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[1]+"GenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      
      fIJBckgHisto2RecCuts   	     = new AliFragFuncIntraJetHistos("Bckg"+title[2]+"RecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto2RecLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[2]+"RecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto2Gen   	     = new AliFragFuncIntraJetHistos("Bckg"+title[2]+"Gen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
      fIJBckgHisto2GenLeading        = new AliFragFuncIntraJetHistos("Bckg"+title[2]+"GenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    } // end: background intra-jet
  } // end: background

  if(fPhiCorrMode){
      fPhiCorrHistosJetArea = new AliFragFuncQATrackHistos("JetArea",  fPhiCorrNBinsPt,  fPhiCorrPtMin,  fPhiCorrPtMax,
	      fPhiCorrNBinsEta, fPhiCorrEtaMin, fPhiCorrEtaMax,
	      fPhiCorrNBinsPhi, fPhiCorrPhiMin, fPhiCorrPhiMax,
	      fQATrackHighPtThreshold);

      fPhiCorrHistosTransverseArea = new AliFragFuncQATrackHistos("TransverseArea", fPhiCorrNBinsPt, fPhiCorrPtMin, fPhiCorrPtMax,
	      fPhiCorrNBinsEta, fPhiCorrEtaMin, fPhiCorrEtaMax,
	      fPhiCorrNBinsPhi, fPhiCorrPhiMin, fPhiCorrPhiMax,
	      fQATrackHighPtThreshold);

      fPhiCorrHistosAwayArea = new AliFragFuncQATrackHistos("AwayArea", fPhiCorrNBinsPt, fPhiCorrPtMin, fPhiCorrPtMax,
	      fPhiCorrNBinsEta, fPhiCorrEtaMin, fPhiCorrEtaMax,
	      fPhiCorrNBinsPhi, fPhiCorrPhiMin, fPhiCorrPhiMax,
	      fQATrackHighPtThreshold);
  } // end: phi correlation

  
  // ____________ define histograms ____________________
  
  if(fQAMode){
    if(fQAMode&1){ // track QA
      fQATrackHistosRec->DefineHistos();
      fQATrackHistosRecCuts->DefineHistos();
      fQATrackHistosGen->DefineHistos();
    }

    if(fQAMode&2){ // jet QA
      fQAJetHistosRec->DefineHistos();
      fQAJetHistosRecCuts->DefineHistos();
      fQAJetHistosRecCutsLeading->DefineHistos();
      fQAJetHistosGen->DefineHistos();
      fQAJetHistosGenLeading->DefineHistos();
      if(fEffMode) fQAJetHistosRecEffLeading->DefineHistos();
    }
  }

  if(fFFMode){
    fFFHistosRecCuts->DefineHistos();
    fFFHistosRecLeading->DefineHistos();
    fFFHistosRecLeadingTrack->DefineHistos();
    fFFHistosGen->DefineHistos();
    fFFHistosGenLeading->DefineHistos();
    fFFHistosGenLeadingTrack->DefineHistos();
  }

  if(fIJMode){
    fIJHistosRecCuts->DefineHistos();
    fIJHistosRecLeading->DefineHistos();
    fIJHistosRecLeadingTrack->DefineHistos();
    fIJHistosGen->DefineHistos();
    fIJHistosGenLeading->DefineHistos();
    fIJHistosGenLeadingTrack->DefineHistos();
  }

  if(fDJMode){
    if(fDJMode&1){
      fFFDiJetHistosRecCuts->DefineDiJetHistos();
      fFFDiJetHistosRecLeading->DefineDiJetHistos();
      fFFDiJetHistosRecLeadingTrack->DefineDiJetHistos();
      fFFDiJetHistosGen->DefineDiJetHistos();
      fFFDiJetHistosGenLeading->DefineDiJetHistos();
      fFFDiJetHistosGenLeadingTrack->DefineDiJetHistos();
    }
    
    if(fDJMode&2){
      fQADiJetHistosRecCuts->DefineQADiJetHistos();
      fQADiJetHistosGen->DefineQADiJetHistos();
    }
  } // end: di-jet

  if(fEffMode){
    if(fQAMode&1){
      fQATrackHistosRecEffGen->DefineHistos();
      fQATrackHistosRecEffRec->DefineHistos(); 
    }
    if(fFFMode){
      fFFHistosRecEffGen->DefineHistos();
      fFFHistosRecEffRec->DefineHistos();
    }
  } // end: efficiency

  // Background
  if(fBckgMode){
    if(fFFMode){
      fFFBckgHisto0RecCuts->DefineHistos();
      fFFBckgHisto0RecLeading->DefineHistos();
      fFFBckgHisto0Gen->DefineHistos();
      fFFBckgHisto0GenLeading->DefineHistos();
      fFFBckgHisto1RecCuts->DefineHistos();
      fFFBckgHisto1RecLeading->DefineHistos();
      fFFBckgHisto1Gen->DefineHistos();
      fFFBckgHisto1GenLeading->DefineHistos();
      fFFBckgHisto2RecCuts->DefineHistos();
      fFFBckgHisto2RecLeading->DefineHistos();
      fFFBckgHisto2Gen->DefineHistos();
      fFFBckgHisto2GenLeading->DefineHistos();
      fFFBckgHisto3RecCuts->DefineHistos();
      fFFBckgHisto3RecLeading->DefineHistos();
      fFFBckgHisto3Gen->DefineHistos();
      fFFBckgHisto3GenLeading->DefineHistos();
    }

    if(fIJMode){
      fIJBckgHisto0RecCuts->DefineHistos();
      fIJBckgHisto0RecLeading->DefineHistos();
      fIJBckgHisto0Gen->DefineHistos();
      fIJBckgHisto0GenLeading->DefineHistos();
      fIJBckgHisto1RecCuts->DefineHistos();
      fIJBckgHisto1RecLeading->DefineHistos();
      fIJBckgHisto1Gen->DefineHistos();
      fIJBckgHisto1GenLeading->DefineHistos();
      fIJBckgHisto2RecCuts->DefineHistos();
      fIJBckgHisto2RecLeading->DefineHistos();
      fIJBckgHisto2Gen->DefineHistos();
      fIJBckgHisto2GenLeading->DefineHistos();
    }

    if(fQAMode&1){
      fQABckgHisto0RecCuts->DefineHistos();
      fQABckgHisto0Gen->DefineHistos();
      fQABckgHisto1RecCuts->DefineHistos();
      fQABckgHisto1Gen->DefineHistos();
      fQABckgHisto2RecCuts->DefineHistos();
      fQABckgHisto2Gen->DefineHistos();
      fQABckgHisto3RecCuts->DefineHistos();
      fQABckgHisto3Gen->DefineHistos();
    }
  } // end: background
  
  if(fPhiCorrMode){
      fPhiCorrHistosJetArea->DefineHistos();
      fPhiCorrHistosTransverseArea->DefineHistos();
      fPhiCorrHistosAwayArea->DefineHistos();
  }

  Bool_t genJets    = (fJetTypeGen != kJetsUndef) ? kTRUE : kFALSE;
  Bool_t genTracks  = (fTrackTypeGen != kTrackUndef) ? kTRUE : kFALSE;
  Bool_t recJetsEff = (fJetTypeRecEff != kJetsUndef) ? kTRUE : kFALSE;

  fCommonHistList->Add(fh1EvtSelection);
  fCommonHistList->Add(fh1EvtMult);
  fCommonHistList->Add(fh1EvtCent);
  fCommonHistList->Add(fh1VertexNContributors);
  fCommonHistList->Add(fh1VertexZ);    
  fCommonHistList->Add(fh1nRecJetsCuts);
  if(genJets && genTracks){
    fCommonHistList->Add(fh1Xsec);
    fCommonHistList->Add(fh1Trials);
    fCommonHistList->Add(fh1PtHard);
    fCommonHistList->Add(fh1PtHardTrials);
    if(genJets) fCommonHistList->Add(fh1nGenJets);
  }

  // FF histograms
  if(fFFMode){
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecLeading->AddToOutput(fCommonHistList);
    fFFHistosRecLeadingTrack->AddToOutput(fCommonHistList);
    if(genJets && genTracks){
      fCommonHistList->Add(fh1Xsec);
      fCommonHistList->Add(fh1Trials);
      fCommonHistList->Add(fh1PtHard);
      fCommonHistList->Add(fh1PtHardTrials);
      if(genJets) fCommonHistList->Add(fh1nGenJets);

      fFFHistosGen->AddToOutput(fCommonHistList);
      fFFHistosGenLeading->AddToOutput(fCommonHistList);
      fFFHistosGenLeadingTrack->AddToOutput(fCommonHistList);
    }
  }

  // Background
  if(fBckgMode){
    if(fFFMode){
      fFFBckgHisto0RecCuts->AddToOutput(fCommonHistList);
      fFFBckgHisto0RecLeading->AddToOutput(fCommonHistList);
      fFFBckgHisto1RecCuts->AddToOutput(fCommonHistList);
      fFFBckgHisto1RecLeading->AddToOutput(fCommonHistList);
      fFFBckgHisto2RecCuts->AddToOutput(fCommonHistList);
      fFFBckgHisto2RecLeading->AddToOutput(fCommonHistList);
      fFFBckgHisto3RecCuts->AddToOutput(fCommonHistList);
      fFFBckgHisto3RecLeading->AddToOutput(fCommonHistList);
      if(genJets && genTracks){
	fFFBckgHisto0Gen->AddToOutput(fCommonHistList);
	fFFBckgHisto0GenLeading->AddToOutput(fCommonHistList);
	fFFBckgHisto1Gen->AddToOutput(fCommonHistList);
	fFFBckgHisto1GenLeading->AddToOutput(fCommonHistList);
	fFFBckgHisto2Gen->AddToOutput(fCommonHistList);
	fFFBckgHisto2GenLeading->AddToOutput(fCommonHistList);
	fFFBckgHisto3Gen->AddToOutput(fCommonHistList);
	fFFBckgHisto3GenLeading->AddToOutput(fCommonHistList);
      }
    }

    if(fQAMode&1){
      fQABckgHisto0RecCuts->AddToOutput(fCommonHistList);
      fQABckgHisto1RecCuts->AddToOutput(fCommonHistList);
      fQABckgHisto2RecCuts->AddToOutput(fCommonHistList);
      fQABckgHisto3RecCuts->AddToOutput(fCommonHistList);
      if(genJets && genTracks){
	fQABckgHisto0Gen->AddToOutput(fCommonHistList);
	fQABckgHisto1Gen->AddToOutput(fCommonHistList);
	fQABckgHisto2Gen->AddToOutput(fCommonHistList);
	fQABckgHisto3Gen->AddToOutput(fCommonHistList);
      }
    }

    if(fBckgType[0]==kBckgOutLJ || fBckgType[1]==kBckgOutLJ || fBckgType[2]==kBckgOutLJ || fBckgType[3]==kBckgOutLJ)
      fCommonHistList->Add(fh1OutLeadingMult);
    if(fBckgType[0]==kBckgOutLJStat || fBckgType[1]==kBckgOutLJStat || fBckgType[2]==kBckgOutLJStat || fBckgType[3]==kBckgOutLJStat)
      fCommonHistList->Add(fh1OutLeadingStatMult);
    if(fBckgType[0]==kBckgPerp || fBckgType[1]==kBckgPerp || fBckgType[2]==kBckgPerp || fBckgType[3]==kBckgPerp )
      fCommonHistList->Add(fh1PerpMult);
    if(fBckgType[0]==kBckgASide || fBckgType[1]==kBckgASide || fBckgType[2]==kBckgASide || fBckgType[3]==kBckgASide)
      fCommonHistList->Add(fh1ASideMult);
    if(fBckgType[0]==kBckgASideWindow || fBckgType[1]==kBckgASideWindow || fBckgType[2]==kBckgASideWindow || fBckgType[3]==kBckgASideWindow)
      fCommonHistList->Add(fh1ASideWindowMult);
    if(fBckgType[0]==kBckgPerpWindow || fBckgType[1]==kBckgPerpWindow || fBckgType[2]==kBckgPerpWindow || fBckgType[3]==kBckgPerpWindow)
      fCommonHistList->Add(fh1PerpWindowMult);
    if(fBckgType[0]==kBckgOut2J || fBckgType[1]==kBckgOut2J || fBckgType[2]==kBckgOut2J || fBckgType[3]==kBckgOut2J)
      fCommonHistList->Add(fh1Out2JetsMult);
    if(fBckgType[0]==kBckgOut3J || fBckgType[1]==kBckgOut3J || fBckgType[2]==kBckgOut3J || fBckgType[3]==kBckgOut3J)
      fCommonHistList->Add(fh1Out3JetsMult);
    if(fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters)
      fCommonHistList->Add(fh1MedianClustersMult);
    if(fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)
      fCommonHistList->Add(fh1OutClustersMult);
  }

  // QA  
  if(fQAMode){
    if(fQAMode&1){ // track QA
      fQATrackHistosRec->AddToOutput(fCommonHistList);
      fQATrackHistosRecCuts->AddToOutput(fCommonHistList);
      if(genTracks) fQATrackHistosGen->AddToOutput(fCommonHistList);
    }

    if(fQAMode&2){ // jet QA
      fQAJetHistosRec->AddToOutput(fCommonHistList);
      fQAJetHistosRecCuts->AddToOutput(fCommonHistList);
      fQAJetHistosRecCutsLeading->AddToOutput(fCommonHistList);
      if(recJetsEff && fEffMode) fQAJetHistosRecEffLeading->AddToOutput(fCommonHistList); 
      if(genJets){
	fQAJetHistosGen->AddToOutput(fCommonHistList);
	fQAJetHistosGenLeading->AddToOutput(fCommonHistList);
      }
    }
  }

  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters ||  fBckgType[3]==kBckgClusters || 
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
      fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)) {
    fCommonHistList->Add(fh1nRecBckgJetsCuts);
    if(genJets) fCommonHistList->Add(fh1nGenBckgJets);
  }
    
  // phi correlation
  if(fPhiCorrMode){
    fPhiCorrHistosJetArea->AddToOutput(fCommonHistList);
    fPhiCorrHistosTransverseArea->AddToOutput(fCommonHistList);
    fPhiCorrHistosAwayArea->AddToOutput(fCommonHistList);
  }

  // intra-jet
  if(fIJMode){
    fIJHistosRecCuts->AddToOutput(fCommonHistList);
    fIJHistosRecLeading->AddToOutput(fCommonHistList);
    fIJHistosRecLeadingTrack->AddToOutput(fCommonHistList);
      
    // Background
    if(fBckgMode){
      fIJBckgHisto0RecCuts->AddToOutput(fCommonHistList);
      fIJBckgHisto0RecLeading->AddToOutput(fCommonHistList);
      fIJBckgHisto1RecCuts->AddToOutput(fCommonHistList);
      fIJBckgHisto1RecLeading->AddToOutput(fCommonHistList);
      fIJBckgHisto2RecCuts->AddToOutput(fCommonHistList);
      fIJBckgHisto2RecLeading->AddToOutput(fCommonHistList);
    }
      
    if(genJets && genTracks){
      fIJHistosGen->AddToOutput(fCommonHistList);
      fIJHistosGenLeading->AddToOutput(fCommonHistList);
      fIJHistosGenLeadingTrack->AddToOutput(fCommonHistList);
      
      // Background
      if(fBckgMode){
	fIJBckgHisto0Gen->AddToOutput(fCommonHistList);
	fIJBckgHisto0GenLeading->AddToOutput(fCommonHistList);
	fIJBckgHisto1Gen->AddToOutput(fCommonHistList);
	fIJBckgHisto1GenLeading->AddToOutput(fCommonHistList);
	fIJBckgHisto2Gen->AddToOutput(fCommonHistList);
	fIJBckgHisto2GenLeading->AddToOutput(fCommonHistList);
      }
    } // end: gen
  } // end: intra-jet
  
  if(fDJMode){
    if(fDJMode&1){
      fFFDiJetHistosRecCuts->AddToOutput(fCommonHistList);
      fFFDiJetHistosRecLeading->AddToOutput(fCommonHistList);
      fFFDiJetHistosRecLeadingTrack->AddToOutput(fCommonHistList);
      if(genJets && genTracks){
	fFFDiJetHistosGen->AddToOutput(fCommonHistList);
        fFFDiJetHistosGenLeading->AddToOutput(fCommonHistList);
        fFFDiJetHistosGenLeadingTrack->AddToOutput(fCommonHistList);
      }
    } // end: di-jet
    if(fDJMode&2){
      fQADiJetHistosRecCuts->AddToOutput(fCommonHistList);
      if(genJets && genTracks){
	fQADiJetHistosGen->AddToOutput(fCommonHistList);
      }
    } // end: di-jet QA
  } // end: di-jet
  
  if(fEffMode && recJetsEff && genTracks){
    if(fQAMode&1){
      fQATrackHistosRecEffGen->AddToOutput(fCommonHistList);
      fQATrackHistosRecEffRec->AddToOutput(fCommonHistList);
      fCommonHistList->Add(fhnResponseSinglePt);
    }
    if(fFFMode){
      fFFHistosRecEffGen->AddToOutput(fCommonHistList);
      fFFHistosRecEffRec->AddToOutput(fCommonHistList);
      fCommonHistList->Add(fhnResponseJetTrackPt);
      fCommonHistList->Add(fhnResponseJetZ);
      fCommonHistList->Add(fhnResponseJetXi);
    }
    fCommonHistList->Add(fh1nRecEffJets);
    fCommonHistList->Add(fh2PtRecVsGenPrim); 
  }
  

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fCommonHistList->GetEntries(); ++i){
    TH1 *h1 = dynamic_cast<TH1*>(fCommonHistList->At(i));
    if (h1) h1->Sumw2();
    else{
      THnSparse *hnSparse = dynamic_cast<THnSparse*>(fCommonHistList->At(i));
      if(hnSparse) hnSparse->Sumw2();
    }
  }
  
  TH1::AddDirectory(oldStatus);
}

//_______________________________________________
void AliAnalysisTaskFragmentationFunction::Init()
{
  // Initialization
  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::Init()");

}

//_____________________________________________________________
void AliAnalysisTaskFragmentationFunction::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  if(fDebug > 1) Printf("AliAnalysisTaskFragmentationFunction::UserExec()");
	

  if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);
  // Trigger selection
 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!(inputHandler->IsEventSelected() & AliVEvent::kMB)){
    if(inputHandler->InheritsFrom("AliESDInputHandler") && fUsePhysicsSelection){ // PhysicsSelection only with ESD input
      fh1EvtSelection->Fill(1.);
      if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
      PostData(1, fCommonHistList);
      return;
    }
  }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD){
    if(fDebug>3) Printf("%s:%d ESDEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  fMCEvent = MCEvent();
  if(!fMCEvent){
    if(fDebug>3) Printf("%s:%d MCEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  // get AOD event from input/ouput
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) {
    fAOD  =  ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD  = ((AliAODHandler*)handler)->GetAOD();
      if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
    }
  }
  
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return;
  }
  
  // event selection **************************************************
  // *** event class ***
  Double_t centPercent = -1;
  if(fEventClass>0){
    Int_t cl = 0;
    if(handler->InheritsFrom("AliAODInputHandler")){ 
      // since it is not supported by the helper task define own classes
      centPercent = fAOD->GetHeader()->GetCentrality();
      cl = 1;
      if(centPercent>10) cl = 2;
      if(centPercent>30) cl = 3;
      if(centPercent>50) cl = 4;
    }
    else {
      cl = AliAnalysisHelperJetTasks::EventClass();
    }
    
    if(cl!=fEventClass){
      // event not in selected event class, reject event
      if (fDebug > 1) Printf("%s:%d event not in selected event class: event REJECTED ...",(char*)__FILE__,__LINE__);
      fh1EvtSelection->Fill(2.);
      PostData(1, fCommonHistList);
      return;
    }
  }

  // *** vertex cut ***
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(!nTracksPrim){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return;
  }
  
  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>fMaxVertexZ){
    if (fDebug > 1) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ()); 
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return; 
  }
  
  TString primVtxName(primVtx->GetName());

  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(5.);
    PostData(1, fCommonHistList);
    return;
  }

  if (fDebug > 1) Printf("%s:%d event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fh1EvtSelection->Fill(0.);
  fh1EvtCent->Fill(centPercent);


  //___ get MC information __________________________________________________________________

  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  
  if(fMCEvent){
     AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    if(genHeader){
     AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
     AliGenHijingEventHeader*  hijingGenHeader = 0x0;

     if(pythiaGenHeader){
	 if(fDebug>3) Printf("%s:%d pythiaGenHeader found", (char*)__FILE__,__LINE__);
	 nTrials = pythiaGenHeader->Trials();
	 ptHard  = pythiaGenHeader->GetPtHard();

	 fh1PtHard->Fill(ptHard);
	 fh1PtHardTrials->Fill(ptHard,nTrials);


     } else { // no pythia, hijing?

	 if(fDebug>3) Printf("%s:%d no pythiaGenHeader found", (char*)__FILE__,__LINE__);

         hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
         if(!hijingGenHeader){
            Printf("%s:%d no pythiaGenHeader or hjingGenHeader found", (char*)__FILE__,__LINE__);
         } else {
	    if(fDebug>3) Printf("%s:%d hijingGenHeader found", (char*)__FILE__,__LINE__);
	 }
     }

     fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
    }
  }
  
  //___ fetch jets __________________________________________________________________________
 
  Int_t nJ = GetListOfJets(fJetsRec, kJetsRec);
  Int_t nRecJets = 0;
  if(nJ>=0) nRecJets = fJetsRec->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec jets: %d %d",(char*)__FILE__,__LINE__,nJ,nRecJets);
  if(nJ != nRecJets) Printf("%s:%d Mismatch Selected Rec Jets: %d %d",(char*)__FILE__,__LINE__,nJ,nRecJets);

  Int_t nJCuts = GetListOfJets(fJetsRecCuts, kJetsRecAcceptance);
  Int_t nRecJetsCuts = 0;
  if(nJCuts>=0) nRecJetsCuts = fJetsRecCuts->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  if(nRecJetsCuts != nJCuts) Printf("%s:%d Mismatch selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  fh1nRecJetsCuts->Fill(nRecJetsCuts);

  if(fJetTypeGen==kJetsKine || fJetTypeGen == kJetsKineAcceptance) fJetsGen->SetOwner(kTRUE); // kine aod jets allocated on heap, delete them with TList::Clear() 
  Int_t nJGen  = GetListOfJets(fJetsGen, fJetTypeGen);
  Int_t nGenJets = 0;
  if(nJGen>=0) nGenJets = fJetsGen->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
  if(nJGen != nGenJets) Printf("%s:%d Mismatch selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
  fh1nGenJets->Fill(nGenJets);


  if(fJetTypeRecEff==kJetsKine || fJetTypeRecEff == kJetsKineAcceptance) fJetsRecEff->SetOwner(kTRUE); // kine aod jets allocated on heap, delete them with TList::Clear() 
  Int_t nJRecEff  = GetListOfJets(fJetsRecEff, fJetTypeRecEff);
  Int_t nRecEffJets = 0;
  if(nJRecEff>=0) nRecEffJets = fJetsRecEff->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected RecEff jets: %d %d",(char*)__FILE__,__LINE__,nJRecEff,nRecEffJets);
  if(nJRecEff != nRecEffJets) Printf("%s:%d Mismatch selected RecEff jets: %d %d",(char*)__FILE__,__LINE__,nJRecEff,nRecEffJets);
  fh1nRecEffJets->Fill(nRecEffJets);

  //____ fetch background jets ___________________________________________________
  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
      fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)){

    Int_t nBJ = GetListOfBckgJets(fBckgJetsRec, kJetsRec);
    Int_t nRecBckgJets = 0;
    if(nBJ>=0) nRecBckgJets = fBckgJetsRec->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);
    if(nBJ != nRecBckgJets) Printf("%s:%d Mismatch Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);

    Int_t nBJCuts = GetListOfBckgJets(fBckgJetsRecCuts, kJetsRecAcceptance);
    Int_t nRecBckgJetsCuts = 0;
    if(nBJCuts>=0) nRecBckgJetsCuts = fBckgJetsRecCuts->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Rec background jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
    if(nRecBckgJetsCuts != nBJCuts) Printf("%s:%d Mismatch selected Rec background jets after cuts: %d %d",(char*)__FILE__,__LINE__,nBJCuts,nRecBckgJetsCuts);
    fh1nRecBckgJetsCuts->Fill(nRecBckgJetsCuts);

    if(0){ // protection OB - not yet implemented 
      if(fJetTypeGen==kJetsKine || fJetTypeGen == kJetsKineAcceptance) fBckgJetsGen->SetOwner(kTRUE); // kine aod jets allocated on heap, delete them with TList::Clear()
      Int_t nBJGen  = GetListOfBckgJets(fBckgJetsGen, fJetTypeGen);
      Int_t nGenBckgJets = 0;
      if(nBJGen>=0) nGenBckgJets = fBckgJetsGen->GetEntries();
      if(fDebug>2)Printf("%s:%d Selected Gen background jets: %d %d",(char*)__FILE__,__LINE__,nBJGen,nGenBckgJets);
      if(nBJGen != nGenBckgJets) Printf("%s:%d Mismatch selected Gen background jets: %d %d",(char*)__FILE__,__LINE__,nBJGen,nGenBckgJets);
      fh1nGenBckgJets->Fill(nGenBckgJets);
    }
  }


  //____ fetch particles __________________________________________________________
  
  Int_t nT = GetListOfTracks(fTracksRec, kTrackAOD);
  Int_t nRecPart = 0;
  if(nT>=0) nRecPart = fTracksRec->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,nRecPart);
  if(nRecPart != nT) Printf("%s:%d Mismatch selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,nRecPart);
  

  Int_t nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);
  Int_t nRecPartCuts = 0;
  if(nTCuts>=0) nRecPartCuts = fTracksRecCuts->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,nRecPartCuts);
  if(nRecPartCuts != nTCuts) Printf("%s:%d Mismatch selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,nRecPartCuts);
  fh1EvtMult->Fill(nRecPartCuts);



  Int_t nTGen = GetListOfTracks(fTracksGen,fTrackTypeGen);
  Int_t nGenPart = 0;
  if(nTGen>=0) nGenPart = fTracksGen->GetEntries();
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
  if(nGenPart != nTGen) Printf("%s:%d Mismatch selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
  
  
  //____ analysis, fill histos ___________________________________________________
  
  if(fQAMode){
    // loop over tracks
    if(fQAMode&1){
      for(Int_t it=0; it<nRecPart; ++it){
	AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRec->At(it));
	if(part)fQATrackHistosRec->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
      }
      for(Int_t it=0; it<nRecPartCuts; ++it){
	AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(it));
	if(part)fQATrackHistosRecCuts->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
      }
      for(Int_t it=0; it<nGenPart; ++it){
	AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksGen->At(it));
	if(part)fQATrackHistosGen->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
      }
    }

    // loop over jets
    
    if(fQAMode&2){
      for(Int_t ij=0; ij<nRecJets; ++ij){
	AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRec->At(ij));
	if(jet)fQAJetHistosRec->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      }
    }
  }

  if(fQAMode || fFFMode || fIJMode || fPhiCorrMode){
    for(Int_t ij=0; ij<nRecJetsCuts; ++ij){
      
      AliAODJet* jet = (AliAODJet*)(fJetsRecCuts->At(ij));
      if(fQAMode&2) fQAJetHistosRecCuts->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      
      if(ij==0){ // leading jet
	
	if(fQAMode&2) fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
	
	TList* jettracklist = new TList();
	Double_t sumPt      = 0.;
	Float_t leadTrackPt = 0.;
	TLorentzVector* leadTrackV = new TLorentzVector();
	
	if(GetFFRadius()<=0){
	  GetJetTracksTrackrefs(jettracklist, jet);
	} else {
	  GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt);
	}
	
	for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	  
	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	  if(!trackVP)continue;
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	  
	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  if(fFFMode) fFFHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt);
	  if(fIJMode) fIJHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	  
	  if(it==0){ // leading track 
	    leadTrackPt = trackPt;
	    leadTrackV->SetPxPyPzE(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	    
	    if(fFFMode) fFFHistosRecLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE);
	    if(fIJMode) fIJHistosRecLeadingTrack->FillIntraJet( leadTrackV, jet->MomentumVector() );
	  }
	  if(fFFMode) fFFHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt);
	  if(fIJMode) fIJHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  
	  delete trackV;
	}
	
	// ff and ij for background study
	if(fBckgMode){
	  if(fBckgType[0]!=-1)
	    FillBckgHistos(fBckgType[0], fTracksRecCuts, fJetsRecCuts, jet, leadTrackPt, leadTrackV, 
			   fFFBckgHisto0RecCuts, fFFBckgHisto0RecLeading,
			   fIJBckgHisto0RecCuts, fIJBckgHisto0RecLeading,
			   fQABckgHisto0RecCuts);
	  if(fBckgType[1]!=-1)
	    FillBckgHistos(fBckgType[1], fTracksRecCuts, fJetsRecCuts, jet, leadTrackPt, leadTrackV,
			   fFFBckgHisto1RecCuts, fFFBckgHisto1RecLeading,
			   fIJBckgHisto1RecCuts, fIJBckgHisto1RecLeading,
			   fQABckgHisto1RecCuts);
	  if(fBckgType[2]!=-1)
	    FillBckgHistos(fBckgType[2], fTracksRecCuts, fJetsRecCuts, jet, leadTrackPt, leadTrackV,
			   fFFBckgHisto2RecCuts, fFFBckgHisto2RecLeading,
			   fIJBckgHisto2RecCuts, fIJBckgHisto2RecLeading,
			   fQABckgHisto2RecCuts);
	  if(fBckgType[3]!=-1)
	    FillBckgHistos(fBckgType[3], fTracksRecCuts, fJetsRecCuts, jet, leadTrackPt, leadTrackV,
			   fFFBckgHisto3RecCuts, fFFBckgHisto3RecLeading,
			   fIJBckgHisto3RecCuts, fIJBckgHisto3RecLeading,
			   fQABckgHisto3RecCuts);
	} // end if(fBckgMode)
	
	delete leadTrackV;
	delete jettracklist;

	// phi correlation
	if(fPhiCorrMode){
	  for(Int_t it=0; it<nRecPartCuts; ++it){
	    AliVParticle *part = (AliVParticle*)(fTracksRecCuts->At(it));
	    
	    Float_t partEta = part->Eta();
	    Float_t partPhi = part->Phi();
	    Float_t partPt  = part->Pt();
	    
	    fPhiCorrHistosJetArea->FillTrackQA( partEta, 
						TVector2::Phi_mpi_pi( jet->Phi() - partPhi ),
						partPt,
						kTRUE);
	    
	    fPhiCorrHistosTransverseArea->FillTrackQA( partEta, 
						       TVector2::Phi_mpi_pi( jet->Phi() - partPhi + TMath::Pi()/2),
						       partPt,
						       kTRUE);
	    
	    fPhiCorrHistosAwayArea->FillTrackQA( partEta, 
						 TVector2::Phi_mpi_pi( jet->Phi() - partPhi + TMath::Pi()),
						 partPt,
						 kTRUE);
	  }
	} // end: phi-correlation
	
      } // end: leading jet
    } // end: rec. jets after cuts
    
    // generated jets
    for(Int_t ij=0; ij<nGenJets; ++ij){
      
      AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsGen->At(ij));
      if(!jet)continue;
      if(fQAMode&2) fQAJetHistosGen->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      
      if(ij==0){ // leading jet
	
	if(fQAMode&2) fQAJetHistosGenLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
	
	TList* jettracklist = new TList();
	Double_t sumPt      = 0.;
	Float_t leadTrackPt = 0.;
	TLorentzVector* leadTrackV = new TLorentzVector();
	
	if(GetFFRadius()<=0){
	  GetJetTracksTrackrefs(jettracklist, jet);
	} else {
	  GetJetTracksPointing(fTracksGen, jettracklist, jet, GetFFRadius(), sumPt);
	}
	
	for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	  
	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	  
	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  if(fFFMode) fFFHistosGen->FillFF( trackPt, jetPt, incrementJetPt );
	  if(fIJMode) fIJHistosGen->FillIntraJet( trackV, jet->MomentumVector() );
	  
	  if(it==0){ // leading track
	    leadTrackPt = trackPt;
	    leadTrackV->SetPxPyPzE(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	    
	    if(fFFMode) fFFHistosGenLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE );	  
	    if(fIJMode) fIJHistosGenLeadingTrack->FillIntraJet( leadTrackV, jet->MomentumVector() );
	  }
	  if(fFFMode) fFFHistosGenLeading->FillFF( trackPt, leadTrackPt, incrementJetPt );
	  if(fIJMode) fIJHistosGenLeading->FillIntraJet( trackV, leadTrackV );
	  
	  delete trackV;
	}
	
	delete leadTrackV;
	delete jettracklist;
      }
    }
  } // end: QA, FF and intra-jet

  //_______ DiJet part _____________________________________________________
  if(fDJMode){
    if (nRecJetsCuts > 1) 
      {
	AliAODJet* jet1 = (AliAODJet*)(fJetsRecCuts->At(0));
	AliAODJet* jet2 = (AliAODJet*)(fJetsRecCuts->At(1));
	
	// DiJet deltaphi calculation
	Double_t phi1 = TVector2::Phi_0_2pi(jet1->Phi());
	Double_t phi2 = TVector2::Phi_0_2pi(jet2->Phi());
	Double_t deltaPhi = TMath::Abs(phi1-phi2); 
	if (deltaPhi > TMath::Pi() && deltaPhi < 2*TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	
	// DiJet CDF cut calculation
	Double_t et1     = TMath::Abs(jet1->E()*TMath::Sin(jet1->Theta()));
	Double_t et2     = TMath::Abs(jet2->E()*TMath::Sin(jet2->Theta()));
	Double_t sumEt   = et1 + et2;
	Double_t normEt1PlusEt2   = TMath::Sqrt(et1*et1+et2*et2+2*et1*et2*TMath::Cos(deltaPhi));
	Double_t ratio = (Double_t)(normEt1PlusEt2/sumEt);
	
	// DiJet events selection
	Bool_t positionCut       = 0;
	Bool_t positionEnergyCut = 0;
	Bool_t cdfCut            = 0; 
	
	// Position cut :
	if (deltaPhi > fDiJetDeltaPhiCut) positionCut = 1;
	// Position-Energy cut :
	if ((deltaPhi > fDiJetDeltaPhiCut) && ((jet2->Pt()) >= fDiJetPtFractionCut*(jet1->Pt()))) positionEnergyCut = 1;
	// CDF cut :
	if (ratio < fDiJetCDFCut) cdfCut = 1;
	
	Int_t go = 0;
	
	if (fDiJetCut == 1 && positionCut == 1) go = 1;
	if (fDiJetCut == 2 && positionEnergyCut == 1) go = 1;
	if (fDiJetCut == 3 && cdfCut == 1) go = 1;
	
	if (go)
	  {
	    Double_t deltaEta      = TMath::Abs(jet1->Eta()-jet2->Eta());
	    Double_t deltaPt       = TMath::Abs(jet1->Pt()-jet2->Pt());
	    Double_t inbal         = (jet1->Pt()-jet2->Pt())/(jet1->Pt()+jet2->Pt());
	    Double_t meanEt        = (Double_t)((et1+et2)/2.);
	    Double_t invariantMass = (Double_t)InvMass(jet1,jet2);
	    
	    Double_t  jetBin = GetDiJetBin(invariantMass, jet1->Pt(), meanEt, fDiJetKindBins);
	    
	    if (jetBin > 0)
	      {
		if(fDJMode&2) fQADiJetHistosRecCuts->FillDiJetQA(invariantMass, deltaPhi, deltaEta, deltaPt, inbal, jetBin);
		
		if(fDJMode&1){
		  TList* jettracklist1 = new TList();
		  Double_t sumPt1      = 0.;
		  Float_t leadTrackPt1 = 0;
		  
		  TList* jettracklist2 = new TList();
		  Double_t sumPt2      = 0.;
		  Float_t leadTrackPt2 = 0;
		  
		  if(GetFFRadius()<=0)
		    {
		      GetJetTracksTrackrefs(jettracklist1, jet1);
		      GetJetTracksTrackrefs(jettracklist2, jet2);
		    }
		  else
		    {
		      GetJetTracksPointing(fTracksRecCuts, jettracklist1, jet1, GetFFRadius(), sumPt1);
		      GetJetTracksPointing(fTracksRecCuts, jettracklist2, jet2, GetFFRadius(), sumPt2);
		    }
		  
		  Int_t nTracks = jettracklist1->GetSize(); 
		  if (jettracklist1->GetSize() < jettracklist2->GetSize()) nTracks = jettracklist2->GetSize();
		  
		  for(Int_t it=0; it<nTracks; ++it)
		    {
		      if (it < jettracklist1->GetSize())
			{ 
			  AliVParticle *vp = (dynamic_cast<AliVParticle*> (jettracklist1->At(it)));
			  Float_t trackPt1 = (vp?vp->Pt():0);
			  Float_t jetPt1   = jet1->Pt();
			  
			  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
			  
			  fFFDiJetHistosRecCuts->FillDiJetFF(1, trackPt1, jetPt1, jetBin, incrementJetPt);
			  fFFDiJetHistosRecCuts->FillDiJetFF(0, trackPt1, jetPt1, jetBin, incrementJetPt);
			  
			  if (it == 0)
			    {
			      leadTrackPt1 = trackPt1;
			      
			      fFFDiJetHistosRecLeadingTrack->FillDiJetFF(1, leadTrackPt1, jetPt1, jetBin, kTRUE); 
			      fFFDiJetHistosRecLeadingTrack->FillDiJetFF(0, leadTrackPt1, jetPt1, jetBin, kTRUE); 
			    }
			  
			  fFFDiJetHistosRecLeading->FillDiJetFF(1, trackPt1, leadTrackPt1, jetBin, incrementJetPt); 
			  fFFDiJetHistosRecLeading->FillDiJetFF(0, trackPt1, leadTrackPt1, jetBin, incrementJetPt); 
			}
		      
		      if (it < jettracklist2->GetSize())
			{ 
			  Float_t trackPt2   = ((AliVParticle*)(jettracklist2->At(it)))->Pt();
			  Float_t jetPt2     = jet2->Pt();
			  
			  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
			  
			  fFFDiJetHistosRecCuts->FillDiJetFF(2, trackPt2, jetPt2, jetBin, incrementJetPt);
			  fFFDiJetHistosRecCuts->FillDiJetFF(0, trackPt2, jetPt2, jetBin, incrementJetPt);
			  
			  if (it == 0)
			    {
			      leadTrackPt2 = trackPt2;
			      
			      fFFDiJetHistosRecLeadingTrack->FillDiJetFF(2, leadTrackPt2, jetPt2, jetBin, kTRUE);
			      fFFDiJetHistosRecLeadingTrack->FillDiJetFF(0, leadTrackPt2, jetPt2, jetBin, kTRUE);
			    }
			  
			  fFFDiJetHistosRecLeading->FillDiJetFF(2, trackPt2, leadTrackPt2, jetBin, incrementJetPt);
			  fFFDiJetHistosRecLeading->FillDiJetFF(0, trackPt2, leadTrackPt2, jetBin, incrementJetPt);
			}
		    } // End loop on tracks
		  
		  delete jettracklist1;
		  delete jettracklist2;
		}
	      } // End if(jetBin > 0)
	    else { Printf("Jet bins for di-jet studies not set !");}
	  } // End if(go)
      } // End if(nRecJets > 1)
    
    if (nGenJets > 1)
      {
	AliAODJet* jet1 = dynamic_cast<AliAODJet*>(fJetsGen->At(0));
	AliAODJet* jet2 = dynamic_cast<AliAODJet*>(fJetsGen->At(1));
	if(jet1&&jet2){
	
	Double_t deltaPhi = 0;
	Double_t phi1 = TVector2::Phi_0_2pi(jet1->Phi());
	Double_t phi2 = TVector2::Phi_0_2pi(jet2->Phi());
	deltaPhi      = TMath::Abs(phi1-phi2); 
	if (deltaPhi > TMath::Pi() && deltaPhi < 2*TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	
	Double_t et1            = TMath::Abs(jet1->E()*TMath::Sin(jet1->Theta()));
	Double_t et2            = TMath::Abs(jet2->E()*TMath::Sin(jet2->Theta()));
	Double_t sumEt          = et1 + et2;
	Double_t normEt1PlusEt2 = TMath::Sqrt(et1*et1+et2*et2+2*et1*et2*TMath::Cos(deltaPhi));
	Double_t ratio          = (Double_t)(normEt1PlusEt2/sumEt);
	
	// DiJet events selection
	Bool_t positionCut       = 0;
	Bool_t positionEnergyCut = 0;
	Bool_t cdfCut            = 0; 
	
	// Position cut :
	if (deltaPhi > fDiJetDeltaPhiCut) positionCut = 1;
	// Position-Energy cut :
	if ((deltaPhi > fDiJetDeltaPhiCut) && ((jet2->Pt()) >= fDiJetPtFractionCut*(jet1->Pt()))) positionEnergyCut = 1;
	// CDF cut :
	if (ratio < fDiJetCDFCut) cdfCut = 1;    
	
	Int_t go = 0;
	
	if (fDiJetCut == 1 && positionCut == 1) go = 1;
	if (fDiJetCut == 2 && positionEnergyCut == 1) go = 1;
	if (fDiJetCut == 3 && cdfCut == 1) go = 1;
	
	if (go)
	  {
	    Double_t deltaEta      = TMath::Abs(jet1->Eta()-jet2->Eta());
	    Double_t deltaPt       = TMath::Abs(jet1->Pt()-jet2->Pt());
	    Double_t inbal         = (jet1->Pt()-jet2->Pt())/(jet1->Pt()+jet2->Pt());
	    Double_t meanEt        = (Double_t)((et1+et2)/2.);
	    Double_t invariantMass = (Double_t)InvMass(jet1,jet2);
	    
	    Double_t jetBin = GetDiJetBin(invariantMass, jet1->Pt(), meanEt, fDiJetKindBins);
	    
	    if(jetBin > 0)
	      {
		if(fDJMode&2) fQADiJetHistosGen->FillDiJetQA(invariantMass, deltaPhi, deltaEta, deltaPt, inbal, jetBin);
		
		if(fDJMode&1){
		  TList* jettracklist1 = new TList();
		  Double_t sumPt1 = 0.;
		  Float_t leadTrackPt1 = 0.;
		  
		  TList* jettracklist2 = new TList();
		  Double_t sumPt2 = 0.;
		  Float_t leadTrackPt2 = 0.;
		  
		  if(GetFFRadius()<=0)
		    {
		      GetJetTracksTrackrefs(jettracklist1, jet1);
		      GetJetTracksTrackrefs(jettracklist2, jet2);
		    }
		  else
		    {
		      GetJetTracksPointing(fTracksGen, jettracklist1, jet1, GetFFRadius(), sumPt1);
		      GetJetTracksPointing(fTracksGen, jettracklist2, jet2, GetFFRadius(), sumPt2);
		    }
		  
		  Int_t nTracks = jettracklist1->GetSize(); 
		  if (jettracklist1->GetSize() < jettracklist2->GetSize()) nTracks = jettracklist2->GetSize();
		  
		  for(Int_t it=0; it<nTracks; ++it)
		    {
		      if (it < jettracklist1->GetSize())
			{ 
			  Float_t trackPt1 = ((AliVParticle*)jettracklist1->At(it))->Pt();
			  Float_t jetPt1 = jet1->Pt();
			  
			  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
			  
			  fFFDiJetHistosGen->FillDiJetFF( 1, trackPt1, jetPt1, jetBin, incrementJetPt);
			  fFFDiJetHistosGen->FillDiJetFF( 0, trackPt1, jetPt1, jetBin, incrementJetPt);
			  
			  if(it==0)
			    { 
			      leadTrackPt1 = trackPt1;
			      
			      fFFDiJetHistosGenLeadingTrack->FillDiJetFF( 1, leadTrackPt1, jetPt1, jetBin, kTRUE);
			      fFFDiJetHistosGenLeadingTrack->FillDiJetFF( 0, leadTrackPt1, jetPt1, jetBin, kTRUE);
			    }
			  
			  fFFDiJetHistosGenLeading->FillDiJetFF( 1, trackPt1, leadTrackPt1, jetBin, incrementJetPt);
			  fFFDiJetHistosGenLeading->FillDiJetFF( 0, trackPt1, leadTrackPt1, jetBin, incrementJetPt);
			}
		      
		      if (it < jettracklist2->GetSize())
			{ 
			  Float_t trackPt2 = (dynamic_cast<AliVParticle*>(jettracklist2->At(it)))->Pt();
			  Float_t jetPt2 = jet2->Pt();
			  
			  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
			  
			  fFFDiJetHistosGen->FillDiJetFF( 2, trackPt2, jetPt2, jetBin, incrementJetPt);
			  fFFDiJetHistosGen->FillDiJetFF( 0, trackPt2, jetPt2, jetBin, incrementJetPt);
			  
			  if (it==0)
			    { 
			      leadTrackPt2 = trackPt2;
			      
			      fFFDiJetHistosGenLeadingTrack->FillDiJetFF( 2, leadTrackPt2, jetPt2, jetBin, kTRUE);
			      fFFDiJetHistosGenLeadingTrack->FillDiJetFF( 0, leadTrackPt2, jetPt2, jetBin, kTRUE);
			    }
			  
			  fFFDiJetHistosGenLeading->FillDiJetFF( 2, trackPt2, leadTrackPt2, jetBin, incrementJetPt);
			  fFFDiJetHistosGenLeading->FillDiJetFF( 0, trackPt2, leadTrackPt2, jetBin, incrementJetPt);
			}
		    } // End loop on tracks
		  
		  delete jettracklist1;
		  delete jettracklist2;
		}
	      } // End if(jetBin > 0)
	    else { Printf("Jet bins for di-jet studies not set !");}
	  } // End if (go)
	}// end if jet1 and jet2
      } // End if(nGenJets > 1)
  } // end: di-jet
  
  // ____ efficiency _______________________________

  if(fEffMode && (fJetTypeRecEff != kJetsUndef)){

    // arrays for generated particles: reconstructed AOD track index, isPrimary flag, are initialized in AssociateGenRec(...) function
    TArrayI indexAODTr; 
    TArrayS isGenPrim; 
  
    // array for reconcstructed AOD tracks: generated particle index, initialized in AssociateGenRec(...) function
    TArrayI indexMCTr; 
  
    Int_t  nTracksAODMCCharged = GetListOfTracks(fTracksAODMCCharged, kTrackAODMCCharged);
    if(fDebug>2)Printf("%s:%d selected AODMC tracks: %d ",(char*)__FILE__,__LINE__,nTracksAODMCCharged);
  
    Int_t  nTracksRecQualityCuts = GetListOfTracks(fTracksRecQualityCuts, kTrackAODQualityCuts);
    if(fDebug>2)Printf("%s:%d selected rec tracks quality after cuts, full acceptance/pt : %d ",(char*)__FILE__,__LINE__,nTracksRecQualityCuts);
  
    // associate gen and rec tracks, store indices in TArrays 
    AssociateGenRec(fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,indexMCTr,isGenPrim);
  
    // single track eff
    if(fQAMode&1) FillSingleTrackRecEffHisto(fQATrackHistosRecEffGen,fQATrackHistosRecEffRec,fTracksAODMCCharged,indexAODTr,isGenPrim);
    if(fQAMode&1) FillSingleTrackResponse(fhnResponseSinglePt,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim);


    // jet track eff
    
    Double_t sumPtGenLeadingJetRecEff = 0;
    Double_t sumPtRecLeadingJetRecEff = 0;
    
    for(Int_t ij=0; ij<nRecEffJets; ++ij){ 
    
      AliAODJet* jet = (AliAODJet*)(fJetsRecEff->At(ij));
    
      if(ij==0){ // leading jet
	
	TList* jettracklistGen = new TList();
	GetJetTracksPointing(fTracksGen, jettracklistGen, jet, GetFFRadius(), sumPtGenLeadingJetRecEff); // for efficiency: gen tracks from pointing with gen/rec jet
	
	TList* jettracklistRec = new TList();
	GetJetTracksPointing(fTracksRecCuts,jettracklistRec, jet, GetFFRadius(), sumPtRecLeadingJetRecEff); // bin efficiency in jet pt bins using rec tracks  
	
	Double_t jetEta   = jet->Eta();
	Double_t jetPhi   = TVector2::Phi_0_2pi(jet->Phi());
	
	if(fQAMode&2) fQAJetHistosRecEffLeading->FillJetQA( jetEta, jetPhi, sumPtGenLeadingJetRecEff ); 
	
	if(fFFMode) FillJetTrackRecEffHisto(fFFHistosRecEffGen,fFFHistosRecEffRec,sumPtGenLeadingJetRecEff,sumPtRecLeadingJetRecEff,
					    jettracklistGen,fTracksAODMCCharged,indexAODTr,isGenPrim,fUseRecEffRecJetPtBins); 
	
	if(fFFMode) FillJetTrackResponse(fhnResponseJetTrackPt,fhnResponseJetZ,fhnResponseJetXi,sumPtGenLeadingJetRecEff,sumPtRecLeadingJetRecEff,
 					 jettracklistGen,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,fUseResponseRecJetPtBins);

	delete jettracklistGen;
	delete jettracklistRec;
      }
    }
    

    // bckgr eff: complementary cones 

    if(0){
      
      for(Int_t ij=0; ij<nRecEffJets; ++ij){ 
	
	AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRecEff->At(ij));
	
	if(ij==0){ // leading jet
	  
	  TList* perpjettracklistGen = new TList();
	  Double_t sumPtGen = 0.; 
	  
	  GetTracksTiltedwrpJetAxis(TMath::Pi()/2.,fTracksGen, perpjettracklistGen, jet, GetFFRadius() , sumPtGen); // for efficiency: gen tracks perp to gen/rec jet
	  
	  // here could be your histos !!! 
	  // FillJetTrackRecEffHisto(fFFBckgrPerpHistosRecEffGen,fFFBckgrPerpHistosRecEffRec,sumPtGenLeadingJetRecEff,sumPtRecLeadingJetRecEff,perpjettracklistGen,
	  //	  	         fTracksAODMCCharged,indexAODTr,isGenPrim,fUseRecEffRecJetPtBins);
	  
	  delete perpjettracklistGen;
	}
      }
    }

    // bgr eff: outside N leading jets

    if(0){
      
      Int_t nCases = 1;
      
      TList* outjettracklistGen = new TList();
      Double_t sumPtGen = 0.;
      
      GetTracksOutOfNJets(nCases, fTracksGen, outjettracklistGen, fJetsRecEff, sumPtGen); // for efficiency: gen tracks outide n gen/rec jets 
      
      // here could be your histos !!! 
      // FillJetTrackRecEffHisto(fFFBckgrOutHistosRecEffGen,fFFBckgrOutHistosRecEffRec,sumPtGenLeadingJetRecEff,sumPtRecLeadingJetRecEff,
      //                         outjettracklistGen,fTracksAODMCCharged,indexAODTr,isGenPrim,fUseRecEffRecJetPtBins);
      
      delete outjettracklistGen;
    }
  }     

  //___________________
  
  fTracksRec->Clear();
  fTracksRecCuts->Clear();
  fTracksGen->Clear();
  fTracksAODMCCharged->Clear();
  fTracksRecQualityCuts->Clear();

  fJetsRec->Clear();
  fJetsRecCuts->Clear();
  fJetsGen->Clear();
  fJetsRecEff->Clear();

  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || 
      fBckgType[2]==kBckgClustersOutLeading || fBckgType[3]==kBckgClustersOutLeading)){

    fBckgJetsRec->Clear();
    fBckgJetsRecCuts->Clear();
    fBckgJetsGen->Clear();
  }

  //Post output data.
  PostData(1, fCommonHistList);
}

//________________________________________________________________________________________
Double_t AliAnalysisTaskFragmentationFunction::InvMass(const AliAODJet* jet1, const AliAODJet* jet2)
{
  // cald DiJet inv mass

  Double_t invMass = 0.;
  invMass = TMath::Sqrt(pow(jet1->E()+jet2->E(),2) - pow(jet1->Px()+jet2->Px(),2) - 
			pow(jet1->Py()+jet2->Py(),2) - pow(jet1->Pz()+jet2->Pz(),2));

  return invMass;

}

//________________________________________________________________________________________
Double_t AliAnalysisTaskFragmentationFunction::GetDiJetBin(Double_t invMass, Double_t leadingJetPt, Double_t EtMean, Int_t kindBins)
{
  // calc DiJet bin according to kindBins parameter

  Double_t jetBinOk = 0.;
  Double_t jetBin = 0.;

  Float_t stepInvMass = (fDiJetJetInvMassMax - fDiJetJetInvMassMin)/fDiJetNBinsJetInvMass;
  Float_t stepPt = (fDiJetJetPtMax - fDiJetJetPtMin)/fDiJetNBinsJetPt;

  if (kindBins == 1)
    {
      for(Int_t i=0; i<fDiJetNBinsJetInvMass; ++i)
	{
	  jetBin = fDiJetJetInvMassMin + i*stepInvMass + stepInvMass/2.;
	  if(((fDiJetJetInvMassMin+i*stepInvMass) <= invMass) &&
	     (fDiJetJetInvMassMin + (i+1)*stepInvMass) > invMass) {jetBinOk = jetBin; break;}
          else jetBinOk = -1.;
	}
    }
  else if (kindBins == 3)
    {
      for(Int_t i=0; i<fDiJetNBinsJetPt; ++i)
	{
	  jetBin = fDiJetJetPtMin + i*stepPt + stepPt/2.;
	  if(((fDiJetJetPtMin+i*stepPt) <= EtMean) &&
	     (fDiJetJetPtMin + (i+1)*stepPt) > EtMean) {jetBinOk = jetBin; break;}
          else jetBinOk = -1.;
	}
    }
  else if (kindBins == 2)
    {
      for(Int_t i=0; i<fDiJetNBinsJetPt; ++i)
	{
	  jetBin = fDiJetJetPtMin + i*stepPt + stepPt/2.;
	  if(((fDiJetJetPtMin+i*stepPt) <= leadingJetPt) &&
	     (fDiJetJetPtMin + (i+1)*stepPt) > leadingJetPt) {jetBinOk = jetBin; break;}
          else jetBinOk = -1.;
	}
    }
  else {Printf("WARNING: kindBins wrongly set ! Please make sure to call SetKindSlices() and set the kind parameter to 1, 2 or 3.\n");}

  return jetBinOk;

}


//______________________________________________________________
void AliAnalysisTaskFragmentationFunction::Terminate(Option_t *) 
{
  // terminated

  if(fDebug > 1) printf("AliAnalysisTaskFragmentationFunction::Terminate() \n");
}  

//_________________________________________________________________________________
Int_t AliAnalysisTaskFragmentationFunction::GetListOfTracks(TList *list, Int_t type)
{
  // fill list of tracks selected according to type

  if(fDebug > 2) Printf("%s:%d Selecting tracks with %d", (char*)__FILE__,__LINE__,type);
  
  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }

  if(!fAOD) return -1;

  if(!fAOD->GetTracks()) return 0;

  if(type==kTrackUndef) return 0;
  
  Int_t iCount = 0;
  if(type==kTrackAODCuts || type==kTrackAODQualityCuts || type==kTrackAOD){

    // all rec. tracks, esd filter mask, eta range

    
    for(Int_t it=0; it<fAOD->GetNumberOfTracks(); ++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      
      if(type == kTrackAODCuts || type==kTrackAODQualityCuts ){
	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))   continue;
	if(type == kTrackAODCuts){
	  if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax) continue;
	  if(tr->Phi() < fTrackPhiMin || tr->Phi() > fTrackPhiMax) continue;
	  if(tr->Pt()  < fTrackPtCut) continue;
	}
      }
      list->Add(tr);
      iCount++;
    }
  }
  else if (type==kTrackKineAll || type==kTrackKineCharged || type==kTrackKineChargedAcceptance){
    // kine particles, all or rather charged
    if(!fMCEvent) return iCount;
    
    for(Int_t it=0; it<fMCEvent->GetNumberOfTracks(); ++it){
      AliMCParticle* part = (AliMCParticle*) fMCEvent->GetTrack(it);
      
      if(type == kTrackKineCharged || type == kTrackKineChargedAcceptance){
	if(part->Charge()==0) continue;
	
	if(type == kTrackKineChargedAcceptance && 
	   (       part->Eta() < fTrackEtaMin
		|| part->Eta() > fTrackEtaMax
		|| part->Phi() < fTrackPhiMin
		|| part->Phi() > fTrackPhiMax 
		|| part->Pt()  < fTrackPtCut)) continue;
      }
      
      list->Add(part);
      iCount++;
    }
  }
  else if (type==kTrackAODMCCharged || type==kTrackAODMCAll || type==kTrackAODMCChargedAcceptance) {
    // MC particles (from AOD), physical primaries, all or rather charged or rather charged within acceptance
    if(!fAOD) return -1;
    
    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    
    for(int it=0; it<tca->GetEntriesFast(); ++it){
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
      if(!part)continue;
      if(!part->IsPhysicalPrimary())continue;
      
      if (type==kTrackAODMCCharged || type==kTrackAODMCChargedAcceptance){
	if(part->Charge()==0) continue;
	if(type==kTrackAODMCChargedAcceptance && 
	   (     part->Eta() > fTrackEtaMax
	      || part->Eta() < fTrackEtaMin
	      || part->Phi() > fTrackPhiMax
	      || part->Phi() < fTrackPhiMin
	      || part->Pt()  < fTrackPtCut)) continue;
      }
      
      list->Add(part);
      iCount++;
    }
  }
  
  list->Sort();
  return iCount;
  
}
// _______________________________________________________________________________
Int_t AliAnalysisTaskFragmentationFunction::GetListOfJets(TList *list, Int_t type)
{
  // fill list of jets selected according to type

  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }

  if(type == kJetsRec || type == kJetsRecAcceptance){ // reconstructed jets

    if(fBranchRecJets.Length()==0){
      Printf("%s:%d no rec jet branch specified", (char*)__FILE__,__LINE__);
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    TClonesArray *aodRecJets = 0; 
    if(fBranchRecJets.Length()) aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRecJets.Data()));
    if(!aodRecJets)             aodRecJets = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(fBranchRecJets.Data()));

    if(!aodRecJets){
      if(fBranchRecJets.Length()) Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fBranchRecJets.Data());
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    // Reorder jet pt and fill new temporary AliAODJet objects
    Int_t nRecJets = 0;
    
    for(Int_t ij=0; ij<aodRecJets->GetEntries(); ++ij){

      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ij));
      if(!tmp) continue;

      if( tmp->Pt() < fJetPtCut ) continue;
      if( type == kJetsRecAcceptance &&
	  (    tmp->Eta() < fJetEtaMin
	    || tmp->Eta() > fJetEtaMax
	    || tmp->Phi() < fJetPhiMin
	    || tmp->Phi() > fJetPhiMax )) continue;

      if(fBckgSubMethod && fBranchRecJets.Contains("B0") &&
	 fBranchRecJets.Contains("KT")) {

	AliAODJet *tmpJet = GetAODBckgSubJet(tmp, fBckgSubMethod);
	
	if(!tmpJet) continue;

	list->Add(tmpJet);

	nRecJets++;
      }
      else {
	list->Add(tmp);
	
	nRecJets++;
      }
    }
    
    list->Sort();
    
    return nRecJets;
  }
  else if(type == kJetsKine || type == kJetsKineAcceptance){
    
    // generated jets
    Int_t nGenJets = 0;
    
    if(!fMCEvent){
      if(fDebug>1) Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
      return 0;
    }
   
    AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenHijingEventHeader*  hijingGenHeader = 0x0;

    if(!pythiaGenHeader){
      hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
      
      if(!hijingGenHeader){
         Printf("%s:%d no pythiaGenHeader or hijingGenHeader found", (char*)__FILE__,__LINE__);
         return 0;
      }else{
         TLorentzVector mom[4];
         AliAODJet* jet[4];
         hijingGenHeader->GetJets(mom[0], mom[1], mom[2], mom[3]);

         for(Int_t i=0; i<2; ++i){
	    if(!mom[i].Pt()) continue;
            jet[i] = new AliAODJet(mom[i]);

            if( type == kJetsKineAcceptance &&
                (    jet[i]->Eta() < fJetEtaMin
                  || jet[i]->Eta() > fJetEtaMax
                  || jet[i]->Phi() < fJetPhiMin
                  || jet[i]->Phi() > fJetPhiMax )) continue;

	    list->Add(jet[i]);
	    nGenJets++;
	 }
	 list->Sort();
         return nGenJets;
      }
    }
    
    // fetch the pythia generated jets
    for(int ip=0; ip<pythiaGenHeader->NTriggerJets(); ++ip){
      
      Float_t p[4];
      AliAODJet *jet = new AliAODJet();
      pythiaGenHeader->TriggerJet(ip, p);
      jet->SetPxPyPzE(p[0], p[1], p[2], p[3]);

      if( type == kJetsKineAcceptance &&
          (    jet->Eta() < fJetEtaMin
            || jet->Eta() > fJetEtaMax
            || jet->Phi() < fJetPhiMin
            || jet->Phi() > fJetPhiMax )) continue;
      
	list->Add(jet);
	nGenJets++;
    }
    list->Sort();
    return nGenJets;
  }
  else if(type == kJetsGen || type == kJetsGenAcceptance ){

    if(fBranchGenJets.Length()==0){
      if(fDebug>1) Printf("%s:%d no gen jet branch specified", (char*)__FILE__,__LINE__);
      return 0;
    }
    
    TClonesArray *aodGenJets = 0;
    if(fBranchGenJets.Length()) aodGenJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchGenJets.Data()));
    if(!aodGenJets)             aodGenJets = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(fBranchGenJets.Data()));

    if(!aodGenJets){
      if(fDebug>0){
	if(fBranchGenJets.Length())         Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGenJets.Data());
      }
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    Int_t nGenJets = 0;
    
    for(Int_t ig=0; ig<aodGenJets->GetEntries(); ++ig){
	  
      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodGenJets->At(ig));
      if(!tmp) continue;
	  
      if( tmp->Pt() < fJetPtCut ) continue;
      if( type == kJetsGenAcceptance &&
	  (    tmp->Eta() < fJetEtaMin
	    || tmp->Eta() > fJetEtaMax
	    || tmp->Phi() < fJetPhiMin
	    || tmp->Phi() > fJetPhiMax )) continue;
      
	list->Add(tmp);
      	nGenJets++;
    }
    list->Sort();
    return nGenJets;
  } 
  else{
    if(fDebug>0)Printf("%s:%d no such type %d",(char*)__FILE__,__LINE__,type);
    return 0;
  }
}

// ___________________________________________________________________________________
Int_t AliAnalysisTaskFragmentationFunction::GetListOfBckgJets(TList *list, Int_t type)  
{
  // fill list of bgr clusters selected according to type

  if(type == kJetsRec || type == kJetsRecAcceptance){ // reconstructed jets

    if(fBranchRecBckgClusters.Length()==0){ 
      Printf("%s:%d no rec jet branch specified", (char*)__FILE__,__LINE__);
      if(fDebug>1)fAOD->Print();
      return 0;
    }
    
    TClonesArray *aodRecJets = 0; 
    if(fBranchRecBckgClusters.Length()) aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRecBckgClusters.Data()));
    if(!aodRecJets)                     aodRecJets = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(fBranchRecBckgClusters.Data()));
    
    if(!aodRecJets){
      if(fBranchRecBckgClusters.Length()) Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fBranchRecBckgClusters.Data());
      if(fDebug>1)fAOD->Print();
      return 0;
    }
    
    // Reorder jet pt and fill new temporary AliAODJet objects
    Int_t nRecJets = 0;
    
    for(Int_t ij=0; ij<aodRecJets->GetEntries(); ++ij){
      
      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ij));
      if(!tmp) continue;

      // if( tmp->Pt() < fJetPtCut ) continue; // no pt cut on bckg clusters !
      if( type == kJetsRecAcceptance &&
       	  (    tmp->Eta() < fJetEtaMin
	       || tmp->Eta() > fJetEtaMax
	       || tmp->Phi() < fJetPhiMin
	       || tmp->Phi() > fJetPhiMax )) continue;
            
      list->Add(tmp);
	
      nRecJets++;
      
    }
    
    list->Sort();
    
    return nRecJets;
  }

  //  /*
  // MC clusters still Under construction
  //  */

  return 0;
} 

// _________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(THnSparse* h,const Int_t dim, const char** labels)
{
  // Set properties of THnSparse 

  for(Int_t i=0; i<dim; i++){
    h->GetAxis(i)->SetTitle(labels[i]);
    h->GetAxis(i)->SetTitleColor(1);
  }
}

// __________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y)
{
  //Set properties of histos (x and y title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
}

// _________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y, const char* z)
{
  //Set properties of histos (x,y and z title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->SetZTitle(z);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->GetZaxis()->SetTitleColor(1);
}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetJetTracksPointing(TList* inputlist, TList* outputlist, const AliAODJet* jet, const Double_t radius,Double_t& sumPt)
{
  // fill list of tracks in cone around jet axis  

  sumPt = 0;

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){

      outputlist->Add(track);
      
      sumPt += track->Pt();
    }
  }
  
  outputlist->Sort();
}

// ___________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetJetTracksTrackrefs(TList* list, const AliAODJet* jet)
{
  // list of jet tracks from trackrefs
  
  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  for (Int_t itrack=0; itrack<nTracks; itrack++) {
    
    AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(itrack));
    if(!track){
      AliError("expected ref track not found ");
      continue;
    }
        
    list->Add(track);
  }
  
  list->Sort();
}

// _ ________________________________________________________________________________________________________________________________
void  AliAnalysisTaskFragmentationFunction::AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,TArrayS& isGenPrim)
{
  // associate generated and reconstructed tracks, fill TArrays of list indices


  Int_t nTracksRec  = tracksRec->GetSize();
  Int_t nTracksGen  = tracksAODMCCharged->GetSize();
  TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));

  if(!nTracksGen) return;
  if(!tca)        return;
  
  // set size
  indexAODTr.Set(nTracksGen);
  indexMCTr.Set(nTracksRec);
  isGenPrim.Set(nTracksGen);

  indexAODTr.Reset(-1);
  indexMCTr.Reset(-1);
  isGenPrim.Reset(0);

  // loop over reconstructed tracks, get generated track 

  for(Int_t iRec=0; iRec<nTracksRec; iRec++){ 
      
    AliAODTrack* rectrack = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)); 
    if(!rectrack)continue;
    Int_t label = TMath::Abs(rectrack->GetLabel());

    // find MC track in our list
    AliAODMCParticle* gentrack = dynamic_cast<AliAODMCParticle*> (tca->At(label));
   
    Int_t listIndex = -1;
    if(gentrack) listIndex = tracksAODMCCharged->IndexOf(gentrack);

    if(listIndex>=0){

      indexAODTr[listIndex] = iRec;
      indexMCTr[iRec]       = listIndex;
    }
  } 


  // define primary sample for reconstruction efficiency

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksAODMCCharged->At(iGen));
    if(!gentrack)continue;
    Int_t pdg = gentrack->GetPdgCode();    

    // 211 - pi, 2212 - proton, 321 - Kaon, 11 - electron, 13 - muon
    if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || 
       TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){
      
      isGenPrim[iGen] = kTRUE;

      Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

      if(iRec>=0){
	Float_t genPt = gentrack->Pt();
	AliAODTrack* vt = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)); 
	if(vt){
	  Float_t recPt = vt->Pt();
	  fh2PtRecVsGenPrim->Fill(genPt,recPt);
	}
      }
    }
  }
}

// _____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillSingleTrackRecEffHisto(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, 
								      const TArrayI& indexAODTr, const TArrayS& isGenPrim){

  // fill QA for single track reconstruction efficiency
  
  Int_t nTracksGen  = tracksGen->GetSize();

  if(!nTracksGen) return;

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    if(isGenPrim[iGen] != 1) continue; // select primaries

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksGen->At(iGen));
    if(!gentrack)continue;
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 

    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtCut) continue;

    trackQAGen->FillTrackQA(etaGen, phiGen, ptGen);

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 
    if(iRec>=0) trackQARec->FillTrackQA(etaGen, phiGen, ptGen);
  }
}

// ______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillJetTrackRecEffHisto(TObject* histGen, TObject* histRec, Double_t jetPtGen, Double_t jetPtRec, TList* jetTrackList, 
								   const TList* tracksGen, const TArrayI& indexAODTr, const TArrayS& isGenPrim, const Bool_t useRecJetPt)
{

  // fill objects for jet track reconstruction efficiency
  // arguments histGen/histRec can be of different type: AliFragFuncHistos*, AliFragFuncIntraJetHistos*, ...

  Int_t nTracksJet = jetTrackList->GetSize(); // list with AODMC tracks

  if(!nTracksJet) return; 

  Bool_t incrementJetPtGenFF = kTRUE; // needed in case we fill FFHistos
  Bool_t incrementJetPtRecFF = kTRUE; // needed in case we fill FFHistos

  for(Int_t iTr=0; iTr<nTracksJet; iTr++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (jetTrackList->At(iTr));
    if(!gentrack)continue;
    // find jet track in gen tracks list
    Int_t iGen = tracksGen->IndexOf(gentrack); 

    if(iGen<0){
      if(fDebug>0) Printf("%s:%d gen jet track not found ",(char*)__FILE__,__LINE__);
      continue;
    }

    if(isGenPrim[iGen] != 1) continue; // select primaries
    
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 

    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtCut) continue;

    Double_t z = ptGen / jetPtGen;
    Double_t xi = 0;
    if(z>0) xi = TMath::Log(1/z);

    Int_t iRec   = indexAODTr[iGen]; // can be -1 if no good reconstructed track 
    Bool_t isRec = (iRec>=0) ? kTRUE : kFALSE; 

    if(dynamic_cast<AliFragFuncHistos*>(histGen) && dynamic_cast<AliFragFuncHistos*>(histRec)){
      
      // after checking can afford normal cast
      AliFragFuncHistos* effFFhistGen = (AliFragFuncHistos*)(histGen); 
      AliFragFuncHistos* effFFhistRec = (AliFragFuncHistos*)(histRec); 

      if(useRecJetPt) effFFhistGen->FillFF( ptGen, jetPtRec, incrementJetPtGenFF );
      else            effFFhistGen->FillFF( ptGen, jetPtGen, incrementJetPtGenFF );

      incrementJetPtGenFF = kFALSE;

      if(isRec){
	if(useRecJetPt) effFFhistRec->FillFF( ptGen, jetPtRec, incrementJetPtRecFF );
	else            effFFhistRec->FillFF( ptGen, jetPtGen, incrementJetPtRecFF );
	
	incrementJetPtRecFF = kFALSE;
      }
    }
    else if(dynamic_cast<AliFragFuncIntraJetHistos*>(histGen) && dynamic_cast<AliFragFuncIntraJetHistos*>(histRec)){
      
      // eff for IJ histos ...
      
    }
  }
}


// _____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillSingleTrackResponse(THnSparse* hnResponse, TList* tracksGen,  TList* tracksRec,
								   const TArrayI& indexAODTr, const TArrayS& isGenPrim)
{
  // fill response matrix for single tracks 
  

  Int_t nTracksGen  = tracksGen->GetSize();

  if(!nTracksGen) return;

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    if(isGenPrim[iGen] != 1) continue; // select primaries

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksGen->At(iGen));
    if(!gentrack)continue;
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 
    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtCut) continue;

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 
    if(iRec>=0){ 
      AliAODTrack* rectrack = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)); 
      if(!rectrack)continue;
      Double_t ptRec = rectrack->Pt();
      
      Double_t entries[2] = {ptRec,ptGen}; // AliCorrFW convention: gen vs rec
      hnResponse->Fill(entries);
    }
  }
}


// ______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillJetTrackResponse(THnSparse* hnResponsePt, THnSparse* hnResponseZ, THnSparse* hnResponseXi, 
								Double_t jetPtGen, Double_t jetPtRec, TList* jetTrackList, 
								const TList* tracksGen, TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isGenPrim,const Bool_t useRecJetPt)
{
  // fill response matrix for tracks in jets
  
  Int_t nTracksJet = jetTrackList->GetSize(); // list with AODMC tracks

  if(!nTracksJet) return; 

  for(Int_t iTr=0; iTr<nTracksJet; iTr++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (jetTrackList->At(iTr));
    if(!gentrack)continue;
    // find jet track in gen tracks list
    Int_t iGen = tracksGen->IndexOf(gentrack); 

    if(iGen<0){
      if(fDebug>0) Printf("%s:%d gen jet track not found ",(char*)__FILE__,__LINE__);
      continue;
    }

    if(isGenPrim[iGen] != 1) continue; // select primaries
    
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 

    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtCut) continue;

    Double_t zGen = ptGen / jetPtRec;
    Double_t xiGen = 0;
    if(zGen>0) xiGen = TMath::Log(1/zGen);

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

    if(iRec>=0){

      AliAODTrack* rectrack = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)); 
      if(rectrack){
	Double_t ptRec = rectrack->Pt();
	
	Double_t zRec    = ptRec / jetPtRec;
	Double_t xiRec   = 0;
	if(zRec>0) xiRec = TMath::Log(1/zRec);
      
	Double_t jetPt = useRecJetPt ? jetPtRec : jetPtGen;
	
	Double_t entriesPt[3] = {jetPt,ptRec,ptGen}; // AliCorrFW convention: gen vs rec
	Double_t entriesZ[3]  = {jetPt,zRec,zGen}; 
	Double_t entriesXi[3] = {jetPt,xiRec,xiGen}; 
	
	hnResponsePt->Fill(entriesPt);
	hnResponseZ->Fill(entriesZ);
	hnResponseXi->Fill(entriesXi);
      }
    } 
  }
}

// _____________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, AliAODJet* jet, Double_t radius,Double_t& sumPt)
{
  // List of tracks in cone perpendicular to the jet azimuthal direction

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);

  TVector3 jet3mom(jetMom);
  // Rotate phi and keep eta unchanged
  Double_t etaTilted = jet3mom.Eta();
  Double_t phiTilted = TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;
  if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t deta = track3mom.Eta() - etaTilted;
    Double_t dphi = TMath::Abs(track3mom.Phi() - phiTilted);
    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
    Double_t dR = TMath::Sqrt(deta * deta + dphi * dphi);

    if(dR<=radius){
      outputlist->Add(track);
      sumPt += track->Pt();
    }
  }

}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetTracksTiltedwrpJetAxisWindow(Float_t alpha, TList* inputlist, TList* outputlist, AliAODJet* jet, Double_t radius,Double_t& sumPt,Double_t &normFactor)
{
  // List of tracks in cone perpendicular to the jet azimuthal direction

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);

  TVector3 jet3mom(jetMom);
  // Rotate phi and keep eta unchanged
  Double_t etaTilted = jet3mom.Eta();
  Double_t phiTilted = TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;
  if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++)
    {

      AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
      if(!track)continue;
      Float_t trackEta = track->Eta();
      Float_t trackPhi = track->Phi();

      if( ( phiTilted-radius >= 0 ) && ( phiTilted+radius <= 2*TMath::Pi()))
	{
	  if((trackPhi<=phiTilted+radius) && 
	     (trackPhi>=phiTilted-radius) &&
	     (trackEta<=fTrackEtaMax)&&(trackEta>=fTrackEtaMin)) // 0.9 and - 0.9
	    {
	      outputlist->Add(track);
	      sumPt += track->Pt();
	    }
	}
      else if( phiTilted-radius < 0 ) 
	{
	  if((( trackPhi < phiTilted+radius ) ||
	      ( trackPhi > 2*TMath::Pi()-(radius-phiTilted) )) &&
	     (( trackEta <= fTrackEtaMax ) && ( trackEta >= fTrackEtaMin ))) 
	    {
	      outputlist->Add(track);
	      sumPt += track->Pt();
	    }
	}
      else if( phiTilted+radius > 2*TMath::Pi() )
	{
	  if((( trackPhi > phiTilted-radius ) ||
              ( trackPhi < phiTilted+radius-2*TMath::Pi() )) &&
	     (( trackEta <= fTrackEtaMax ) && ( trackEta >= fTrackEtaMin ))) 
	    {
	      outputlist->Add(track);
	      sumPt += track->Pt();
	    }
	}
    }

  // Jet area - Temporarily added should be modified with the proper jet area value
  Float_t areaJet = CalcJetArea(etaTilted,radius);
  Float_t areaTilted = 2*radius*(fTrackEtaMax-fTrackEtaMin);

  normFactor = (Float_t) 1. / (areaJet / areaTilted);

}


// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetTracksOutOfNJets(Int_t nCases, TList* inputlist, TList* outputlist, TList* jetlist, Double_t& sumPt)
{
  // List of tracks outside cone around N jet axis  
  // Particles taken randomly

  sumPt = 0;
  //  Int_t   nj  = jetlist->GetSize();
  Float_t rc  = GetFFRadius();
  Float_t rcl = GetFFBckgRadius();

  // Estimate jet and background areas
  Float_t* areaJet = new Float_t[nCases];
  memset(areaJet, 0, sizeof(Float_t) * nCases);
  Float_t* areaJetLarge = new Float_t[nCases];
  memset(areaJetLarge, 0, sizeof(Float_t) * nCases);
  Float_t areaFull = (fTrackEtaMax-fTrackEtaMin)*(fTrackPhiMax-fTrackPhiMin);
  Float_t areaOut = areaFull;

  //estimate jets and background areas
  Int_t nOut = 0;
  Int_t ijet = 0;
  TList* templist = new TList();
  TClonesArray *vect3Jet = new TClonesArray("TVector3",nCases);

  for(Int_t ij=0; ij<nCases; ++ij) 
    {
      // Get jet information
      AliAODJet* jet = dynamic_cast<AliAODJet*>(jetlist->At(ij));
      if(!jet)continue;
      TVector3 jet3mom;
      jet3mom.SetPtEtaPhi(jet->Pt(),jet->Eta(),jet->Phi());
      new((*vect3Jet)[ijet]) TVector3((TVector3)jet3mom);
      Float_t etaJet = (Float_t)((TVector3*) vect3Jet->At(ij))->Eta();

      // Jet area
      areaJet[ij] = CalcJetArea(etaJet,rc);

      // Area jet larger angle
      areaJetLarge[ij] = CalcJetArea(etaJet,rcl);

      // Outside jet area
      areaOut = areaOut - areaJetLarge[ij];
      ijet++;
    }

  // List of all tracks outside jet areas
  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){
    
    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);
    
    Double_t *dR = new Double_t[nCases];
    for(Int_t ij=0; ij<nCases; ij++)
	dR[ij] = (Double_t)((TVector3*) vect3Jet->At(ij))->DeltaR(track3mom);

    if((nCases==1 && (dR[0]>rcl)) ||
       (nCases==2 && (dR[0]>rcl && dR[1]>rcl)) ||
       (nCases==3 && (dR[0]>rcl && dR[1]>rcl && dR[2]>rcl)))
      {
	templist->Add(track);
	nOut++;
      }
    delete [] dR;
  }

  // Take tracks randomly
  Int_t nScaled = (Int_t) (nOut * areaJet[0] / areaOut + 0.5);
  TArrayI* ar = new TArrayI(nOut);

  for(Int_t init=0; init<nOut; init++)
    (*ar)[init] = init;

  Int_t *randIndex = new Int_t[nScaled];
  for(Int_t init2=0; init2<nScaled; init2++)
    randIndex[init2] = -1;

  // Select nScaled different random numbers in nOut
  for(Int_t i=0; i<nScaled; i++)
    {
      Int_t* tmpArr = new Int_t[nOut-i];
      Int_t temp = fRandom->Integer(nOut-i);
      for(Int_t ind = 0; ind< ar->GetSize()-1; ind++)
        {
          if(ind<temp) tmpArr[ind] = (*ar)[ind];
          else tmpArr[ind] = (*ar)[ind+1];
        }
      randIndex[i] = (*ar)[temp];

      ar->Set(nOut-i-1,tmpArr);

      delete [] tmpArr;

    }

  for(Int_t ipart=0; ipart<nScaled; ipart++)
    {
      AliVParticle* track = (AliVParticle*)(templist->At(randIndex[ipart]));
      outputlist->Add(track);
      sumPt += track->Pt();
    }

  outputlist->Sort();

  delete vect3Jet;
  delete templist;
  delete [] areaJetLarge;
  delete [] areaJet;
  delete ar;
  delete [] randIndex;

}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetTracksOutOfNJetsStat(Int_t nCases, TList* inputlist, TList* outputlist, TList* jetlist, Double_t& sumPt, Double_t &normFactor)
{
  // List of tracks outside cone around N jet axis  
  // All particles taken + final scaling factor 

  sumPt = 0;
  Float_t rc  = GetFFRadius();
  Float_t rcl = GetFFBckgRadius();

  // Estimate jet and background areas
  Float_t* areaJet = new Float_t[nCases];
  memset(areaJet, 0, sizeof(Float_t) * nCases);
  Float_t* areaJetLarge = new Float_t[nCases];
  memset(areaJetLarge, 0, sizeof(Float_t) * nCases);
  Float_t areaFull = (fTrackEtaMax-fTrackEtaMin)*(fTrackPhiMax-fTrackPhiMin);
  Float_t areaOut = areaFull;

  //estimate jets and background areas
  Int_t nOut = 0;
  Int_t ijet = 0;
  TClonesArray *vect3Jet = new TClonesArray("TVector3",nCases);

  for(Int_t ij=0; ij<nCases; ++ij) 
    {
      // Get jet information
      AliAODJet* jet = dynamic_cast<AliAODJet*>(jetlist->At(ij));
      if(!jet)continue;
      TVector3 jet3mom;
      jet3mom.SetPtEtaPhi(jet->Pt(),jet->Eta(),jet->Phi());
      new((*vect3Jet)[ijet]) TVector3((TVector3)jet3mom);
      Float_t etaJet = (Float_t)((TVector3*) vect3Jet->At(ij))->Eta();

      // Jet area
      areaJet[ij] = CalcJetArea(etaJet,rc);

      // Area jet larger angle
      areaJetLarge[ij] = CalcJetArea(etaJet,rcl);

      // Outside jets area
      areaOut = areaOut - areaJetLarge[ij];
      ijet++;
    }

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){
    
    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);
    
    Double_t *dR = new Double_t[nCases];
    for(Int_t ij=0; ij<nCases; ij++)
	dR[ij] = (Double_t)((TVector3*) vect3Jet->At(ij))->DeltaR(track3mom);

    if((nCases==0) ||
       (nCases==1 && (dR[0]>rcl)) ||
       (nCases==2 && (dR[0]>rcl && dR[1]>rcl)) ||
       (nCases==3 && (dR[0]>rcl && dR[1]>rcl && dR[2]>rcl)))
      {
	outputlist->Add(track);
	sumPt += track->Pt();
	nOut++;
      }
    delete [] dR;
  }

  if(nCases==0) areaJet[0] = TMath::Pi()*rc*rc;
  normFactor = (Float_t) 1./(areaJet[0] / areaOut); 

  outputlist->Sort();

  delete vect3Jet;
  delete [] areaJetLarge;
  delete [] areaJet;

}

// ______________________________________________________________________________________________________________________________________________________
Float_t AliAnalysisTaskFragmentationFunction::CalcJetArea(const Float_t etaJet, const Float_t rc) const
{
  // calculate area of jet with eta etaJet and radius rc

  Float_t detamax = etaJet + rc;
  Float_t detamin = etaJet - rc;
  Float_t accmax = 0.0; Float_t accmin = 0.0;
  if(detamax > fTrackEtaMax){ // sector outside etamax
    Float_t h = fTrackEtaMax - etaJet;
    accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
  }
  if(detamin < fTrackEtaMin){ // sector outside etamin
    Float_t h = fTrackEtaMax + etaJet;
    accmin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
  }
  Float_t areaJet = rc*rc*TMath::Pi() - accmax - accmin;
  
  return areaJet;

}

// ___________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetClusterTracksOutOf1Jet(AliAODJet* jet, TList* outputlist, Double_t &normFactor)
{
  // fill tracks from bckgCluster branch in list, 
  // for all clusters outside jet cone 
  // sum up total area of clusters

  Double_t rc  = GetFFRadius();
  Double_t rcl = GetFFBckgRadius();
    
  Double_t areaTotal   = 0;
  Double_t sumPtTotal  = 0;

  for(Int_t ij=0; ij<fBckgJetsRec->GetEntries(); ++ij){
      
    AliAODJet* bgrCluster = (AliAODJet*)(fBckgJetsRec->At(ij)); // not 'recCuts': use all clusters in full eta range
    
    Double_t dR = jet->DeltaR(bgrCluster);  
	 
    if(dR<rcl) continue;
	 
    Double_t clusterPt = bgrCluster->Pt();
    Double_t area      = bgrCluster->EffectiveAreaCharged();
    areaTotal  += area;
    sumPtTotal += clusterPt;
    
    Int_t nTracksJet = bgrCluster->GetRefTracks()->GetEntries();

    for(Int_t it = 0; it<nTracksJet; it++){
	
      AliVParticle*   track = dynamic_cast<AliVParticle*>(bgrCluster->GetTrack(it));
	
      Float_t trackPt  = track->Pt();
      Float_t trackEta = track->Eta();
      Float_t trackPhi = TVector2::Phi_0_2pi(track->Phi());
	
      if(trackEta < fTrackEtaMin || trackEta > fTrackEtaMax) continue;
      if(trackPhi < fTrackPhiMin || trackPhi > fTrackPhiMax) continue;
      if(trackPt  < fTrackPtCut) continue;
	
      outputlist->Add(track);
    }
  }
    
  Double_t areaJet = TMath::Pi()*rc*rc;
  if(areaTotal) normFactor = (Float_t) 1./(areaJet / areaTotal); 

  outputlist->Sort();
}    

// _______________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetClusterTracksMedian(TList* outputlist, Double_t &normFactor)
{
  // fill tracks from bckgCluster branch, 
  // using cluster with median density (odd number of clusters) 
  // or picking randomly one of the two closest to median (even number)
  
  Int_t nBckgClusters = fBckgJetsRec->GetEntries(); // not 'recCuts': use all clusters in full eta range

  Double_t* bgrDensity = new Double_t[nBckgClusters];
  Int_t*    indices    = new Int_t[nBckgClusters];
    
  for(Int_t ij=0; ij<nBckgClusters; ++ij){
      
    AliAODJet* bgrCluster = (AliAODJet*)(fBckgJetsRec->At(ij));
    Double_t clusterPt    = bgrCluster->Pt();
    Double_t area         = bgrCluster->EffectiveAreaCharged();
      
    Double_t density = 0;
    if(area>0) density = clusterPt/area;

    bgrDensity[ij] = density;
    indices[ij]    = ij;
  }
   
  TMath::Sort(nBckgClusters, bgrDensity, indices); 
  
  // get median cluster

  AliAODJet* medianCluster = 0;
  Double_t   medianDensity = 0;

  if(TMath::Odd(nBckgClusters)){

    Int_t medianIndex = indices[(Int_t) (0.5*(nBckgClusters-1))];
    medianCluster = (AliAODJet*)(fBckgJetsRec->At(medianIndex));
    
    Double_t clusterPt = medianCluster->Pt();
    Double_t area      = medianCluster->EffectiveAreaCharged();
    
    if(area>0) medianDensity = clusterPt/area;
  }
  else{

    Int_t medianIndex1 = indices[(Int_t) (0.5*nBckgClusters-1)];
    Int_t medianIndex2 = indices[(Int_t) (0.5*nBckgClusters)];

    AliAODJet* medianCluster1 = (AliAODJet*)(fBckgJetsRec->At(medianIndex1));
    AliAODJet* medianCluster2 = (AliAODJet*)(fBckgJetsRec->At(medianIndex2));
    
    Double_t density1 = 0;
    Double_t clusterPt1 = medianCluster1->Pt();
    Double_t area1      = medianCluster1->EffectiveAreaCharged();
    if(area1>0) density1 = clusterPt1/area1;

    Double_t density2 = 0;
    Double_t clusterPt2 = medianCluster2->Pt();
    Double_t area2      = medianCluster2->EffectiveAreaCharged();
    if(area2>0) density2 = clusterPt2/area2;
    
    medianDensity = 0.5*(density1+density2);
    
    medianCluster = ( (fRandom->Rndm()>0.5) ? medianCluster1 : medianCluster2 );  // select one randomly to avoid adding areas
  }
  
  Int_t nTracksJet = medianCluster->GetRefTracks()->GetEntries();

  for(Int_t it = 0; it<nTracksJet; it++){
	
    AliVParticle* track = dynamic_cast<AliVParticle*>(medianCluster->GetTrack(it));
	
    Float_t trackPt  = track->Pt();
    Float_t trackEta = track->Eta();
    Float_t trackPhi = TVector2::Phi_0_2pi(track->Phi());
    
    if(trackEta < fTrackEtaMin || trackEta > fTrackEtaMax) continue;
    if(trackPhi < fTrackPhiMin || trackPhi > fTrackPhiMax) continue;
    if(trackPt  < fTrackPtCut) continue;
	
    outputlist->Add(track);
  }	
    
  Double_t areaMedian = medianCluster->EffectiveAreaCharged();
  Double_t areaJet = TMath::Pi()*GetFFRadius()*GetFFRadius();
  
  if(areaMedian) normFactor = (Float_t) 1./(areaJet / areaMedian); 
  
  outputlist->Sort();

  delete[] bgrDensity;
  delete[] indices; 
}    

// ______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillBckgHistos(Int_t type, TList* inputtracklist, TList* inputjetlist, AliAODJet* jet, Float_t leadTrackPt, TLorentzVector* leadTrackV,
							  AliFragFuncHistos* ffbckghistocuts, AliFragFuncHistos* ffbckghistoleading, 
							  AliFragFuncIntraJetHistos* ijbckghistocuts, AliFragFuncIntraJetHistos* ijbckghistoleading,   
							  AliFragFuncQATrackHistos* qabckghistocuts){

  // List of tracks outside jets for background study
  TList* tracklistout2jets     = new TList();
  TList* tracklistout3jets     = new TList();
  TList* tracklistout2jetsStat = new TList();
  TList* tracklistout3jetsStat = new TList();
  Double_t sumPtOut2Jets       = 0.;
  Double_t sumPtOut3Jets       = 0.;
  Double_t sumPtOut2JetsStat   = 0.;
  Double_t sumPtOut3JetsStat   = 0.;
  Double_t normFactor2Jets     = 0.;
  Double_t normFactor3Jets     = 0.;

  Int_t nRecJetsCuts = fJetsRecCuts->GetEntries();

  if(nRecJetsCuts>1) {
    GetTracksOutOfNJets(2,inputtracklist, tracklistout2jets, inputjetlist, sumPtOut2Jets);
    GetTracksOutOfNJetsStat(2,inputtracklist, tracklistout2jetsStat, inputjetlist,sumPtOut2JetsStat, normFactor2Jets);

  }
  if(nRecJetsCuts>2) {
    GetTracksOutOfNJets(3,inputtracklist, tracklistout3jets, inputjetlist, sumPtOut3Jets);
    GetTracksOutOfNJetsStat(3,inputtracklist, tracklistout3jetsStat, inputjetlist, sumPtOut3JetsStat, normFactor3Jets);
  }

  if(type==kBckgOutLJ || type==kBckgOutAJ)
    {
      TList* tracklistoutleading       = new TList();
      Double_t sumPtOutLeading     = 0.; 
      GetTracksOutOfNJets(1,inputtracklist, tracklistoutleading, inputjetlist, sumPtOutLeading);
      if(type==kBckgOutLJ) fh1OutLeadingMult->Fill(tracklistoutleading->GetSize());
      
      for(Int_t it=0; it<tracklistoutleading->GetSize(); ++it){
	
	AliVParticle* trackVP   = (AliVParticle*)(tracklistoutleading->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOutLJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	    
	    // Fill track QA for background
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==1 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistoutleading->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOutLJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
        // All cases included
        if(nRecJetsCuts==1 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
      }
      delete tracklistoutleading;
    }
  if(type==kBckgOutLJStat || type==kBckgOutAJStat)
    { 
      TList* tracklistoutleadingStat   = new TList();
      Double_t sumPtOutLeadingStat = 0.; 
      Double_t normFactorLeading   = 0.;

      GetTracksOutOfNJetsStat(1,inputtracklist, tracklistoutleadingStat, inputjetlist, sumPtOutLeadingStat, normFactorLeading);
      if(type==kBckgOutLJStat) fh1OutLeadingStatMult->Fill(tracklistoutleadingStat->GetSize());

      for(Int_t it=0; it<tracklistoutleadingStat->GetSize(); ++it){
	
	AliVParticle* trackVP   = dynamic_cast<AliVParticle*>(tracklistoutleadingStat->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	// Stat plots
	if(type==kBckgOutLJStat)
	  {
	    if(fFFMode)ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorLeading);
	    
	    if(fFFMode)ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorLeading);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactorLeading);
	    
	    // Fill track QA for background
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt); // OB added bgr QA
	  }

	// All cases included
	if(nRecJetsCuts==1 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorLeading );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorLeading);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactorLeading );

	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt ); // OB added bgr QA 

	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistoutleadingStat->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOutLJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorLeading);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactorLeading);
          }
        // All cases included
        if(nRecJetsCuts==1 && type==kBckgOutLJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorLeading);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactorLeading );
          }
      }

      delete tracklistoutleadingStat;
    }

  if(type==kBckgPerp)
    {
      Double_t sumPtPerp = 0.;
      TList* tracklistperp = new TList();
      GetTracksTiltedwrpJetAxis(TMath::Pi()/2., inputtracklist,tracklistperp,jet,GetFFRadius(),sumPtPerp);
      fh1PerpMult->Fill(tracklistperp->GetSize());
      
      for(Int_t it=0; it<tracklistperp->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistperp->At(it));
	if(!trackVP)continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	
	if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	
	// Fill track QA for background
	if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistperp->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
         if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
      }

      delete tracklistperp;
    }

 if(type==kBckgASide)
    {
      Double_t sumPtASide = 0.;
      TList* tracklistaside = new TList();
      GetTracksTiltedwrpJetAxis(TMath::Pi(),inputtracklist,tracklistaside,jet,GetFFRadius(),sumPtASide);
      fh1ASideMult->Fill(tracklistaside->GetSize());

      for(Int_t it=0; it<tracklistaside->GetSize(); ++it){

        AliVParticle*   trackVP = (AliVParticle*)(tracklistaside->At(it));
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();

        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
        if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );

        if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
        if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);

        delete trackV;
      }
      if(tracklistaside->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
         if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
      }

      delete tracklistaside;
    }

  if(type==kBckgASideWindow)
    {
      Double_t normFactorASide = 0.;
      Double_t sumPtASideW = 0.;
      TList* tracklistasidew = new TList();
      GetTracksTiltedwrpJetAxisWindow(TMath::Pi(),inputtracklist,tracklistasidew,jet,GetFFRadius(),sumPtASideW,normFactorASide);
      fh1ASideWindowMult->Fill(tracklistasidew->GetSize());

      for(Int_t it=0; it<tracklistasidew->GetSize(); ++it){

        AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistasidew->At(it));
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();
        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorASide);
        if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorASide);

        if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorASide );
        if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactorASide );

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt, kFALSE, normFactorASide);

        delete trackV;
      }
      if(tracklistasidew->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorASide);
         if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactorASide);
      }

      delete tracklistasidew;
    }

  if(type==kBckgPerpWindow)
    {
      Double_t normFactorPerp = 0.;
      Double_t sumPtPerpW = 0.;
      TList* tracklistperpw = new TList();
      GetTracksTiltedwrpJetAxisWindow(TMath::Pi()/2.,inputtracklist,tracklistperpw,jet,GetFFRadius(),sumPtPerpW,normFactorPerp);
      fh1PerpWindowMult->Fill(tracklistperpw->GetSize());

      for(Int_t it=0; it<tracklistperpw->GetSize(); ++it){

        AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistperpw->At(it));
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();
        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorPerp);
        if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorPerp);

        if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorPerp );
        if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactorPerp );

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt, kFALSE, normFactorPerp);

        delete trackV;
      }
      if(tracklistperpw->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorPerp);
         if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactorPerp);
      }

      delete tracklistperpw;
    }


  if(type==kBckgOut2J || type==kBckgOutAJ)
    {
      if(type==kBckgOut2J) fh1Out2JetsMult->Fill(tracklistout2jets->GetSize());
      for(Int_t it=0; it<tracklistout2jets->GetSize(); ++it){
    
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jets->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOut2J)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	    
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==2 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistout2jets->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOut2J)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
        // All cases included
        if(nRecJetsCuts==2 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
      }

    }

  if(type==kBckgOut2JStat || type==kBckgOutAJStat)
    {
      for(Int_t it=0; it<tracklistout2jetsStat->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jetsStat->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(type==kBckgOut2JStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor2Jets );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor2Jets);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactor2Jets );
	    
	    //if(fQAMode&1) 	qabckghistocuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==2 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor2Jets );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor2Jets);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactor2Jets );

	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt ); // OB added bgr QA 

	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistout2jetsStat->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOut2JStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor2Jets);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactor2Jets);
          }
        // All cases included
        if(nRecJetsCuts==2 && type==kBckgOutAJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor2Jets);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactor2Jets );
          }
      }

    }

  if(type==kBckgOut3J || type==kBckgOutAJ)
    {
      if(type==kBckgOut3J) fh1Out3JetsMult->Fill(tracklistout3jets->GetSize());
      
      for(Int_t it=0; it<tracklistout3jets->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jets->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(type==kBckgOut3J)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	    
	    qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==3 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV );
	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistout3jets->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOut3J)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
        // All cases included
        if(nRecJetsCuts==3 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt );
          }
      }


    }

  if(type==kBckgOut3JStat || type==kBckgOutAJStat)
    {
      for(Int_t it=0; it<tracklistout3jetsStat->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jetsStat->At(it));
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOut3JStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor3Jets);
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor3Jets);
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor3Jets);
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactor3Jets);
	    
	    //if(fQAMode&1) 	qabckghistocuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==3 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor3Jets );
	    if(fIJMode) ijbckghistocuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor3Jets);
	    
	    if(fFFMode) ffbckghistoleading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor3Jets );
	    if(fIJMode) ijbckghistoleading->FillIntraJet( trackV, leadTrackV, normFactor3Jets );

	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt ); // OB added bgr QA

	  }
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistout3jetsStat->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(type==kBckgOut3JStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor3Jets);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactor3Jets);
          }
        // All cases included
        if(nRecJetsCuts==3 && type==kBckgOutAJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor3Jets);
            if(fFFMode) ffbckghistoleading->FillFF( -1., leadTrackPt , incrementJetPt, normFactor3Jets );
          }
      }

    }

  if(type==kBckgClustersOutLeading){ // clusters bgr: all tracks in clusters out of leading jet
    
    TList* tracklistClustersOutLeading = new TList();
    Double_t normFactorClusters = 0;
    Float_t jetPt   = jet->Pt();
    
    GetClusterTracksOutOf1Jet(jet, tracklistClustersOutLeading, normFactorClusters);
    fh1OutClustersMult->Fill(tracklistClustersOutLeading->GetSize());
    
    for(Int_t it=0; it<tracklistClustersOutLeading->GetSize(); ++it){
      
      AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistClustersOutLeading->At(it));
      TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
      
      Float_t trackPt = trackVP->Pt();
      
      Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
      
      if(fFFMode)   ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorClusters );
      if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt ); 

      delete trackV;
    }
    
    delete tracklistClustersOutLeading;
    
  }
  
  if(type == kBckgClusters){ // clusters bgr: all tracks in 'median cluster' 
    
    TList* tracklistClustersMedian = new TList();
    Double_t normFactorClusters = 0;
    Float_t jetPt = jet->Pt();
    
    GetClusterTracksMedian(tracklistClustersMedian, normFactorClusters); 
    fh1MedianClustersMult->Fill(tracklistClustersMedian->GetSize());
    
    for(Int_t it=0; it<tracklistClustersMedian->GetSize(); ++it){
      
      AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistClustersMedian->At(it));
      TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
      
      Float_t trackPt = trackVP->Pt();
      
      Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
      
      if(fFFMode)   ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorClusters );
      if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt );
      
      delete trackV;
    }
    
    delete tracklistClustersMedian;
  }
  
  delete tracklistout2jets;
  delete tracklistout3jets;
  delete tracklistout2jetsStat;
  delete tracklistout3jetsStat;
  
}

// ______________________________________________________________________________________________________________________________________________________
AliAODJet* AliAnalysisTaskFragmentationFunction::GetAODBckgSubJet(AliAODJet* jet, Int_t method)
{
  // correct jet pt for (mean bgr energy density) x (jet area) 

  if(!fBranchRecBackJets.Length()){
    if(fDebug>0)Printf("%s:%d background branch name not set",(char*)__FILE__,__LINE__);
    return 0;
  }

  static AliAODJetEventBackground*   externalBackground =  (AliAODJetEventBackground*)(fAOD->FindListObject(fBranchRecBackJets.Data()));

  if(!externalBackground){
    if(fDebug>0)Printf("%s:%d no external background object found in %s",(char*)__FILE__,__LINE__,fBranchRecBackJets.Data());
    return 0;
  }

  Float_t rho = externalBackground->GetBackground(method);

  // Calculate background and subtract it from jet pt
  Float_t ptBack = rho*jet->EffectiveAreaCharged();
  Float_t ptSub = jet->Pt()-ptBack;

  // Get px, py, pz from eta, phi, pt
  TLorentzVector vecSub;
  AliAODJet *tmpJet = 0;
  if(ptSub>=0){
    vecSub.SetPtEtaPhiM(ptSub,jet->Eta(),jet->Phi(),0.);
    tmpJet = new AliAODJet(vecSub.Px(),vecSub.Py(),vecSub.Pz(),vecSub.P());
  }

  return tmpJet;

}


