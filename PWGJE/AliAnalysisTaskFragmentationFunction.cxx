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
#include "TAxis.h"

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
using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliAnalysisTaskFragmentationFunction)

//____________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction()
   : AliAnalysisTaskSE()
   ,fESD(0)
   ,fAOD(0)
   ,fAODJets(0)  
   ,fAODExtension(0)
   ,fNonStdFile("")
   ,fBranchRecJets("jets")
   ,fBranchRecBckgClusters("")
   ,fBranchGenJets("")
   ,fBranchEmbeddedJets("")
   ,fTrackTypeGen(0)
   ,fJetTypeGen(0)
   ,fJetTypeRecEff(0)
   ,fUseAODInputJets(kTRUE)
   ,fFilterMask(0)
   ,fUsePhysicsSelection(kTRUE)
   ,fEvtSelectionMask(0)
   ,fEventClass(0)
   ,fMaxVertexZ(10)
   ,fTrackPtCut(0)
   ,fTrackEtaMin(0)
   ,fTrackEtaMax(0)
   ,fTrackPhiMin(0)
   ,fTrackPhiMax(0)
   ,fUseExtraTracks(0)
   ,fUseExtraTracksBgr(0)
   ,fCutFractionPtEmbedded(0)
   ,fUseEmbeddedJetAxis(0)
   ,fUseEmbeddedJetPt(0)
   ,fJetPtCut(0)
   ,fJetEtaMin(0)
   ,fJetEtaMax(0)
   ,fJetPhiMin(0)
   ,fJetPhiMax(0)
   ,fFFRadius(0)
   ,fFFMinLTrackPt(-1)
   ,fFFMaxTrackPt(-1)
   ,fFFMinnTracks(0)
   ,fFFBckgRadius(0)
   ,fBckgMode(0)
   ,fQAMode(0)
   ,fFFMode(0)
   ,fEffMode(0)
   ,fJSMode(0)
   ,fAvgTrials(0)
   ,fTracksRecCuts(0)
   ,fTracksGen(0)
   ,fTracksAODMCCharged(0)
   ,fTracksAODMCChargedSecNS(0)
   ,fTracksAODMCChargedSecS(0)
   ,fTracksRecQualityCuts(0)
   ,fJetsRec(0)
   ,fJetsRecCuts(0)
   ,fJetsGen(0)
   ,fJetsRecEff(0)
   ,fJetsEmbedded(0)
   ,fBckgJetsRec(0)
   ,fBckgJetsRecCuts(0)
   ,fBckgJetsGen(0)
   ,fQATrackHistosRecCuts(0)
   ,fQATrackHistosGen(0)
   ,fQAJetHistosRec(0)
   ,fQAJetHistosRecCuts(0)
   ,fQAJetHistosRecCutsLeading(0)
   ,fQAJetHistosGen(0)
   ,fQAJetHistosGenLeading(0)
   ,fQAJetHistosRecEffLeading(0)
   ,fFFHistosRecCuts(0)
   ,fFFHistosGen(0)
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
   ,fh1nEmbeddedJets(0)
   ,fh1nRecBckgJetsCuts(0)
   ,fh1nGenBckgJets(0)
   ,fh2PtRecVsGenPrim(0)
   ,fh2PtRecVsGenSec(0)
   ,fQATrackHistosRecEffGen(0)  
   ,fQATrackHistosRecEffRec(0)
   ,fQATrackHistosSecRecNS(0)   
   ,fQATrackHistosSecRecS(0)   
   ,fQATrackHistosSecRecSsc(0)   
   ,fFFHistosRecEffRec(0)
   ,fFFHistosSecRecNS(0)
   ,fFFHistosSecRecS(0)
   ,fFFHistosSecRecSsc(0)
   // Background 
   ,fh1BckgMult0(0)
   ,fh1BckgMult1(0)
   ,fh1BckgMult2(0)
   ,fh1BckgMult3(0)
   ,fh1BckgMult4(0)
   ,fh1FractionPtEmbedded(0)
   ,fh1IndexEmbedded(0)
   ,fh2DeltaPtVsJetPtEmbedded(0)
   ,fh2DeltaPtVsRecJetPtEmbedded(0)
   ,fh1DeltaREmbedded(0)
   ,fQABckgHisto0RecCuts(0)  
   ,fQABckgHisto0Gen(0)      
   ,fQABckgHisto1RecCuts(0)  
   ,fQABckgHisto1Gen(0)      
   ,fQABckgHisto2RecCuts(0)  
   ,fQABckgHisto2Gen(0)
   ,fQABckgHisto3RecCuts(0)
   ,fQABckgHisto3Gen(0)
   ,fQABckgHisto4RecCuts(0)
   ,fQABckgHisto4Gen(0)
   ,fFFBckgHisto0RecCuts(0)
   ,fFFBckgHisto0Gen(0)       
   ,fFFBckgHisto1RecCuts(0)
   ,fFFBckgHisto1Gen(0)       
   ,fFFBckgHisto2RecCuts(0)
   ,fFFBckgHisto2Gen(0)       
   ,fFFBckgHisto3RecCuts(0)
   ,fFFBckgHisto3Gen(0)       
   ,fFFBckgHisto4RecCuts(0)
   ,fFFBckgHisto4Gen(0)       
   ,fFFBckgHisto0RecEffRec(0)
   ,fFFBckgHisto0SecRecNS(0)  
   ,fFFBckgHisto0SecRecS(0)   
   ,fFFBckgHisto0SecRecSsc(0)
    // jet shape   
   ,fProNtracksLeadingJet(0)
   ,fProDelR80pcPt(0)
   ,fProNtracksLeadingJetGen(0)
   ,fProDelR80pcPtGen(0)
   ,fProNtracksLeadingJetBgrPerp2(0)
   ,fProNtracksLeadingJetRecPrim(0)  
   ,fProDelR80pcPtRecPrim(0)
   ,fProNtracksLeadingJetRecSecNS(0) 
   ,fProNtracksLeadingJetRecSecS(0)  
   ,fProNtracksLeadingJetRecSecSsc(0)

   ,fRandom(0)
{
   // default constructor
  fBckgType[0] = 0;
  fBckgType[1] = 0;
  fBckgType[2] = 0;
  fBckgType[3] = 0;
  fBckgType[4] = 0;

  for(Int_t ii=0; ii<5; ii++){
    fProDelRPtSum[ii]          = 0;
    fProDelRPtSumGen[ii]       = 0;
    fProDelRPtSumBgrPerp2[ii]  = 0;
    fProDelRPtSumRecPrim[ii]   = 0;
    fProDelRPtSumRecSecNS[ii]  = 0;
    fProDelRPtSumRecSecS[ii]   = 0;
    fProDelRPtSumRecSecSsc[ii] = 0;
  }
}

//_______________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fAODJets(0)  
  ,fAODExtension(0)
  ,fNonStdFile("")
  ,fBranchRecJets("jets")
  ,fBranchRecBckgClusters("")
  ,fBranchGenJets("")
  ,fBranchEmbeddedJets("")
  ,fTrackTypeGen(0)
  ,fJetTypeGen(0)
  ,fJetTypeRecEff(0)
  ,fUseAODInputJets(kTRUE)
  ,fFilterMask(0)
  ,fUsePhysicsSelection(kTRUE)
  ,fEvtSelectionMask(0)
  ,fEventClass(0)
  ,fMaxVertexZ(10)
  ,fTrackPtCut(0)
  ,fTrackEtaMin(0)
  ,fTrackEtaMax(0)
  ,fTrackPhiMin(0)
  ,fTrackPhiMax(0)
  ,fUseExtraTracks(0)
  ,fUseExtraTracksBgr(0)
  ,fCutFractionPtEmbedded(0)
  ,fUseEmbeddedJetAxis(0)
  ,fUseEmbeddedJetPt(0)  
  ,fJetPtCut(0)
  ,fJetEtaMin(0)
  ,fJetEtaMax(0)
  ,fJetPhiMin(0)
  ,fJetPhiMax(0)
  ,fFFRadius(0)
  ,fFFMinLTrackPt(-1)
  ,fFFMaxTrackPt(-1)
  ,fFFMinnTracks(0)
  ,fFFBckgRadius(0)
  ,fBckgMode(0)
  ,fQAMode(0)
  ,fFFMode(0)
  ,fEffMode(0)
  ,fJSMode(0)
  ,fAvgTrials(0)
  ,fTracksRecCuts(0)
  ,fTracksGen(0)
  ,fTracksAODMCCharged(0)
  ,fTracksAODMCChargedSecNS(0)
  ,fTracksAODMCChargedSecS(0)
  ,fTracksRecQualityCuts(0)
  ,fJetsRec(0)
  ,fJetsRecCuts(0)
  ,fJetsGen(0)
  ,fJetsRecEff(0)
  ,fJetsEmbedded(0)
  ,fBckgJetsRec(0)
  ,fBckgJetsRecCuts(0)
  ,fBckgJetsGen(0)
  ,fQATrackHistosRecCuts(0)
  ,fQATrackHistosGen(0)
  ,fQAJetHistosRec(0)
  ,fQAJetHistosRecCuts(0)
  ,fQAJetHistosRecCutsLeading(0)
  ,fQAJetHistosGen(0)
  ,fQAJetHistosGenLeading(0)
  ,fQAJetHistosRecEffLeading(0)
  ,fFFHistosRecCuts(0)
  ,fFFHistosGen(0)
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
  ,fh1nEmbeddedJets(0)
  ,fh1nRecBckgJetsCuts(0)
  ,fh1nGenBckgJets(0)
  ,fh2PtRecVsGenPrim(0)
  ,fh2PtRecVsGenSec(0)
  ,fQATrackHistosRecEffGen(0)  
  ,fQATrackHistosRecEffRec(0)
  ,fQATrackHistosSecRecNS(0) 
  ,fQATrackHistosSecRecS(0) 
  ,fQATrackHistosSecRecSsc(0) 
  ,fFFHistosRecEffRec(0)
  ,fFFHistosSecRecNS(0)
  ,fFFHistosSecRecS(0)
  ,fFFHistosSecRecSsc(0)
  // Background
  ,fh1BckgMult0(0)
  ,fh1BckgMult1(0)
  ,fh1BckgMult2(0)
  ,fh1BckgMult3(0)
  ,fh1BckgMult4(0)
  ,fh1FractionPtEmbedded(0)
  ,fh1IndexEmbedded(0)
  ,fh2DeltaPtVsJetPtEmbedded(0)
  ,fh2DeltaPtVsRecJetPtEmbedded(0)
  ,fh1DeltaREmbedded(0)
  ,fQABckgHisto0RecCuts(0)  
  ,fQABckgHisto0Gen(0)      
  ,fQABckgHisto1RecCuts(0)  
  ,fQABckgHisto1Gen(0)      
  ,fQABckgHisto2RecCuts(0)  
  ,fQABckgHisto2Gen(0) 
  ,fQABckgHisto3RecCuts(0)  
  ,fQABckgHisto3Gen(0)
  ,fQABckgHisto4RecCuts(0)  
  ,fQABckgHisto4Gen(0)
  ,fFFBckgHisto0RecCuts(0)
  ,fFFBckgHisto0Gen(0)       
  ,fFFBckgHisto1RecCuts(0)
  ,fFFBckgHisto1Gen(0)       
  ,fFFBckgHisto2RecCuts(0)
  ,fFFBckgHisto2Gen(0)       
  ,fFFBckgHisto3RecCuts(0)
  ,fFFBckgHisto3Gen(0)       
  ,fFFBckgHisto4RecCuts(0)
  ,fFFBckgHisto4Gen(0)       
  ,fFFBckgHisto0RecEffRec(0)
  ,fFFBckgHisto0SecRecNS(0)  
  ,fFFBckgHisto0SecRecS(0)   
  ,fFFBckgHisto0SecRecSsc(0) 
  // jet shape   
  ,fProNtracksLeadingJet(0)
  ,fProDelR80pcPt(0)
  ,fProNtracksLeadingJetGen(0)
  ,fProDelR80pcPtGen(0)
  ,fProNtracksLeadingJetBgrPerp2(0)
  ,fProNtracksLeadingJetRecPrim(0)
  ,fProDelR80pcPtRecPrim(0)
  ,fProNtracksLeadingJetRecSecNS(0) 
  ,fProNtracksLeadingJetRecSecS(0)  
  ,fProNtracksLeadingJetRecSecSsc(0)
  ,fRandom(0)
{
  // constructor
  fBckgType[0] = 0;
  fBckgType[1] = 0;
  fBckgType[2] = 0;
  fBckgType[3] = 0;
  fBckgType[4] = 0;
 
  for(Int_t ii=0; ii<5; ii++){
    fProDelRPtSum[ii]          = 0;
    fProDelRPtSumGen[ii]       = 0;
    fProDelRPtSumBgrPerp2[ii]  = 0;
    fProDelRPtSumRecPrim[ii]   = 0;
    fProDelRPtSumRecSecNS[ii]  = 0;
    fProDelRPtSumRecSecS[ii]   = 0;
    fProDelRPtSumRecSecSsc[ii] = 0;
  }
  
  DefineOutput(1,TList::Class());
}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const  AliAnalysisTaskFragmentationFunction &copy)
  : AliAnalysisTaskSE()
  ,fESD(copy.fESD)
  ,fAOD(copy.fAOD)
  ,fAODJets(copy.fAODJets)  
  ,fAODExtension(copy.fAODExtension)
  ,fNonStdFile(copy.fNonStdFile)
  ,fBranchRecJets(copy.fBranchRecJets)
  ,fBranchRecBckgClusters(copy.fBranchRecBckgClusters)
  ,fBranchGenJets(copy.fBranchGenJets)
  ,fBranchEmbeddedJets(copy.fBranchEmbeddedJets)
  ,fTrackTypeGen(copy.fTrackTypeGen)
  ,fJetTypeGen(copy.fJetTypeGen)
  ,fJetTypeRecEff(copy.fJetTypeRecEff)
  ,fUseAODInputJets(copy.fUseAODInputJets)
  ,fFilterMask(copy.fFilterMask)
  ,fUsePhysicsSelection(copy.fUsePhysicsSelection)
  ,fEvtSelectionMask(copy.fEvtSelectionMask)
  ,fEventClass(copy.fEventClass)
  ,fMaxVertexZ(copy.fMaxVertexZ)
  ,fTrackPtCut(copy.fTrackPtCut)
  ,fTrackEtaMin(copy.fTrackEtaMin)
  ,fTrackEtaMax(copy.fTrackEtaMax)
  ,fTrackPhiMin(copy.fTrackPhiMin)
  ,fTrackPhiMax(copy.fTrackPhiMax)
  ,fUseExtraTracks(copy.fUseExtraTracks)
  ,fUseExtraTracksBgr(copy.fUseExtraTracksBgr)
  ,fCutFractionPtEmbedded(copy.fCutFractionPtEmbedded)
  ,fUseEmbeddedJetAxis(copy.fUseEmbeddedJetAxis)
  ,fUseEmbeddedJetPt(copy.fUseEmbeddedJetPt)
  ,fJetPtCut(copy.fJetPtCut)
  ,fJetEtaMin(copy.fJetEtaMin)
  ,fJetEtaMax(copy.fJetEtaMax)
  ,fJetPhiMin(copy.fJetPhiMin)
  ,fJetPhiMax(copy.fJetPhiMax)
  ,fFFRadius(copy.fFFRadius)
  ,fFFMinLTrackPt(copy.fFFMinLTrackPt)
  ,fFFMaxTrackPt(copy.fFFMaxTrackPt)
  ,fFFMinnTracks(copy.fFFMinnTracks)
  ,fFFBckgRadius(copy.fFFBckgRadius)
  ,fBckgMode(copy.fBckgMode)
  ,fQAMode(copy.fQAMode)
  ,fFFMode(copy.fFFMode)
  ,fEffMode(copy.fEffMode)
  ,fJSMode(copy.fJSMode)
  ,fAvgTrials(copy.fAvgTrials)
  ,fTracksRecCuts(copy.fTracksRecCuts)
  ,fTracksGen(copy.fTracksGen)
  ,fTracksAODMCCharged(copy.fTracksAODMCCharged)
  ,fTracksAODMCChargedSecNS(copy.fTracksAODMCChargedSecNS)
  ,fTracksAODMCChargedSecS(copy.fTracksAODMCChargedSecS)
  ,fTracksRecQualityCuts(copy.fTracksRecQualityCuts)
  ,fJetsRec(copy.fJetsRec)
  ,fJetsRecCuts(copy.fJetsRecCuts)
  ,fJetsGen(copy.fJetsGen)
  ,fJetsRecEff(copy.fJetsRecEff)
  ,fJetsEmbedded(copy.fJetsEmbedded)
  ,fBckgJetsRec(copy.fBckgJetsRec)
  ,fBckgJetsRecCuts(copy.fBckgJetsRecCuts)
  ,fBckgJetsGen(copy.fBckgJetsGen)
  ,fQATrackHistosRecCuts(copy.fQATrackHistosRecCuts)
  ,fQATrackHistosGen(copy.fQATrackHistosGen)
  ,fQAJetHistosRec(copy.fQAJetHistosRec)
  ,fQAJetHistosRecCuts(copy.fQAJetHistosRecCuts)
  ,fQAJetHistosRecCutsLeading(copy.fQAJetHistosRecCutsLeading)
  ,fQAJetHistosGen(copy.fQAJetHistosGen)
  ,fQAJetHistosGenLeading(copy.fQAJetHistosGenLeading)
  ,fQAJetHistosRecEffLeading(copy.fQAJetHistosRecEffLeading)
  ,fFFHistosRecCuts(copy.fFFHistosRecCuts)
  ,fFFHistosGen(copy.fFFHistosGen)
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
  ,fh1nEmbeddedJets(copy.fh1nEmbeddedJets)
  ,fh1nRecBckgJetsCuts(copy.fh1nRecBckgJetsCuts)
  ,fh1nGenBckgJets(copy.fh1nGenBckgJets)
  ,fh2PtRecVsGenPrim(copy.fh2PtRecVsGenPrim)
  ,fh2PtRecVsGenSec(copy.fh2PtRecVsGenSec)
  ,fQATrackHistosRecEffGen(copy.fQATrackHistosRecEffGen)  
  ,fQATrackHistosRecEffRec(copy.fQATrackHistosRecEffRec)  
  ,fQATrackHistosSecRecNS(copy.fQATrackHistosSecRecNS)  
  ,fQATrackHistosSecRecS(copy.fQATrackHistosSecRecS)  
  ,fQATrackHistosSecRecSsc(copy.fQATrackHistosSecRecSsc)  
  ,fFFHistosRecEffRec(copy.fFFHistosRecEffRec)  
  ,fFFHistosSecRecNS(copy.fFFHistosSecRecNS)   
  ,fFFHistosSecRecS(copy.fFFHistosSecRecS)   
  ,fFFHistosSecRecSsc(copy.fFFHistosSecRecSsc)   
  // Background
  ,fh1BckgMult0(copy.fh1BckgMult0)
  ,fh1BckgMult1(copy.fh1BckgMult1)
  ,fh1BckgMult2(copy.fh1BckgMult2)
  ,fh1BckgMult3(copy.fh1BckgMult3)
  ,fh1BckgMult4(copy.fh1BckgMult4)
  ,fh1FractionPtEmbedded(copy.fh1FractionPtEmbedded)
  ,fh1IndexEmbedded(copy.fh1IndexEmbedded)
  ,fh2DeltaPtVsJetPtEmbedded(copy.fh2DeltaPtVsJetPtEmbedded)
  ,fh2DeltaPtVsRecJetPtEmbedded(copy.fh2DeltaPtVsRecJetPtEmbedded)
  ,fh1DeltaREmbedded(copy.fh1DeltaREmbedded)
  ,fQABckgHisto0RecCuts(copy.fQABckgHisto0RecCuts)  
  ,fQABckgHisto0Gen(copy.fQABckgHisto0Gen)      
  ,fQABckgHisto1RecCuts(copy.fQABckgHisto1RecCuts)  
  ,fQABckgHisto1Gen(copy.fQABckgHisto1Gen)      
  ,fQABckgHisto2RecCuts(copy.fQABckgHisto2RecCuts)  
  ,fQABckgHisto2Gen(copy.fQABckgHisto2Gen)
  ,fQABckgHisto3RecCuts(copy.fQABckgHisto3RecCuts)  
  ,fQABckgHisto3Gen(copy.fQABckgHisto3Gen)     
  ,fQABckgHisto4RecCuts(copy.fQABckgHisto4RecCuts)  
  ,fQABckgHisto4Gen(copy.fQABckgHisto4Gen)     
  ,fFFBckgHisto0RecCuts(copy.fFFBckgHisto0RecCuts)
  ,fFFBckgHisto0Gen(copy.fFFBckgHisto0Gen)       
  ,fFFBckgHisto1RecCuts(copy.fFFBckgHisto1RecCuts)
  ,fFFBckgHisto1Gen(copy.fFFBckgHisto1Gen)       
  ,fFFBckgHisto2RecCuts(copy.fFFBckgHisto2RecCuts)
  ,fFFBckgHisto2Gen(copy.fFFBckgHisto2Gen)       
  ,fFFBckgHisto3RecCuts(copy.fFFBckgHisto3RecCuts)
  ,fFFBckgHisto3Gen(copy.fFFBckgHisto3Gen)       
  ,fFFBckgHisto4RecCuts(copy.fFFBckgHisto4RecCuts)
  ,fFFBckgHisto4Gen(copy.fFFBckgHisto4Gen)       
  ,fFFBckgHisto0RecEffRec(copy.fFFBckgHisto0RecEffRec)
  ,fFFBckgHisto0SecRecNS(copy.fFFBckgHisto0SecRecNS)  
  ,fFFBckgHisto0SecRecS(copy.fFFBckgHisto0SecRecS)   
  ,fFFBckgHisto0SecRecSsc(copy.fFFBckgHisto0SecRecSsc) 
  // jet shape   
  ,fProNtracksLeadingJet(copy.fProNtracksLeadingJet)                
  ,fProDelR80pcPt(copy.fProDelR80pcPt)                       
  ,fProNtracksLeadingJetGen(copy.fProNtracksLeadingJetGen)             
  ,fProDelR80pcPtGen(copy.fProDelR80pcPtGen)                    
  ,fProNtracksLeadingJetBgrPerp2(copy.fProNtracksLeadingJetBgrPerp2)        
  ,fProNtracksLeadingJetRecPrim(copy.fProNtracksLeadingJetRecPrim)  
  ,fProDelR80pcPtRecPrim(copy.fProDelR80pcPtRecPrim)
  ,fProNtracksLeadingJetRecSecNS(copy.fProNtracksLeadingJetRecSecNS) 
  ,fProNtracksLeadingJetRecSecS(copy.fProNtracksLeadingJetRecSecS)  
  ,fProNtracksLeadingJetRecSecSsc(copy.fProNtracksLeadingJetRecSecSsc)
  ,fRandom(copy.fRandom)
{
  // copy constructor
  fBckgType[0] = copy.fBckgType[0];
  fBckgType[1] = copy.fBckgType[1];
  fBckgType[2] = copy.fBckgType[2];
  fBckgType[3] = copy.fBckgType[3];
  fBckgType[4] = copy.fBckgType[4];


  for(Int_t ii=0; ii<5; ii++){
    fProDelRPtSum[ii]          = copy.fProDelRPtSum[ii];
    fProDelRPtSumGen[ii]       = copy.fProDelRPtSumGen[ii];
    fProDelRPtSumBgrPerp2[ii]  = copy.fProDelRPtSumBgrPerp2[ii];
    fProDelRPtSumRecPrim[ii]   = copy.fProDelRPtSumRecPrim[ii];
    fProDelRPtSumRecSecNS[ii]  = copy.fProDelRPtSumRecSecNS[ii];
    fProDelRPtSumRecSecS[ii]   = copy.fProDelRPtSumRecSecS[ii];
    fProDelRPtSumRecSecSsc[ii] = copy.fProDelRPtSumRecSecSsc[ii];
  }
}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction& AliAnalysisTaskFragmentationFunction::operator=(const AliAnalysisTaskFragmentationFunction& o)
{
  // assignment
  
  if(this!=&o){

    AliAnalysisTaskSE::operator=(o);
    fESD                           = o.fESD;
    fAOD                           = o.fAOD;
    fAODJets                       = o.fAODJets;  
    fAODExtension                  = o.fAODExtension;
    fNonStdFile                    = o.fNonStdFile;
    fBranchRecJets                 = o.fBranchRecJets;
    fBranchRecBckgClusters         = o.fBranchRecBckgClusters;
    fBranchGenJets                 = o.fBranchGenJets;
    fBranchEmbeddedJets            = o.fBranchEmbeddedJets;
    fTrackTypeGen                  = o.fTrackTypeGen;
    fJetTypeGen                    = o.fJetTypeGen;
    fJetTypeRecEff                 = o.fJetTypeRecEff;
    fUseAODInputJets               = o.fUseAODInputJets;
    fFilterMask                    = o.fFilterMask;
    fUsePhysicsSelection           = o.fUsePhysicsSelection;
    fEvtSelectionMask              = o.fEvtSelectionMask;
    fEventClass                    = o.fEventClass;
    fMaxVertexZ                    = o.fMaxVertexZ;
    fTrackPtCut                    = o.fTrackPtCut;
    fTrackEtaMin                   = o.fTrackEtaMin;
    fTrackEtaMax                   = o.fTrackEtaMax;
    fTrackPhiMin                   = o.fTrackPhiMin;
    fTrackPhiMax                   = o.fTrackPhiMax;
    fUseExtraTracks                = o.fUseExtraTracks;
    fUseExtraTracksBgr             = o.fUseExtraTracksBgr;
    fCutFractionPtEmbedded         = o.fCutFractionPtEmbedded;
    fUseEmbeddedJetAxis            = o.fUseEmbeddedJetAxis;
    fUseEmbeddedJetPt              = o.fUseEmbeddedJetPt;
    fJetPtCut                      = o.fJetPtCut;
    fJetEtaMin                     = o.fJetEtaMin;
    fJetEtaMax                     = o.fJetEtaMax;
    fJetPhiMin                     = o.fJetPhiMin;
    fJetPhiMax                     = o.fJetPhiMin;
    fFFRadius                      = o.fFFRadius;
    fFFMinLTrackPt                 = o.fFFMinLTrackPt;
    fFFMaxTrackPt                  = o.fFFMaxTrackPt;
    fFFMinnTracks                  = o.fFFMinnTracks;
    fFFBckgRadius                  = o.fFFBckgRadius;
    fBckgMode                      = o.fBckgMode;
    fQAMode                        = o.fQAMode;
    fFFMode                        = o.fFFMode;
    fEffMode                       = o.fEffMode;
    fJSMode                        = o.fJSMode;
    fBckgType[0]                   = o.fBckgType[0];
    fBckgType[1]                   = o.fBckgType[1];
    fBckgType[2]                   = o.fBckgType[2];
    fBckgType[3]                   = o.fBckgType[3];
    fBckgType[4]                   = o.fBckgType[4];
    fAvgTrials                     = o.fAvgTrials;
    fTracksRecCuts                 = o.fTracksRecCuts;
    fTracksGen                     = o.fTracksGen;
    fTracksAODMCCharged            = o.fTracksAODMCCharged;
    fTracksAODMCChargedSecNS       = o.fTracksAODMCChargedSecNS;
    fTracksAODMCChargedSecS        = o.fTracksAODMCChargedSecS;
    fTracksRecQualityCuts          = o.fTracksRecQualityCuts;
    fJetsRec                       = o.fJetsRec;
    fJetsRecCuts                   = o.fJetsRecCuts;
    fJetsGen                       = o.fJetsGen;
    fJetsRecEff                    = o.fJetsRecEff;
    fJetsEmbedded                  = o.fJetsEmbedded;
    fBckgJetsRec                   = o.fBckgJetsRec;
    fBckgJetsRecCuts               = o.fBckgJetsRecCuts;
    fBckgJetsGen                   = o.fBckgJetsGen;
    fQATrackHistosRecCuts          = o.fQATrackHistosRecCuts;
    fQATrackHistosGen              = o.fQATrackHistosGen;
    fQAJetHistosRec                = o.fQAJetHistosRec;
    fQAJetHistosRecCuts            = o.fQAJetHistosRecCuts;
    fQAJetHistosRecCutsLeading     = o.fQAJetHistosRecCutsLeading;
    fQAJetHistosGen                = o.fQAJetHistosGen;
    fQAJetHistosGenLeading         = o.fQAJetHistosGenLeading;
    fQAJetHistosRecEffLeading      = o.fQAJetHistosRecEffLeading;
    fFFHistosRecCuts               = o.fFFHistosRecCuts;
    fFFHistosGen                   = o.fFFHistosGen;
    fQATrackHighPtThreshold        = o.fQATrackHighPtThreshold; 
    fFFNBinsJetPt                  = o.fFFNBinsJetPt;    
    fFFJetPtMin                    = o.fFFJetPtMin; 
    fFFJetPtMax                    = o.fFFJetPtMax;
    fFFNBinsPt                     = o.fFFNBinsPt;      
    fFFPtMin                       = o.fFFPtMin;        
    fFFPtMax                       = o.fFFPtMax;        
    fFFNBinsXi                     = o.fFFNBinsXi;      
    fFFXiMin                       = o.fFFXiMin;        
    fFFXiMax                       = o.fFFXiMax;        
    fFFNBinsZ                      = o.fFFNBinsZ;       
    fFFZMin                        = o.fFFZMin;         
    fFFZMax                        = o.fFFZMax;         
    fQAJetNBinsPt                  = o.fQAJetNBinsPt;   
    fQAJetPtMin                    = o.fQAJetPtMin;     
    fQAJetPtMax                    = o.fQAJetPtMax;     
    fQAJetNBinsEta                 = o.fQAJetNBinsEta;  
    fQAJetEtaMin                   = o.fQAJetEtaMin;    
    fQAJetEtaMax                   = o.fQAJetEtaMax;    
    fQAJetNBinsPhi                 = o.fQAJetNBinsPhi;  
    fQAJetPhiMin                   = o.fQAJetPhiMin;    
    fQAJetPhiMax                   = o.fQAJetPhiMax;    
    fQATrackNBinsPt                = o.fQATrackNBinsPt; 
    fQATrackPtMin                  = o.fQATrackPtMin;   
    fQATrackPtMax                  = o.fQATrackPtMax;   
    fQATrackNBinsEta               = o.fQATrackNBinsEta;
    fQATrackEtaMin                 = o.fQATrackEtaMin;  
    fQATrackEtaMax                 = o.fQATrackEtaMax;  
    fQATrackNBinsPhi               = o.fQATrackNBinsPhi;
    fQATrackPhiMin                 = o.fQATrackPhiMin;  
    fQATrackPhiMax                 = o.fQATrackPhiMax;  
    fCommonHistList                = o.fCommonHistList;
    fh1EvtSelection                = o.fh1EvtSelection;
    fh1VertexNContributors         = o.fh1VertexNContributors;
    fh1VertexZ                     = o.fh1VertexZ;
    fh1EvtMult                     = o.fh1EvtMult;
    fh1EvtCent                     = o.fh1EvtCent;
    fh1Xsec                        = o.fh1Xsec;
    fh1Trials                      = o.fh1Trials;
    fh1PtHard                      = o.fh1PtHard;
    fh1PtHardTrials                = o.fh1PtHardTrials;
    fh1nRecJetsCuts                = o.fh1nRecJetsCuts;
    fh1nGenJets                    = o.fh1nGenJets; 
    fh1nRecEffJets                 = o.fh1nRecEffJets;
    fh1nEmbeddedJets               = o.fh1nEmbeddedJets;
    fh2PtRecVsGenPrim              = o.fh2PtRecVsGenPrim;
    fh2PtRecVsGenSec               = o.fh2PtRecVsGenSec; 
    fQATrackHistosRecEffGen        = o.fQATrackHistosRecEffGen;  
    fQATrackHistosRecEffRec        = o.fQATrackHistosRecEffRec;  
    fQATrackHistosSecRecNS         = o.fQATrackHistosSecRecNS;  
    fQATrackHistosSecRecS          = o.fQATrackHistosSecRecS;  
    fQATrackHistosSecRecSsc        = o.fQATrackHistosSecRecSsc;  
    fFFHistosRecEffRec             = o.fFFHistosRecEffRec;  
    fFFHistosSecRecNS              = o.fFFHistosSecRecNS;   
    fFFHistosSecRecS               = o.fFFHistosSecRecS;   
    fFFHistosSecRecSsc             = o.fFFHistosSecRecSsc;   
    // Background
    fh1BckgMult0                   = o.fh1BckgMult0;
    fh1BckgMult1                   = o.fh1BckgMult1;
    fh1BckgMult2                   = o.fh1BckgMult2;
    fh1BckgMult3                   = o.fh1BckgMult3;
    fh1BckgMult4                   = o.fh1BckgMult4;
    fh1FractionPtEmbedded          = o.fh1FractionPtEmbedded;
    fh1IndexEmbedded               = o.fh1IndexEmbedded;
    fh2DeltaPtVsJetPtEmbedded      = o.fh2DeltaPtVsJetPtEmbedded;
    fh2DeltaPtVsRecJetPtEmbedded   = o.fh2DeltaPtVsRecJetPtEmbedded;
    fh1DeltaREmbedded              = o.fh1DeltaREmbedded;
    fQABckgHisto0RecCuts           = o.fQABckgHisto0RecCuts;  
    fQABckgHisto0Gen               = o.fQABckgHisto0Gen;      
    fQABckgHisto1RecCuts           = o.fQABckgHisto1RecCuts;  
    fQABckgHisto1Gen               = o.fQABckgHisto1Gen;      
    fQABckgHisto2RecCuts           = o.fQABckgHisto2RecCuts;  
    fQABckgHisto2Gen               = o.fQABckgHisto2Gen;  
    fQABckgHisto3RecCuts           = o.fQABckgHisto3RecCuts;  
    fQABckgHisto3Gen               = o.fQABckgHisto3Gen;  
    fQABckgHisto4RecCuts           = o.fQABckgHisto4RecCuts;  
    fQABckgHisto4Gen               = o.fQABckgHisto4Gen;  
    fFFBckgHisto0RecCuts           = o.fFFBckgHisto0RecCuts;
    fFFBckgHisto0Gen               = o.fFFBckgHisto0Gen;       
    fFFBckgHisto1RecCuts           = o.fFFBckgHisto1RecCuts;
    fFFBckgHisto1Gen               = o.fFFBckgHisto1Gen;       
    fFFBckgHisto2RecCuts           = o.fFFBckgHisto2RecCuts;
    fFFBckgHisto2Gen               = o.fFFBckgHisto2Gen;       
    fFFBckgHisto3RecCuts           = o.fFFBckgHisto3RecCuts;
    fFFBckgHisto3Gen               = o.fFFBckgHisto3Gen;       
    fFFBckgHisto4RecCuts           = o.fFFBckgHisto4RecCuts;
    fFFBckgHisto4Gen               = o.fFFBckgHisto4Gen;       
    fFFBckgHisto0RecEffRec         = o.fFFBckgHisto0RecEffRec;
    fFFBckgHisto0SecRecNS          = o.fFFBckgHisto0SecRecNS;  
    fFFBckgHisto0SecRecS           = o.fFFBckgHisto0SecRecS;  
    fFFBckgHisto0SecRecSsc         = o.fFFBckgHisto0SecRecSsc; 
    fProNtracksLeadingJet          = o.fProNtracksLeadingJet;
    fProDelR80pcPt                 = o.fProDelR80pcPt;                       
    fProNtracksLeadingJetGen       = o.fProNtracksLeadingJetGen;             
    fProDelR80pcPtGen              = o.fProDelR80pcPtGen;                    
    fProNtracksLeadingJetBgrPerp2  = o.fProNtracksLeadingJetBgrPerp2;        
    fProNtracksLeadingJetRecPrim   = o.fProNtracksLeadingJetRecPrim;  
    fProDelR80pcPtRecPrim          = o.fProDelR80pcPtRecPrim;
    fProNtracksLeadingJetRecSecNS  = o.fProNtracksLeadingJetRecSecNS; 
    fProNtracksLeadingJetRecSecS   = o.fProNtracksLeadingJetRecSecS;  
    fProNtracksLeadingJetRecSecSsc = o.fProNtracksLeadingJetRecSecSsc;
    fRandom                        = o.fRandom;

    for(Int_t ii=0; ii<5; ii++){
      fProDelRPtSum[ii]           = o.fProDelRPtSum[ii];
      fProDelRPtSumGen[ii]        = o.fProDelRPtSumGen[ii];
      fProDelRPtSumBgrPerp2[ii]   = o.fProDelRPtSumBgrPerp2[ii];
      fProDelRPtSumRecPrim[ii]    = o.fProDelRPtSumRecPrim[ii];
      fProDelRPtSumRecSecNS[ii]   = o.fProDelRPtSumRecSecNS[ii];
      fProDelRPtSumRecSecS[ii]    = o.fProDelRPtSumRecSecS[ii];
      fProDelRPtSumRecSecSsc[ii]  = o.fProDelRPtSumRecSecSsc[ii];
    }
  }
  
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskFragmentationFunction::~AliAnalysisTaskFragmentationFunction()
{
  // destructor
  
  if(fTracksRecCuts)           delete fTracksRecCuts;
  if(fTracksGen)               delete fTracksGen;
  if(fTracksAODMCCharged)      delete fTracksAODMCCharged;  
  if(fTracksAODMCChargedSecNS) delete fTracksAODMCChargedSecNS;  
  if(fTracksAODMCChargedSecS)  delete fTracksAODMCChargedSecS;  
  if(fTracksRecQualityCuts)    delete fTracksRecQualityCuts; 
  if(fJetsRec)                 delete fJetsRec;
  if(fJetsRecCuts)             delete fJetsRecCuts;
  if(fJetsGen)                 delete fJetsGen;
  if(fJetsRecEff)              delete fJetsRecEff;
  if(fJetsEmbedded)            delete fJetsEmbedded;

  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
      fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading)){

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
							 Int_t nZ , Float_t zMin , Float_t zMax)
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
  fh2Z       = new TH2F(Form("fh2FFZ%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsZ, fZMin, fZMax);
  fh2Xi      = new TH2F(Form("fh2FFXi%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsXi, fXiMin, fXiMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{T} [GeV/c]", "entries"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2TrackPt,"jet p_{T} [GeV/c]","p_{T} [GeV/c]","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Xi,"jet p_{T} [GeV/c]","#xi", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Z,"jet p_{T} [GeV/c]","z","entries");
}

//_______________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncHistos::FillFF(Float_t trackPt, Float_t jetPt, Bool_t incrementJetPt, Float_t norm, 
									Bool_t scaleStrangeness, Float_t scaleFacStrangeness)
{
  // fill FF

  if(incrementJetPt && norm) fh1JetPt->Fill(jetPt,1/norm);
  else if(incrementJetPt) fh1JetPt->Fill(jetPt);

 // Added for proper normalization of FF background estimation
  // when zero track are found in the background region
  if((int)trackPt==-1) return;
 
  if(norm)fh2TrackPt->Fill(jetPt,trackPt,1/norm);
  else if(scaleStrangeness) fh2TrackPt->Fill(jetPt,trackPt,scaleFacStrangeness);
  else fh2TrackPt->Fill(jetPt,trackPt);
 
  Double_t z = 0.;
  if(jetPt>0) z = trackPt / jetPt;
  Double_t xi = 0;
  if(z>0) xi = TMath::Log(1/z);
  
  if(trackPt>(1-1e-06)*jetPt && trackPt<(1+1e-06)*jetPt){ // case z=1 : move entry to last histo bin <1
    z  = 1-1e-06;
    xi = 1e-06;
  }


  if(norm){
    fh2Xi->Fill(jetPt,xi,1/norm);
    fh2Z->Fill(jetPt,z,1/norm);
  }
  else if(scaleStrangeness){
    fh2Xi->Fill(jetPt,xi,scaleFacStrangeness);
    fh2Z->Fill(jetPt,z,scaleFacStrangeness);
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
    fh2PhiPt        = new TH2F(Form("fh2TrackQAPhiPt%s", fNameQAT.Data()), Form("%s: #phi - p_{T} distribution", fNameQAT.Data()), fNBinsPhi, fPhiMin, fPhiMax, fNBinsPt, fPtMin, fPtMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2HighPtEtaPhi, "#eta", "#phi");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2PhiPt, "#phi", "p_{T} [GeV/c]"); 
}

//________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::FillTrackQA(Float_t eta, Float_t phi, Float_t pt, Bool_t weightPt, Float_t norm, 
										    Bool_t scaleStrangeness, Float_t scaleFacStrangeness)
{
  // fill track QA histos
  Float_t weight = 1.;
  if(weightPt) weight = pt;  
  fh2EtaPhi->Fill( eta, phi, weight);
  if(scaleStrangeness) fh2EtaPhi->Fill( eta, phi, scaleFacStrangeness);
  if(pt > fHighPtThreshold) fh2HighPtEtaPhi->Fill( eta, phi, weight);
  if(pt > fHighPtThreshold && scaleStrangeness) fh2HighPtEtaPhi->Fill( eta, phi, weight);
  if(norm) fh1Pt->Fill( pt, 1/norm );
  else if(scaleStrangeness) fh1Pt->Fill(pt,scaleFacStrangeness); 
  else  fh1Pt->Fill( pt );

  if(scaleFacStrangeness) fh2PhiPt->Fill(phi, pt, scaleFacStrangeness);
  else fh2PhiPt->Fill(phi, pt);
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

  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE);  

  fTracksGen = new TList();
  fTracksGen->SetOwner(kFALSE);

  fTracksAODMCCharged = new TList();
  fTracksAODMCCharged->SetOwner(kFALSE);
    
  fTracksAODMCChargedSecNS = new TList();
  fTracksAODMCChargedSecNS->SetOwner(kFALSE);

  fTracksAODMCChargedSecS = new TList();
  fTracksAODMCChargedSecS->SetOwner(kFALSE);

  fTracksRecQualityCuts = new TList(); 
  fTracksRecQualityCuts->SetOwner(kFALSE);

  fJetsRec = new TList();
  fJetsRec->SetOwner(kFALSE);

  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);

  fJetsGen = new TList();
  fJetsGen->SetOwner(kFALSE);

  fJetsRecEff = new TList();
  fJetsRecEff->SetOwner(kFALSE);

  fJetsEmbedded = new TList();
  fJetsEmbedded->SetOwner(kFALSE);


  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters ||  fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
      fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading)){
    
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
  fCommonHistList->SetOwner(kTRUE);

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
  fh1nEmbeddedJets           = new TH1F("fh1nEmbeddedJets","embedded jets per event",10,-0.5,9.5);

  fh2PtRecVsGenPrim          = new TH2F("fh2PtRecVsGenPrim","rec vs gen pt",fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax,fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax);
  fh2PtRecVsGenSec           = new TH2F("fh2PtRecVsGenSec","rec vs gen pt",fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax,fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax);
  
  // embedding
  if(fBranchEmbeddedJets.Length()){
    fh1FractionPtEmbedded         = new TH1F("fh1FractionPtEmbedded","",200,0,2);
    fh1IndexEmbedded              = new TH1F("fh1IndexEmbedded","",11,-1,10);
    fh2DeltaPtVsJetPtEmbedded     = new TH2F("fh2DeltaPtVsJetPtEmbedded","",250,0,250,200,-100,100);
    fh2DeltaPtVsRecJetPtEmbedded  = new TH2F("fh2DeltaPtVsRecJetPtEmbedded","",250,0,250,200,-100,100);
    fh1DeltaREmbedded             = new TH1F("fh1DeltaREmbedded","",50,0,0.5);
    fh1nEmbeddedJets              = new TH1F("fh1nEmbeddedJets","embedded jets per event",10,-0.5,9.5);
  }


  if(fQAMode){
    if(fQAMode&1){ // track QA
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
						     fFFNBinsZ , fFFZMin , fFFZMax );

    fFFHistosGen   	     = new AliFragFuncHistos("Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  } // end: FF
  
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

      fQATrackHistosSecRecNS   = new AliFragFuncQATrackHistos("SecRecNS", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							     fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							     fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							     fQATrackHighPtThreshold);

      fQATrackHistosSecRecS    = new AliFragFuncQATrackHistos("SecRecS", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							     fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							     fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							     fQATrackHighPtThreshold);

      fQATrackHistosSecRecSsc = new AliFragFuncQATrackHistos("SecRecSsc", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);

    }
    if(fFFMode){
      fFFHistosRecEffRec      = new AliFragFuncHistos("RecEffRec", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);

      fFFHistosSecRecNS       = new AliFragFuncHistos("SecRecNS", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFHistosSecRecS        = new AliFragFuncHistos("SecRecS", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFHistosSecRecSsc      = new AliFragFuncHistos("SecRecSsc", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
    }
  } // end: efficiency

  // Background
  if(fBckgMode){
    if(fBckgType[0]==kBckgNone){
      AliError("no bgr method selected !");
    }  
    
    TString title[5];
    for(Int_t i=0; i<5; i++){
      if(fBckgType[i]==kBckgPerp) title[i]="Perp";
      else if(fBckgType[i]==kBckgPerp2) title[i]="Perp2";
      else if(fBckgType[i]==kBckgPerp2Area) title[i]="Perp2Area";
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
      else if(fBckgType[i]==kBckgNone)  title[i]="";
      else printf("Please chose background method number %d!",i);
    }


    if(fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters || 
       fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
       fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading){
      
      fh1nRecBckgJetsCuts        = new TH1F("fh1nRecBckgJetsCuts","reconstructed background jets per event",10,-0.5,9.5);
      fh1nGenBckgJets            = new TH1F("fh1nGenBckgJets","generated background jets per event",10,-0.5,9.5);
    }


    fh1BckgMult0 = new TH1F("fh1BckgMult0","bckg mult "+title[0],500,0,500);
    if(fBckgType[1] != kBckgNone) fh1BckgMult1 = new TH1F("fh1BckgMult1","bckg mult "+title[1],500,0,500);
    if(fBckgType[2] != kBckgNone) fh1BckgMult2 = new TH1F("fh1BckgMult2","bckg mult "+title[2],500,0,500);
    if(fBckgType[3] != kBckgNone) fh1BckgMult3 = new TH1F("fh1BckgMult3","bckg mult "+title[3],500,0,500);
    if(fBckgType[4] != kBckgNone) fh1BckgMult4 = new TH1F("fh1BckgMult4","bckg mult "+title[4],500,0,500);
    
    
    if(fQAMode&1){
      fQABckgHisto0RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[0]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      fQABckgHisto0Gen          = new AliFragFuncQATrackHistos("Bckg"+title[0]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							       fQATrackHighPtThreshold);
      
      if(fBckgType[1] != kBckgNone){
	fQABckgHisto1RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[1]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
	fQABckgHisto1Gen          = new AliFragFuncQATrackHistos("Bckg"+title[1]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
      }
      if(fBckgType[2] != kBckgNone){
	fQABckgHisto2RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[2]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
	fQABckgHisto2Gen          = new AliFragFuncQATrackHistos("Bckg"+title[2]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
      }
      if(fBckgType[3] != kBckgNone){
	fQABckgHisto3RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[3]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
	fQABckgHisto3Gen          = new AliFragFuncQATrackHistos("Bckg"+title[3]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
      }
      if(fBckgType[4] != kBckgNone){
	fQABckgHisto4RecCuts      = new AliFragFuncQATrackHistos("Bckg"+title[4]+"RecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
	fQABckgHisto4Gen          = new AliFragFuncQATrackHistos("Bckg"+title[4]+"Gen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								 fQATrackHighPtThreshold);
      }
    } // end: background QA
    
    if(fFFMode){
      fFFBckgHisto0RecCuts    = new AliFragFuncHistos("Bckg"+title[0]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
      
      fFFBckgHisto0Gen        = new AliFragFuncHistos("Bckg"+title[0]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
     
      if(fBckgType[1] != kBckgNone){
	fFFBckgHisto1RecCuts    = new AliFragFuncHistos("Bckg"+title[1]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
	fFFBckgHisto1Gen        = new AliFragFuncHistos("Bckg"+title[1]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
      }
      if(fBckgType[2] != kBckgNone){      
	fFFBckgHisto2RecCuts    = new AliFragFuncHistos("Bckg"+title[2]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto2Gen        = new AliFragFuncHistos("Bckg"+title[2]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
      }
      if(fBckgType[3] != kBckgNone){
	fFFBckgHisto3RecCuts    = new AliFragFuncHistos("Bckg"+title[3]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto3Gen        = new AliFragFuncHistos("Bckg"+title[3]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
      }
      if(fBckgType[4] != kBckgNone){
	fFFBckgHisto4RecCuts    = new AliFragFuncHistos("Bckg"+title[4]+"RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto4Gen        = new AliFragFuncHistos("Bckg"+title[4]+"Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
      }
      if(fEffMode){	
	fFFBckgHisto0RecEffRec      = new AliFragFuncHistos("Bckg"+title[0]+"RecEffRec", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							    fFFNBinsPt, fFFPtMin, fFFPtMax, 
							    fFFNBinsXi, fFFXiMin, fFFXiMax,  
							    fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto0SecRecNS       = new AliFragFuncHistos("Bckg"+title[0]+"SecRecNS", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							    fFFNBinsPt, fFFPtMin, fFFPtMax, 
							    fFFNBinsXi, fFFXiMin, fFFXiMax,  
							    fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto0SecRecS        = new AliFragFuncHistos("Bckg"+title[0]+"SecRecS", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							    fFFNBinsPt, fFFPtMin, fFFPtMax, 
							    fFFNBinsXi, fFFXiMin, fFFXiMax,  
							    fFFNBinsZ , fFFZMin , fFFZMax);
	
	fFFBckgHisto0SecRecSsc      = new AliFragFuncHistos("Bckg"+title[0]+"SecRecSsc", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							    fFFNBinsPt, fFFPtMin, fFFPtMax, 
							    fFFNBinsXi, fFFXiMin, fFFXiMax,  
							    fFFNBinsZ , fFFZMin , fFFZMax);

      }
    } // end: background FF


  } // end: background
  
 
  // ____________ define histograms ____________________
  
  if(fQAMode){
    if(fQAMode&1){ // track QA
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
    fFFHistosGen->DefineHistos();
  }
  
  if(fEffMode){
    if(fQAMode&1){
      fQATrackHistosRecEffGen->DefineHistos();
      fQATrackHistosRecEffRec->DefineHistos(); 
      fQATrackHistosSecRecNS->DefineHistos(); 
      fQATrackHistosSecRecS->DefineHistos(); 
      fQATrackHistosSecRecSsc->DefineHistos(); 
    }
    if(fFFMode){
      fFFHistosRecEffRec->DefineHistos();
      fFFHistosSecRecNS->DefineHistos();
      fFFHistosSecRecS->DefineHistos();
      fFFHistosSecRecSsc->DefineHistos();
    }
  } // end: efficiency

  // Background
  if(fBckgMode){
    if(fFFMode){
      fFFBckgHisto0RecCuts->DefineHistos();
      fFFBckgHisto0Gen->DefineHistos();      
      if(fBckgType[1] != kBckgNone) fFFBckgHisto1RecCuts->DefineHistos();
      if(fBckgType[1] != kBckgNone) fFFBckgHisto1Gen->DefineHistos();
      if(fBckgType[2] != kBckgNone) fFFBckgHisto2RecCuts->DefineHistos();
      if(fBckgType[2] != kBckgNone) fFFBckgHisto2Gen->DefineHistos();
      if(fBckgType[3] != kBckgNone) fFFBckgHisto3RecCuts->DefineHistos();
      if(fBckgType[3] != kBckgNone) fFFBckgHisto3Gen->DefineHistos();
      if(fBckgType[4] != kBckgNone) fFFBckgHisto4RecCuts->DefineHistos();
      if(fBckgType[4] != kBckgNone) fFFBckgHisto4Gen->DefineHistos();

     if(fEffMode){
	fFFBckgHisto0RecEffRec->DefineHistos();	
	fFFBckgHisto0SecRecNS->DefineHistos();
        fFFBckgHisto0SecRecS->DefineHistos();
        fFFBckgHisto0SecRecSsc->DefineHistos();
      }
    }

    if(fQAMode&1){
      fQABckgHisto0RecCuts->DefineHistos();
      fQABckgHisto0Gen->DefineHistos();
      if(fBckgType[1] != kBckgNone) fQABckgHisto1RecCuts->DefineHistos();
      if(fBckgType[1] != kBckgNone) fQABckgHisto1Gen->DefineHistos();
      if(fBckgType[2] != kBckgNone) fQABckgHisto2RecCuts->DefineHistos();
      if(fBckgType[2] != kBckgNone) fQABckgHisto2Gen->DefineHistos();
      if(fBckgType[3] != kBckgNone) fQABckgHisto3RecCuts->DefineHistos();
      if(fBckgType[3] != kBckgNone) fQABckgHisto3Gen->DefineHistos();
      if(fBckgType[4] != kBckgNone) fQABckgHisto4RecCuts->DefineHistos();
      if(fBckgType[4] != kBckgNone) fQABckgHisto4Gen->DefineHistos();
    }
  } // end: background
  

  Bool_t genJets    = (fJetTypeGen != kJetsUndef) ? kTRUE : kFALSE;
  Bool_t genTracks  = (fTrackTypeGen != kTrackUndef) ? kTRUE : kFALSE;
  Bool_t recJetsEff = (fJetTypeRecEff != kJetsUndef) ? kTRUE : kFALSE;

  fCommonHistList->Add(fh1EvtSelection);
  fCommonHistList->Add(fh1EvtMult);
  fCommonHistList->Add(fh1EvtCent);
  fCommonHistList->Add(fh1VertexNContributors);
  fCommonHistList->Add(fh1VertexZ);    
  fCommonHistList->Add(fh1nRecJetsCuts);
  fCommonHistList->Add(fh1Xsec);
  fCommonHistList->Add(fh1Trials);
  fCommonHistList->Add(fh1PtHard);
  fCommonHistList->Add(fh1PtHardTrials);
 
  if(genJets) fCommonHistList->Add(fh1nGenJets);

  // FF histograms
  if(fFFMode){
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    if(genJets && genTracks){
      fFFHistosGen->AddToOutput(fCommonHistList);
    }
  }

  // Background
  if(fBckgMode){
    if(fFFMode){
      fFFBckgHisto0RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[1] != kBckgNone) fFFBckgHisto1RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[2] != kBckgNone) fFFBckgHisto2RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[3] != kBckgNone) fFFBckgHisto3RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[4] != kBckgNone) fFFBckgHisto4RecCuts->AddToOutput(fCommonHistList);

      if(genJets && genTracks){
	fFFBckgHisto0Gen->AddToOutput(fCommonHistList);
	if(fBckgType[1] != kBckgNone) fFFBckgHisto1Gen->AddToOutput(fCommonHistList);
	if(fBckgType[2] != kBckgNone) fFFBckgHisto2Gen->AddToOutput(fCommonHistList);
	if(fBckgType[3] != kBckgNone) fFFBckgHisto3Gen->AddToOutput(fCommonHistList);
	if(fBckgType[4] != kBckgNone) fFFBckgHisto4Gen->AddToOutput(fCommonHistList);
      }

      if(fEffMode){
        fFFBckgHisto0RecEffRec->AddToOutput(fCommonHistList);
	fFFBckgHisto0SecRecNS->AddToOutput(fCommonHistList);
        fFFBckgHisto0SecRecS->AddToOutput(fCommonHistList);
        fFFBckgHisto0SecRecSsc->AddToOutput(fCommonHistList);
      }
    }

    if(fQAMode&1){
      fQABckgHisto0RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[1] != kBckgNone) fQABckgHisto1RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[2] != kBckgNone) fQABckgHisto2RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[3] != kBckgNone) fQABckgHisto3RecCuts->AddToOutput(fCommonHistList);
      if(fBckgType[4] != kBckgNone) fQABckgHisto4RecCuts->AddToOutput(fCommonHistList);
      if(genJets && genTracks){
	fQABckgHisto0Gen->AddToOutput(fCommonHistList);
	if(fBckgType[1] != kBckgNone) fQABckgHisto1Gen->AddToOutput(fCommonHistList);
	if(fBckgType[2] != kBckgNone) fQABckgHisto2Gen->AddToOutput(fCommonHistList);
	if(fBckgType[3] != kBckgNone) fQABckgHisto3Gen->AddToOutput(fCommonHistList);
	if(fBckgType[4] != kBckgNone) fQABckgHisto4Gen->AddToOutput(fCommonHistList);
      }
    }
    
    if(fh1BckgMult0) fCommonHistList->Add(fh1BckgMult0);
    if(fBckgType[1] != kBckgNone)  fCommonHistList->Add(fh1BckgMult1);
    if(fBckgType[2] != kBckgNone)  fCommonHistList->Add(fh1BckgMult2);
    if(fBckgType[3] != kBckgNone)  fCommonHistList->Add(fh1BckgMult3);
    if(fBckgType[4] != kBckgNone)  fCommonHistList->Add(fh1BckgMult4);
  }
  

  if(fBranchEmbeddedJets.Length()){ 
    fCommonHistList->Add(fh1FractionPtEmbedded);
    fCommonHistList->Add(fh1IndexEmbedded);
    fCommonHistList->Add(fh2DeltaPtVsJetPtEmbedded);      
    fCommonHistList->Add(fh2DeltaPtVsRecJetPtEmbedded);      
    fCommonHistList->Add(fh1DeltaREmbedded);
    fCommonHistList->Add(fh1nEmbeddedJets);  
  }


  // QA  
  if(fQAMode){
    if(fQAMode&1){ // track QA
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
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
      fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading)) {
    fCommonHistList->Add(fh1nRecBckgJetsCuts);
    if(genJets) fCommonHistList->Add(fh1nGenBckgJets);
  }
    
   
  if(fEffMode && recJetsEff && genTracks){
    if(fQAMode&1){
      fQATrackHistosRecEffGen->AddToOutput(fCommonHistList);
      fQATrackHistosRecEffRec->AddToOutput(fCommonHistList);
      fQATrackHistosSecRecNS->AddToOutput(fCommonHistList);
      fQATrackHistosSecRecS->AddToOutput(fCommonHistList);
      fQATrackHistosSecRecSsc->AddToOutput(fCommonHistList);
    }
    if(fFFMode){
      fFFHistosRecEffRec->AddToOutput(fCommonHistList);
      fFFHistosSecRecNS->AddToOutput(fCommonHistList);
      fFFHistosSecRecS->AddToOutput(fCommonHistList);
      fFFHistosSecRecSsc->AddToOutput(fCommonHistList);
    }
    fCommonHistList->Add(fh1nRecEffJets);
    fCommonHistList->Add(fh2PtRecVsGenPrim); 
    fCommonHistList->Add(fh2PtRecVsGenSec); 
  }
  
  // jet shape 
  if(fJSMode){

    fProNtracksLeadingJet          = new TProfile("AvgNoOfTracksLeadingJet","AvgNoOfTracksLeadingJet",100,0,250,0,50); 
    fProDelR80pcPt                 = new TProfile("AvgdelR80pcPt","AvgdelR80pcPt",100,0,250,0,1); 

    if(genJets && genTracks){
      fProNtracksLeadingJetGen       = new TProfile("AvgNoOfTracksLeadingJetGen","AvgNoOfTracksLeadingJetGen",100,0,250,0,50); 
      fProDelR80pcPtGen              = new TProfile("AvgdelR80pcPtGen","AvgdelR80pcPtGen",100,0,250,0,1); 
    }

    if(fBckgMode)
      fProNtracksLeadingJetBgrPerp2  = new TProfile("AvgNoOfTracksLeadingJetBgrPerp2","AvgNoOfTracksLeadingJetBgrPerp2",100,0,250,0,50); 
    
    if(fEffMode){
      fProNtracksLeadingJetRecPrim   = new TProfile("AvgNoOfTracksLeadingJetRecPrim","AvgNoOfTracksLeadingJetRecPrim",100,0,250,0,50); 
      fProDelR80pcPtRecPrim          = new TProfile("AvgdelR80pcPtRecPrim","AvgdelR80pcPtRecPrim",100,0,250,0,1); 
      fProNtracksLeadingJetRecSecNS  = new TProfile("AvgNoOfTracksLeadingJetRecSecNS","AvgNoOfTracksLeadingJetRecSecNS",100,0,250,0,50); 
      fProNtracksLeadingJetRecSecS   = new TProfile("AvgNoOfTracksLeadingJetRecSecS","AvgNoOfTracksLeadingJetRecSecS",100,0,250,0,50); 
      fProNtracksLeadingJetRecSecSsc = new TProfile("AvgNoOfTracksLeadingJetRecSecSsc","AvgNoOfTracksLeadingJetRecSecSsc",100,0,250,0,50); 
    }

    TString strTitJS;	
    for(Int_t ii=0; ii<5; ii++){
      if(ii==0)strTitJS = "_JetPt20to30";
      if(ii==1)strTitJS = "_JetPt30to40";
      if(ii==2)strTitJS = "_JetPt40to60";
      if(ii==3)strTitJS = "_JetPt60to80";
      if(ii==4)strTitJS = "_JetPt80to100";
      
      fProDelRPtSum[ii]            = new TProfile(Form("AvgPtSumDelR%s",strTitJS.Data()),Form("AvgPtSumDelR%s",strTitJS.Data()),100,0,1,0,250);
      if(genJets && genTracks) 
	fProDelRPtSumGen[ii]       = new TProfile(Form("AvgPtSumDelRGen%s",strTitJS.Data()),Form("AvgPtSumDelRGen%s",strTitJS.Data()),100,0,1,0,250);
      if(fBckgMode) 
	fProDelRPtSumBgrPerp2[ii]  = new TProfile(Form("AvgPtSumDelRBgrPerp2%s",strTitJS.Data()),Form("AvgPtSumDelRBgrPerp2%s",strTitJS.Data()),100,0,1,0,250);
      if(fEffMode){
	fProDelRPtSumRecPrim[ii]   = new TProfile(Form("AvgPtSumDelRRecPrim%s",strTitJS.Data()),Form("AvgPtSumDelRRecPrim%s",strTitJS.Data()),100,0,1,0,250);
	fProDelRPtSumRecSecNS[ii]  = new TProfile(Form("AvgPtSumDelRRecSecNS%s",strTitJS.Data()),Form("AvgPtSumDelRRecSecNS%s",strTitJS.Data()),100,0,1,0,250);
	fProDelRPtSumRecSecS[ii]   = new TProfile(Form("AvgPtSumDelRRecSecS%s",strTitJS.Data()),Form("AvgPtSumDelRRecSecS%s",strTitJS.Data()),100,0,1,0,250);
	fProDelRPtSumRecSecSsc[ii] = new TProfile(Form("AvgPtSumDelRRecSecSsc%s",strTitJS.Data()),Form("AvgPtSumDelRRecSecSsc%s",strTitJS.Data()),100,0,1,0,250);
      }
    }
    
    fCommonHistList->Add(fProNtracksLeadingJet);
    fCommonHistList->Add(fProDelR80pcPt);
    for(int ii=0; ii<5; ii++) fCommonHistList->Add(fProDelRPtSum[ii]);

    if(genJets && genTracks){
      fCommonHistList->Add(fProNtracksLeadingJetGen);
      fCommonHistList->Add(fProDelR80pcPtGen);
      for(Int_t ii=0; ii<5; ii++)  fCommonHistList->Add(fProDelRPtSumGen[ii]);
    }
    
    if(fBckgMode){ 
      fCommonHistList->Add(fProNtracksLeadingJetBgrPerp2);
      for(Int_t ii=0; ii<5; ii++) fCommonHistList->Add(fProDelRPtSumBgrPerp2[ii]);
    }

    if(fEffMode){
      fCommonHistList->Add(fProNtracksLeadingJetRecPrim);
      fCommonHistList->Add(fProDelR80pcPtRecPrim);
      for(Int_t ii=0; ii<5; ii++)  fCommonHistList->Add(fProDelRPtSumRecPrim[ii]);
      
      fCommonHistList->Add(fProNtracksLeadingJetRecSecNS);
      for(Int_t ii=0; ii<5; ii++)  fCommonHistList->Add(fProDelRPtSumRecSecNS[ii]);

      fCommonHistList->Add(fProNtracksLeadingJetRecSecS);
      for(Int_t ii=0; ii<5; ii++)  fCommonHistList->Add(fProDelRPtSumRecSecS[ii]);
      
      fCommonHistList->Add(fProNtracksLeadingJetRecSecSsc);
      for(Int_t ii=0; ii<5; ii++)  fCommonHistList->Add(fProDelRPtSumRecSecSsc[ii]);
    }
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

  PostData(1, fCommonHistList);
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
  
  if(!(inputHandler->IsEventSelected() & fEvtSelectionMask)){
    fh1EvtSelection->Fill(1.);
    if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
    PostData(1, fCommonHistList);
    return;
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
    if(fUseAODInputJets) fAODJets = fAOD;
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      fAODJets = fAOD;
      if (fDebug > 1)  Printf("%s:%d AOD event from output", (char*)__FILE__,__LINE__);
    }
  }
  
  if(!fAODJets && !fUseAODInputJets){ // case we have AOD in input & output and want jets from output
    TObject* outHandler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( outHandler && outHandler->InheritsFrom("AliAODHandler") ) {
      fAODJets = ((AliAODHandler*)outHandler)->GetAOD();
      if (fDebug > 1)  Printf("%s:%d jets from output AOD", (char*)__FILE__,__LINE__);
    }
  }
  
  if(fNonStdFile.Length()!=0){
    // case we have an AOD extension - fetch the jets from the extended output
    
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);    
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension not found for %s",fNonStdFile.Data());
    }
  }
  
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return;
  }
  if(!fAODJets){
    Printf("%s:%d AODEvent with jet branch not found", (char*)__FILE__,__LINE__);
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
      if(fESD) centPercent = fESD->GetCentrality()->GetCentralityPercentile("V0M"); // retrieve value 'by hand'
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

  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 

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
      
      //fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
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
  
  
  Int_t nEmbeddedJets =  0; 
  TArrayI iEmbeddedMatchIndex; 
  TArrayF fEmbeddedPtFraction; 
  

  if(fBranchEmbeddedJets.Length()){ 
    Int_t nJEmbedded = GetListOfJets(fJetsEmbedded, kJetsEmbedded);
    if(nJEmbedded>=0) nEmbeddedJets = fJetsEmbedded->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Embedded jets: %d %d",(char*)__FILE__,__LINE__,nJEmbedded,nEmbeddedJets);
    if(nJEmbedded != nEmbeddedJets) Printf("%s:%d Mismatch Selected Embedded Jets: %d %d",(char*)__FILE__,__LINE__,nJEmbedded,nEmbeddedJets);
    fh1nEmbeddedJets->Fill(nEmbeddedJets);
    
    Float_t maxDist = 0.3;

    iEmbeddedMatchIndex.Set(nEmbeddedJets); 
    fEmbeddedPtFraction.Set(nEmbeddedJets); 
    
    iEmbeddedMatchIndex.Reset(-1);
    fEmbeddedPtFraction.Reset(0);
    
    AliAnalysisHelperJetTasks::GetJetMatching(fJetsEmbedded, nEmbeddedJets, 
					      fJetsRecCuts, nRecJetsCuts, 
					      iEmbeddedMatchIndex, fEmbeddedPtFraction,
					      fDebug, maxDist);
    
  }
  
  //____ fetch background clusters ___________________________________________________
  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
      fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading)){

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
  
  Int_t nTCuts;
  if(fUseExtraTracks ==  1)      nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODExtraCuts);
  else if(fUseExtraTracks == -1) nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODExtraonlyCuts);
  else                           nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);
  
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
      for(Int_t it=0; it<nRecPartCuts; ++it){
	AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(it));
	if(part)fQATrackHistosRecCuts->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt() );
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
  
  if(fQAMode || fFFMode){
    for(Int_t ij=0; ij<nRecJetsCuts; ++ij){
      
      AliAODJet* jet = (AliAODJet*)(fJetsRecCuts->At(ij));
      if(fQAMode&2) fQAJetHistosRecCuts->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      
      if(ij==0){ // leading jet
	
	if(fQAMode&2) fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );

	Double_t ptFractionEmbedded = 0; 
	AliAODJet* embeddedJet = 0; 

	if(fBranchEmbeddedJets.Length()){ // find embedded jet

	  Int_t indexEmbedded = -1;
	  for(Int_t i=0; i<nEmbeddedJets; i++){
	    if(iEmbeddedMatchIndex[i] == ij){
	      indexEmbedded      = i;
	      ptFractionEmbedded = fEmbeddedPtFraction[i];
	    }
	  }

	  fh1IndexEmbedded->Fill(indexEmbedded);
	  fh1FractionPtEmbedded->Fill(ptFractionEmbedded);
	  
	  if(indexEmbedded>-1){ 
	    
	    embeddedJet = dynamic_cast<AliAODJet*>(fJetsEmbedded->At(indexEmbedded));
	    if(!embeddedJet) continue;

	    Double_t deltaPt = jet->Pt() - embeddedJet->Pt();
	    Double_t deltaR  = jet->DeltaR((AliVParticle*) (embeddedJet)); 
	    
	    fh2DeltaPtVsJetPtEmbedded->Fill(embeddedJet->Pt(),deltaPt);
	    fh2DeltaPtVsRecJetPtEmbedded->Fill(jet->Pt(),deltaPt);
	    fh1DeltaREmbedded->Fill(deltaR);
	  }
	}

	// get tracks in jet
	TList* jettracklist = new TList();
	Double_t sumPt      = 0.;
	Bool_t isBadJet     = kFALSE;

	if(GetFFRadius()<=0){
	  GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	} else {
	  if(fUseEmbeddedJetAxis){
	    if(embeddedJet) GetJetTracksPointing(fTracksRecCuts, jettracklist, embeddedJet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	  }
	  else              GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	}
	
	if(GetFFMinNTracks()>0 && jettracklist->GetSize()<=GetFFMinNTracks()) isBadJet = kTRUE;
	
	if(isBadJet) continue; 

	if(ptFractionEmbedded>=fCutFractionPtEmbedded){ // if no embedding: ptFraction = cutFraction = 0
	  
	  for(Int_t it=0; it<jettracklist->GetSize(); ++it){

	    AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	    if(!trackVP)continue;

            AliAODTrack * aodtrack  = dynamic_cast<AliAODTrack*>(jettracklist->At(it));
            if(!aodtrack) continue;

	    TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	    
	    Float_t jetPt   = jet->Pt();
	    if(fUseEmbeddedJetPt){
	      if(embeddedJet) jetPt = embeddedJet->Pt();
	      else jetPt = 0;
	    }
	    Float_t trackPt = trackV->Pt();


	    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    
	    if(fFFMode) fFFHistosRecCuts->FillFF(trackPt, jetPt, incrementJetPt);

	    delete trackV;
	  }
	  
	  // background ff
	  if(fBckgMode){
	    if(fBckgType[0]!=kBckgNone)
	      FillBckgHistos(fBckgType[0], fTracksRecCuts, fJetsRecCuts, jet,  
			     fFFBckgHisto0RecCuts,fQABckgHisto0RecCuts, fh1BckgMult0);
	    if(fBckgType[1]!=kBckgNone)
	      FillBckgHistos(fBckgType[1], fTracksRecCuts, fJetsRecCuts, jet,
			     fFFBckgHisto1RecCuts,fQABckgHisto1RecCuts, fh1BckgMult1);
	    if(fBckgType[2]!=kBckgNone)
	      FillBckgHistos(fBckgType[2], fTracksRecCuts, fJetsRecCuts, jet, 
			     fFFBckgHisto2RecCuts,fQABckgHisto2RecCuts, fh1BckgMult2);
	    if(fBckgType[3]!=kBckgNone)
	      FillBckgHistos(fBckgType[3], fTracksRecCuts, fJetsRecCuts, jet, 
			     fFFBckgHisto3RecCuts,fQABckgHisto3RecCuts, fh1BckgMult3);
	    if(fBckgType[4]!=kBckgNone)
	      FillBckgHistos(fBckgType[4], fTracksRecCuts, fJetsRecCuts, jet, 
			     fFFBckgHisto4RecCuts,fQABckgHisto4RecCuts, fh1BckgMult4);
	  } // end if(fBckgMode)
	 

	  if(fJSMode) FillJetShape(jet, jettracklist, fProNtracksLeadingJet, fProDelRPtSum, fProDelR80pcPt);
	   
	  delete jettracklist;	

	} // end: cut embedded ratio
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
	Bool_t isBadJet     = kFALSE;
	
	if(GetFFRadius()<=0){
	  GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	} else {
	  GetJetTracksPointing(fTracksGen, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	}
	
	if(GetFFMinNTracks()>0 && jettracklist->GetSize()<=GetFFMinNTracks()) isBadJet = kTRUE;;
	if(isBadJet) continue; 

	for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	  
	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	  if(!trackVP)continue;
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	  
	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  if(fFFMode) fFFHistosGen->FillFF( trackPt, jetPt, incrementJetPt );
	  
	  delete trackV;
	}

	if(fBckgMode){
	  if(fBckgType[0]!=kBckgNone)
	    FillBckgHistos(fBckgType[0], fTracksGen, fJetsGen, jet,
			   fFFBckgHisto0Gen, fQABckgHisto0Gen);
	  if(fBckgType[1]!=kBckgNone)
	    FillBckgHistos(fBckgType[1], fTracksGen, fJetsGen, jet,
			   fFFBckgHisto1Gen, fQABckgHisto1Gen);
	  if(fBckgType[2]!=kBckgNone)
	    FillBckgHistos(fBckgType[2], fTracksGen, fJetsGen, jet,
			   fFFBckgHisto2Gen, fQABckgHisto2Gen);
	  if(fBckgType[3]!=kBckgNone)
	    FillBckgHistos(fBckgType[3], fTracksGen, fJetsGen, jet,
			   fFFBckgHisto3Gen, fQABckgHisto3Gen);
	  if(fBckgType[4]!=kBckgNone)
	    FillBckgHistos(fBckgType[4], fTracksGen, fJetsGen, jet,
			   fFFBckgHisto4Gen, fQABckgHisto4Gen);
	} // end if(fBckgMode)
	

	if(fJSMode) FillJetShape(jet, jettracklist, fProNtracksLeadingJetGen, fProDelRPtSumGen, fProDelR80pcPtGen);

	delete jettracklist;
      }
    }
  } // end: QA, FF and intra-jet

      
  // ____ efficiency _______________________________

  if(fEffMode && (fJetTypeRecEff != kJetsUndef)){

    // arrays holding for each generated particle the reconstructed AOD track index & isPrimary flag, are initialized in AssociateGenRec(...) function
    TArrayI indexAODTr; 
    TArrayS isGenPrim; 

    // array holding for each reconstructed AOD track generated particle index, initialized in AssociateGenRec(...) function
    TArrayI indexMCTr; 

    // ... and another set for secondaries from strange/non strange mothers (secondary MC tracks are stored in different lists)
    TArrayI indexAODTrSecNS; 
    TArrayS isGenSecNS; 
    TArrayI indexMCTrSecNS; 
   
    TArrayI indexAODTrSecS; 
    TArrayS isGenSecS; 
    TArrayI indexMCTrSecS; 

    Int_t  nTracksAODMCCharged = GetListOfTracks(fTracksAODMCCharged, kTrackAODMCCharged);
    if(fDebug>2)Printf("%s:%d selected AODMC tracks: %d ",(char*)__FILE__,__LINE__,nTracksAODMCCharged);
  
    Int_t  nTracksAODMCChargedSecNS = GetListOfTracks(fTracksAODMCChargedSecNS, kTrackAODMCChargedSecNS);
    if(fDebug>2)Printf("%s:%d selected AODMC secondary tracks NS: %d ",(char*)__FILE__,__LINE__,nTracksAODMCChargedSecNS);
  
    Int_t  nTracksAODMCChargedSecS = GetListOfTracks(fTracksAODMCChargedSecS, kTrackAODMCChargedSecS);
    if(fDebug>2)Printf("%s:%d selected AODMC secondary tracks S: %d ",(char*)__FILE__,__LINE__,nTracksAODMCChargedSecS);

    Int_t  nTracksRecQualityCuts = GetListOfTracks(fTracksRecQualityCuts, kTrackAODQualityCuts);
    if(fDebug>2)Printf("%s:%d selected rec tracks quality after cuts, full acceptance/pt : %d ",(char*)__FILE__,__LINE__,nTracksRecQualityCuts);
  
    // associate gen and rec tracks, store indices in TArrays 
    AssociateGenRec(fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,indexMCTr,isGenPrim,fh2PtRecVsGenPrim); 
    AssociateGenRec(fTracksAODMCChargedSecNS,fTracksRecQualityCuts,indexAODTrSecNS,indexMCTrSecNS,isGenSecNS,fh2PtRecVsGenSec);
    AssociateGenRec(fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,indexMCTrSecS,isGenSecS,fh2PtRecVsGenSec);
  
    // single track eff
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGen,fQATrackHistosRecEffRec,fTracksAODMCCharged,indexAODTr,isGenPrim);

    // secondaries
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecNS,fTracksAODMCChargedSecNS,indexAODTrSecNS,isGenSecNS);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecS,fTracksAODMCChargedSecS,indexAODTrSecS,isGenSecS);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecSsc,fTracksAODMCChargedSecS,indexAODTrSecS,isGenSecS,kTRUE);


    // jet track eff    
    Double_t sumPtGenLeadingJetRecEff = 0;
    Double_t sumPtGenLeadingJetSec    = 0;
    Double_t sumPtRecLeadingJetRecEff = 0;
    
    for(Int_t ij=0; ij<nRecEffJets; ++ij){ // jet loop 
    
      AliAODJet* jet = (AliAODJet*)(fJetsRecEff->At(ij));

      Bool_t isBadJetGenPrim = kFALSE;
      Bool_t isBadJetGenSec  = kFALSE;
      Bool_t isBadJetRec     = kFALSE;
    

      if(ij==0){ // leading jet
	
	// for efficiency: gen tracks from pointing with gen/rec jet
	TList* jettracklistGenPrim = new TList();
	
	// if radius<0 -> trackRefs: collect gen tracks in wide radius + fill FF recEff rec histos with tracks contained in track refs
        // note : FF recEff gen histos will be somewhat useless in this approach

	if(GetFFRadius() >0)
	  GetJetTracksPointing(fTracksAODMCCharged, jettracklistGenPrim, jet, GetFFRadius(), sumPtGenLeadingJetRecEff, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetGenPrim); 
	else
	  GetJetTracksPointing(fTracksAODMCCharged, jettracklistGenPrim, jet, TMath::Abs(GetFFRadius())+0.2, sumPtGenLeadingJetRecEff, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetGenPrim); 

	TList* jettracklistGenSecNS = new TList();
	if(GetFFRadius() >0)
	  GetJetTracksPointing(fTracksAODMCChargedSecNS, jettracklistGenSecNS, jet, GetFFRadius(), sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 
        else
	  GetJetTracksPointing(fTracksAODMCChargedSecNS, jettracklistGenSecNS, jet, TMath::Abs(GetFFRadius())+0.2, sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 

	TList* jettracklistGenSecS = new TList();
	if(GetFFRadius() >0)
	  GetJetTracksPointing(fTracksAODMCChargedSecS, jettracklistGenSecS, jet, GetFFRadius(), sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 
        else
	  GetJetTracksPointing(fTracksAODMCChargedSecS, jettracklistGenSecS, jet, TMath::Abs(GetFFRadius())+0.2, sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 


	// bin efficiency in jet pt bins using rec tracks  
	TList* jettracklistRec = new TList();
	if(GetFFRadius() >0) GetJetTracksPointing(fTracksRecCuts,jettracklistRec, jet, GetFFRadius(), sumPtRecLeadingJetRecEff, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetRec); 
	else                 GetJetTracksTrackrefs(jettracklistRec, jet, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetRec); 
	

	Double_t jetEta   = jet->Eta();
	Double_t jetPhi   = TVector2::Phi_0_2pi(jet->Phi());
	
	if(GetFFMinNTracks()>0 && jettracklistGenPrim->GetSize()<=GetFFMinNTracks())   isBadJetGenPrim = kTRUE;
	if(GetFFMinNTracks()>0 && jettracklistGenSecNS->GetSize()<=GetFFMinNTracks())  isBadJetGenSec  = kTRUE;
	if(GetFFMinNTracks()>0 && jettracklistRec->GetSize()<=GetFFMinNTracks())       isBadJetRec     = kTRUE;

	if(isBadJetRec) continue;

	if(fQAMode&2) fQAJetHistosRecEffLeading->FillJetQA( jetEta, jetPhi, sumPtGenLeadingJetRecEff ); 
	
	if(fFFMode){
	  
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRec,jet,
						    jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim, 
						    0,kFALSE,fJSMode,fProNtracksLeadingJetRecPrim,fProDelRPtSumRecPrim,fProDelR80pcPtRecPrim); 
	  
          else                FillJetTrackHistosRec(fFFHistosRecEffRec,jet,
						    jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,
						    jettracklistRec,kFALSE,fJSMode,fProNtracksLeadingJetRecPrim,fProDelRPtSumRecPrim,fProDelR80pcPtRecPrim); 
	  

	  // secondaries: use jet pt from primaries 
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecNS,jet,
						    jettracklistGenSecNS,fTracksAODMCChargedSecNS,fTracksRecQualityCuts, indexAODTrSecNS,isGenSecNS,
						    0,kFALSE,fJSMode,fProNtracksLeadingJetRecSecNS,fProDelRPtSumRecSecNS); 
	  
	  else                FillJetTrackHistosRec(fFFHistosSecRecNS,jet,
						    jettracklistGenSecNS,fTracksAODMCChargedSecNS,fTracksRecQualityCuts,indexAODTrSecNS,isGenSecNS,
						    jettracklistRec,kFALSE,fJSMode,fProNtracksLeadingJetRecSecNS,fProDelRPtSumRecSecNS);  
	  
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecS,jet,
						    jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,
						    0,kFALSE,fJSMode,fProNtracksLeadingJetRecSecS,fProDelRPtSumRecSecS); 

 	  else                FillJetTrackHistosRec(fFFHistosSecRecS,jet,
						    jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,
						    jettracklistRec,kFALSE,fJSMode,fProNtracksLeadingJetRecSecS,fProDelRPtSumRecSecS);  
	  
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecSsc,jet,
						    jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,
						    0,kTRUE,fJSMode,fProNtracksLeadingJetRecSecSsc,fProDelRPtSumRecSecSsc); 
	  
	  else                FillJetTrackHistosRec(fFFHistosSecRecSsc,jet,
						    jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,
						    jettracklistRec,kTRUE,fJSMode,fProNtracksLeadingJetRecSecSsc,fProDelRPtSumRecSecSsc);
        }
	
	delete jettracklistGenPrim;
	delete jettracklistGenSecNS;
	delete jettracklistGenSecS;
	delete jettracklistRec;
      
      
        if(fBckgMode && fFFMode){ 

	  TList* perpjettracklistGen  = new TList();
	  TList* perpjettracklistGen1 = new TList();
	  TList* perpjettracklistGen2 = new TList();

	  Double_t sumPtGenPerp  = 0.;
	  Double_t sumPtGenPerp1 = 0.;
	  Double_t sumPtGenPerp2 = 0.;
	  GetTracksTiltedwrpJetAxis(TMath::Pi()/2.,fTracksAODMCCharged, perpjettracklistGen1, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerp1); 
	  GetTracksTiltedwrpJetAxis(-1*TMath::Pi()/2.,fTracksAODMCCharged, perpjettracklistGen2, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerp2); 

	  perpjettracklistGen->AddAll(perpjettracklistGen1);
	  perpjettracklistGen->AddAll(perpjettracklistGen2);
	  sumPtGenPerp = 0.5*(sumPtGenPerp1+sumPtGenPerp2);

	  TList* perpjettracklistGenSecNS  = new TList();
	  TList* perpjettracklistGenSecNS1 = new TList();
	  TList* perpjettracklistGenSecNS2 = new TList();

          Double_t sumPtGenPerpNS;
          Double_t sumPtGenPerpNS1;
          Double_t sumPtGenPerpNS2;
          GetTracksTiltedwrpJetAxis(TMath::Pi()/2.,fTracksAODMCChargedSecNS, perpjettracklistGenSecNS1, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerpNS1); 
          GetTracksTiltedwrpJetAxis(-1*TMath::Pi()/2.,fTracksAODMCChargedSecNS, perpjettracklistGenSecNS2, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerpNS2); 

	  perpjettracklistGenSecNS->AddAll(perpjettracklistGenSecNS1);
	  perpjettracklistGenSecNS->AddAll(perpjettracklistGenSecNS2);
	  sumPtGenPerpNS = 0.5*(sumPtGenPerpNS1+sumPtGenPerpNS2);


	  TList* perpjettracklistGenSecS  = new TList();
	  TList* perpjettracklistGenSecS1 = new TList();
	  TList* perpjettracklistGenSecS2 = new TList();

          Double_t sumPtGenPerpS;
          Double_t sumPtGenPerpS1;
          Double_t sumPtGenPerpS2;
          GetTracksTiltedwrpJetAxis(TMath::Pi()/2.,fTracksAODMCChargedSecS, perpjettracklistGenSecS1, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerpS1); 
          GetTracksTiltedwrpJetAxis(-1*TMath::Pi()/2.,fTracksAODMCChargedSecS, perpjettracklistGenSecS2, jet, TMath::Abs(GetFFRadius()) , sumPtGenPerpS2); 

	  perpjettracklistGenSecS->AddAll(perpjettracklistGenSecS1);
	  perpjettracklistGenSecS->AddAll(perpjettracklistGenSecS2);
	  sumPtGenPerpS = 0.5*(sumPtGenPerpS1+sumPtGenPerpS2);


          if(perpjettracklistGen->GetSize() != perpjettracklistGen1->GetSize() + perpjettracklistGen2->GetSize()){
	    cout<<" ERROR: perpjettracklistGen size "<<perpjettracklistGen->GetSize()<<" perp1 "<<perpjettracklistGen1->GetSize()
	        <<" perp2 "<<perpjettracklistGen2->GetSize()<<endl;
            exit(0); 
          }

          if(perpjettracklistGenSecNS->GetSize() != perpjettracklistGenSecNS1->GetSize() + perpjettracklistGenSecNS2->GetSize()){
	    cout<<" ERROR: perpjettracklistGenSecNS size "<<perpjettracklistGenSecNS->GetSize()<<" perp1 "<<perpjettracklistGenSecNS1->GetSize()
	        <<" perp2 "<<perpjettracklistGenSecNS2->GetSize()<<endl;
            exit(0); 
          }

          if(perpjettracklistGenSecS->GetSize() != perpjettracklistGenSecS1->GetSize() + perpjettracklistGenSecS2->GetSize()){
	    cout<<" ERROR: perpjettracklistGenSecS size "<<perpjettracklistGenSecS->GetSize()<<" perp1 "<<perpjettracklistGenSecS1->GetSize()
	        <<" perp2 "<<perpjettracklistGenSecS2->GetSize()<<endl;
            exit(0); 
          }


	  FillJetTrackHistosRec(fFFBckgHisto0RecEffRec,jet,
				perpjettracklistGen,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim); 
	  
	  FillJetTrackHistosRec(fFFBckgHisto0SecRecNS,jet,
				perpjettracklistGenSecNS,fTracksAODMCChargedSecNS,fTracksRecQualityCuts,indexAODTrSecNS,isGenSecNS); 
	  
	  FillJetTrackHistosRec(fFFBckgHisto0SecRecS,jet,
				perpjettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS); 
	  
	  FillJetTrackHistosRec(fFFBckgHisto0SecRecSsc,jet,
				perpjettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,0,kTRUE); 
	  
	  
	  delete perpjettracklistGen;
	  delete perpjettracklistGen1;
	  delete perpjettracklistGen2;

	  delete perpjettracklistGenSecNS;
	  delete perpjettracklistGenSecNS1;
	  delete perpjettracklistGenSecNS2;

	  delete perpjettracklistGenSecS;
	  delete perpjettracklistGenSecS1;
	  delete perpjettracklistGenSecS2;
	  
	}
      }
    }
  }
  
  //___________________
  
  fTracksRecCuts->Clear();
  fTracksGen->Clear();
  fTracksAODMCCharged->Clear();
  fTracksAODMCChargedSecNS->Clear();
  fTracksAODMCChargedSecS->Clear();
  fTracksRecQualityCuts->Clear();

  fJetsRec->Clear();
  fJetsRecCuts->Clear();
  fJetsGen->Clear();
  fJetsRecEff->Clear();
  fJetsEmbedded->Clear();


  if(fBckgMode && 
     (fBckgType[0]==kBckgClusters || fBckgType[1]==kBckgClusters || fBckgType[2]==kBckgClusters || fBckgType[3]==kBckgClusters || fBckgType[4]==kBckgClusters ||
      fBckgType[0]==kBckgClustersOutLeading || fBckgType[1]==kBckgClustersOutLeading || fBckgType[2]==kBckgClustersOutLeading || 
      fBckgType[3]==kBckgClustersOutLeading || fBckgType[4]==kBckgClustersOutLeading)){
    
    fBckgJetsRec->Clear();
    fBckgJetsRecCuts->Clear();
    fBckgJetsGen->Clear();
  }

  
  //Post output data.
  PostData(1, fCommonHistList);
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

  if(type==kTrackAODExtraCuts || type==kTrackAODExtraonlyCuts || type==kTrackAODExtra || type==kTrackAODExtraonly){
    
    TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(fAOD->FindListObject("aodExtraTracks"));
    if(!aodExtraTracks)return iCount;
    for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
      AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
      if (!track) continue;
      
      AliAODTrack *tr = dynamic_cast<AliAODTrack*> (track);
      if(!tr)continue;

      if(type==kTrackAODExtraCuts || type==kTrackAODExtraonlyCuts){

	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))   continue;
	
	if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax) continue;
	if(tr->Phi() < fTrackPhiMin || tr->Phi() > fTrackPhiMax) continue;
	if(tr->Pt()  < fTrackPtCut) continue;
      }    

      list->Add(tr);
      iCount++;
    }
  }

  if(type==kTrackAODCuts || type==kTrackAODQualityCuts || type==kTrackAOD || type==kTrackAODExtraCuts || type==kTrackAODExtra){

    // all rec. tracks, esd filter mask, eta range
    
    for(Int_t it=0; it<fAOD->GetNumberOfTracks(); ++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      
      if(type == kTrackAODCuts || type==kTrackAODQualityCuts || type==kTrackAODExtraCuts){

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
  else if (type==kTrackAODMCCharged || type==kTrackAODMCAll || type==kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS)  {
    // MC particles (from AOD), physical primaries, all or rather charged or rather charged within acceptance
    if(!fAOD) return -1;
    
    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    
    for(int it=0; it<tca->GetEntriesFast(); ++it){
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
      if(!part)continue;
      if(type != kTrackAODMCChargedSecNS && type != kTrackAODMCChargedSecS  && !part->IsPhysicalPrimary())continue;
      if((type == kTrackAODMCChargedSecNS || type == kTrackAODMCChargedSecS) && part->IsPhysicalPrimary())continue;

      if (type==kTrackAODMCCharged || type==kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS){
	if(part->Charge()==0) continue;

	if(type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS){
	  Bool_t isFromStrange = kFALSE;
	  Int_t iMother = part->GetMother();
	  if(iMother >= 0){
	    AliAODMCParticle *partM = dynamic_cast<AliAODMCParticle*>(tca->At(iMother));
	    if(!partM) continue;

	    Int_t codeM =  TMath::Abs(partM->GetPdgCode());
	    Int_t mfl = Int_t (codeM/ TMath::Power(10, Int_t(TMath::Log10(codeM))));
	    if  (mfl == 3 && codeM != 3) isFromStrange = kTRUE;
	    
	    // if(mfl ==3){
	    //   cout<<" mfl "<<mfl<<" codeM "<<partM->GetPdgCode()<<" code this track "<<part->GetPdgCode()<<endl; 
	    //   cout<<" index this track "<<it<<" index daughter 0 "<<partM->GetDaughter(0)<<" 1 "<<partM->GetDaughter(1)<<endl; 
	    // }

	    if(type==kTrackAODMCChargedSecNS && isFromStrange) continue;
	    if(type==kTrackAODMCChargedSecS  && !isFromStrange) continue;
	  }
	}

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
    if(fBranchRecJets.Length())      aodRecJets = dynamic_cast<TClonesArray*>(fAODJets->FindListObject(fBranchRecJets.Data()));
    if(!aodRecJets)                  aodRecJets = dynamic_cast<TClonesArray*>(fAODJets->GetList()->FindObject(fBranchRecJets.Data()));
    if(fAODExtension&&!aodRecJets)   aodRecJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchRecJets.Data()));

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

 
      list->Add(tmp);
      nRecJets++; 
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
    if(fBranchGenJets.Length()) aodGenJets = dynamic_cast<TClonesArray*>(fAODJets->FindListObject(fBranchGenJets.Data()));
    if(!aodGenJets)             aodGenJets = dynamic_cast<TClonesArray*>(fAODJets->GetList()->FindObject(fBranchGenJets.Data()));
    if(fAODExtension&&!aodGenJets)   aodGenJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchGenJets.Data()));

    if(!aodGenJets){
      if(fDebug>0){
	if(fBranchGenJets.Length()) Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGenJets.Data());
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
  else if(type == kJetsEmbedded){ // embedded jets

    if(fBranchEmbeddedJets.Length()==0){
      Printf("%s:%d no embedded jet branch specified", (char*)__FILE__,__LINE__);
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    TClonesArray *aodEmbeddedJets = 0; 
    if(fBranchEmbeddedJets.Length())      aodEmbeddedJets = dynamic_cast<TClonesArray*>(fAODJets->FindListObject(fBranchEmbeddedJets.Data()));
    if(!aodEmbeddedJets)                  aodEmbeddedJets = dynamic_cast<TClonesArray*>(fAODJets->GetList()->FindObject(fBranchEmbeddedJets.Data()));
    if(fAODExtension&&!aodEmbeddedJets)   aodEmbeddedJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchEmbeddedJets.Data()));

    if(!aodEmbeddedJets){
      if(fBranchEmbeddedJets.Length()) Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fBranchEmbeddedJets.Data());
      if(fDebug>1)fAOD->Print();
      return 0;
    }

    // Reorder jet pt and fill new temporary AliAODJet objects
    Int_t nEmbeddedJets = 0;
    
    for(Int_t ij=0; ij<aodEmbeddedJets->GetEntries(); ++ij){

      AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodEmbeddedJets->At(ij));
      if(!tmp) continue;

      if( tmp->Pt() < fJetPtCut ) continue;
      if(    tmp->Eta() < fJetEtaMin
	  || tmp->Eta() > fJetEtaMax
	  || tmp->Phi() < fJetPhiMin
          || tmp->Phi() > fJetPhiMax ) continue;
      
      list->Add(tmp);
      nEmbeddedJets++;
    }
    
    list->Sort();
    
    return nEmbeddedJets;
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
    if(fBranchRecBckgClusters.Length()) aodRecJets = dynamic_cast<TClonesArray*>(fAODJets->FindListObject(fBranchRecBckgClusters.Data()));
    if(!aodRecJets)                     aodRecJets = dynamic_cast<TClonesArray*>(fAODJets->GetList()->FindObject(fBranchRecBckgClusters.Data()));
    if(fAODExtension&&!aodRecJets)      aodRecJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fBranchRecBckgClusters.Data()));    

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
void AliAnalysisTaskFragmentationFunction::GetJetTracksPointing(TList* inputlist, TList* outputlist, const AliAODJet* jet, 
								   const Double_t radius, Double_t& sumPt, const Double_t minPtL, const Double_t maxPt, Bool_t& isBadPt)
{
  // fill list of tracks in cone around jet axis  

  sumPt = 0;
  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;

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

      if(maxPt>0  && track->Pt()>maxPt)  isBadMaxPt = kTRUE;
      if(minPtL>0 && track->Pt()>minPtL) isBadMinPt = kFALSE;
    }
  }
  
  isBadPt = kFALSE; 
  if(minPtL>0 && isBadMinPt) isBadPt = kTRUE;  
  if(maxPt>0  && isBadMaxPt) isBadPt = kTRUE;  
  
  outputlist->Sort();
}

// _________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetJetTracksTrackrefs(TList* list, const AliAODJet* jet, const Double_t minPtL, const Double_t maxPt, Bool_t& isBadPt)
{
  // list of jet tracks from trackrefs
  
  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;

  for(Int_t itrack=0; itrack<nTracks; itrack++) {
    
    AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(itrack));
    if(!track){
      AliError("expected ref track not found ");
      continue;
    }
    
    if(track->Pt()  < fTrackPtCut) continue; // track refs may contain low pt cut (bug in AliFastJetInput) 
    if(maxPt>0 && track->Pt()>maxPt)   isBadMaxPt = kTRUE;
    if(minPtL>0 && track->Pt()>minPtL) isBadMinPt = kFALSE;

    list->Add(track);
  }
  
  isBadPt = kFALSE; 
  if(minPtL>0 && isBadMinPt) isBadPt = kTRUE;  
  if(maxPt>0 && isBadMaxPt)  isBadPt = kTRUE;  

  list->Sort();
}

// _ ________________________________________________________________________________________________________________________________
void  AliAnalysisTaskFragmentationFunction::AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,
							    TArrayS& isRefGen,TH2F* fh2PtRecVsGen)
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
  isRefGen.Set(nTracksGen);

  indexAODTr.Reset(-1);
  indexMCTr.Reset(-1);
  isRefGen.Reset(0);

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


  // define reference sample of primaries/secondaries (for reconstruction efficiency / contamination)

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksAODMCCharged->At(iGen));
    if(!gentrack)continue;
    Int_t pdg = gentrack->GetPdgCode();    

    // 211 - pi, 2212 - proton, 321 - Kaon, 11 - electron, 13 - muon
    if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || 
       TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){
      
      isRefGen[iGen] = kTRUE;

      Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

      if(iRec>=0){
	Float_t genPt = gentrack->Pt();
	AliAODTrack* vt = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)); 
	if(vt){
	  Float_t recPt = vt->Pt();
	  fh2PtRecVsGen->Fill(genPt,recPt);
	}
      }
    }
  }
}

// _____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::FillSingleTrackHistosRecGen(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, 
								       const TArrayI& indexAODTr, const TArrayS& isRefGen, const Bool_t scaleStrangeness){

  // fill QA for single track reconstruction efficiency
  
  Int_t nTracksGen  = tracksGen->GetSize();

  if(!nTracksGen) return;

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    if(isRefGen[iGen] != 1) continue; // select primaries

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksGen->At(iGen));
    if(!gentrack) continue;
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 

    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtCut) continue;

    if(trackQAGen) trackQAGen->FillTrackQA(etaGen, phiGen, ptGen);

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

    if(iRec>=0 && trackQARec){
      if(scaleStrangeness){ 
	Double_t weight = GetMCStrangenessFactor(ptGen);
	trackQARec->FillTrackQA(etaGen, phiGen, ptGen, kFALSE, 0, kTRUE, weight);
      }
      else trackQARec->FillTrackQA(etaGen, phiGen, ptGen);
    }
  }
}

// ______________________________________________________________________________________________________________________________________________________

void  AliAnalysisTaskFragmentationFunction::FillJetTrackHistosRec(AliFragFuncHistos* ffhistRec, AliAODJet* jet, 
								  TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr,
								  const TArrayS& isRefGen, TList* jetTrackListTR, const Bool_t scaleStrangeness,
								  Bool_t fillJS, TProfile* hProNtracksLeadingJet, TProfile** hProDelRPtSum, TProfile* hProDelR80pcPt)
{
  // fill objects for jet track reconstruction efficiency or secondaries contamination 
  // arguments histGen/histRec can be of different type: AliFragFuncHistos*, AliFragFuncIntraJetHistos*, ...
  // jetTrackListTR pointer: track refs if not NULL  
  

  // ensure proper normalization, even for secondaries
  Double_t jetPtRec = jet->Pt();
  ffhistRec->FillFF(-1, jetPtRec, kTRUE);

  Int_t nTracksJet = jetTrackList->GetSize(); // list with AODMC tracks
  if(nTracksJet == 0) return; 
  
  TList* listRecTracks = new TList(); 
  listRecTracks->Clear();
  
  for(Int_t iTr=0; iTr<nTracksJet; iTr++){ // jet tracks loop
    
    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (jetTrackList->At(iTr));
    if(!gentrack)continue;
    // find jet track in gen tracks list
    Int_t iGen = tracksGen->IndexOf(gentrack); 
    
    if(iGen<0){
      if(fDebug>0) Printf("%s:%d gen jet track not found ",(char*)__FILE__,__LINE__);
      continue;
    }
    
    if(isRefGen[iGen] != 1) continue; // select primaries
    
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // gen level acc & pt cuts - skip in case of track refs  
    if(!jetTrackListTR && (etaGen < fTrackEtaMin || etaGen > fTrackEtaMax)) continue;
    if(!jetTrackListTR && (phiGen < fTrackPhiMin || phiGen > fTrackPhiMax)) continue;
    if(!jetTrackListTR &&  ptGen  < fTrackPtCut) continue;
   

    Double_t ptRec = -1;	

    Int_t iRec   = indexAODTr[iGen]; // can be -1 if no good reconstructed track 
    Bool_t isRec = (iRec>=0) ? kTRUE : kFALSE; 

    Bool_t isJetTrack = kFALSE;
    if(!jetTrackListTR) isJetTrack = kTRUE; // skip trackRefs check for tracks in ideal cone 

    if(isRec){
      
      AliAODTrack* rectrack = dynamic_cast<AliAODTrack*> (tracksRec->At(iRec));
      if(!rectrack) continue;

      ptRec = rectrack->Pt();	
      
      if(jetTrackListTR){ 
        Int_t iRecTR = jetTrackListTR->IndexOf(rectrack); 
        if(iRecTR >=0 ) isJetTrack = kTRUE; // rec tracks assigned to jet 
      }
    
      if(isJetTrack){
	
	Double_t trackPt = ptRec;
	Bool_t incrementJetPt = kFALSE; 
	
	if(scaleStrangeness){
	  Double_t weight = GetMCStrangenessFactor(ptGen);	  
	  ffhistRec->FillFF( trackPt, jetPtRec, incrementJetPt, 0, kTRUE, weight );
	}
	else{
	  ffhistRec->FillFF( trackPt, jetPtRec, incrementJetPt );
	}

	listRecTracks->Add(rectrack);
	
      }
    }
  }


  if(fillJS) FillJetShape(jet,listRecTracks,hProNtracksLeadingJet, hProDelRPtSum, hProDelR80pcPt,0,0,scaleStrangeness); 

  delete listRecTracks;

}

// _____________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius,Double_t& sumPt)
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

    // embedded tracks
    if( fUseExtraTracksBgr != 1){
      if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (inputlist->At(itrack))){
	if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
      }
    }
    
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
void AliAnalysisTaskFragmentationFunction::GetTracksTiltedwrpJetAxisWindow(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius,Double_t& sumPt,Double_t &normFactor)
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

      // embedded tracks
      if( fUseExtraTracksBgr != 1){
	if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (inputlist->At(itrack))){
	  if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	  if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	}
      }
      
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
  Float_t rc  = TMath::Abs(GetFFRadius());
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
    
    // embedded tracks
    if( fUseExtraTracksBgr != 1){
      if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (inputlist->At(itrack))){
	if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
      }
    }

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
  Float_t rc  = TMath::Abs(GetFFRadius());
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
    
    // embedded tracks
    if( fUseExtraTracksBgr != 1){
      if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (inputlist->At(itrack))){
	if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
      }
    }

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
	
      // embedded tracks - note: using ref tracks here, fBranchRecBckgClusters has to be consistent
      if( fUseExtraTracksBgr != 1){
	if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (bgrCluster->GetTrack(it))){
	  if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	  if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	}
      }

      AliVParticle*   track = dynamic_cast<AliVParticle*>(bgrCluster->GetTrack(it));
      if(!track) continue;
	
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
  
  normFactor = 0;

  Int_t nBckgClusters = fBckgJetsRec->GetEntries(); // not 'recCuts': use all clusters in full eta range

  if(nBckgClusters<3) return; // need at least 3 clusters (skipping 2 highest)

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
	
    // embedded tracks - note: using ref tracks here, fBranchRecBckgClusters has to be consistent
    if( fUseExtraTracksBgr != 1){
      if(AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*> (medianCluster->GetTrack(it))){
	if(fUseExtraTracksBgr == 0  &&  ((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
	if(fUseExtraTracksBgr == -1 && !((trackAOD->GetFlags() & AliESDtrack::kEmbedded)>0)) continue; 
      }
    }

    AliVParticle* track = dynamic_cast<AliVParticle*>(medianCluster->GetTrack(it));
    if(!track) continue;
	
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
void AliAnalysisTaskFragmentationFunction::FillBckgHistos(Int_t type, TList* inputtracklist, TList* inputjetlist, AliAODJet* jet, 
							     AliFragFuncHistos* ffbckghistocuts, AliFragFuncQATrackHistos* qabckghistocuts, TH1F* fh1Mult){

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

  Int_t nRecJetsCuts = inputjetlist->GetEntries(); 

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
      TList* tracklistoutleading   = new TList();
      Double_t sumPtOutLeading     = 0.; 
      GetTracksOutOfNJets(1,inputtracklist, tracklistoutleading, inputjetlist, sumPtOutLeading);
      if(type==kBckgOutLJ && fh1Mult) fh1Mult->Fill(tracklistoutleading->GetSize());
      
      for(Int_t it=0; it<tracklistoutleading->GetSize(); ++it){

	AliVParticle* trackVP   = (AliVParticle*)(tracklistoutleading->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOutLJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt);
	    
	    // Fill track QA for background
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==1 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
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
          }
        // All cases included
        if(nRecJetsCuts==1 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
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
      if(type==kBckgOutLJStat && fh1Mult) fh1Mult->Fill(tracklistoutleadingStat->GetSize());

      for(Int_t it=0; it<tracklistoutleadingStat->GetSize(); ++it){

	AliVParticle* trackVP   = dynamic_cast<AliVParticle*>(tracklistoutleadingStat->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	// Stat plots
	if(type==kBckgOutLJStat)
	  {
	    if(fFFMode)ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);

	    // Fill track QA for background
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt); // OB added bgr QA
	  }

	// All cases included
	if(nRecJetsCuts==1 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);
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
          }
        // All cases included
        if(nRecJetsCuts==1 && type==kBckgOutLJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorLeading);
          }
      }

      delete tracklistoutleadingStat;
    }

  if(type==kBckgPerp || type==kBckgPerp2 || type==kBckgPerp2Area)
    {
      Double_t sumPtPerp1 = 0.;
      Double_t sumPtPerp2 = 0.;
      TList* tracklistperp  = new TList();
      TList* tracklistperp1 = new TList();
      TList* tracklistperp2 = new TList();

      Double_t norm = 1;
      if(type == kBckgPerp2)     norm = 2; // in FillFF() scaleFac = 1/norm = 0.5 - account for double area; 
      if(type == kBckgPerp2Area) norm = 2*TMath::Pi()*TMath::Abs(GetFFRadius())*TMath::Abs(GetFFRadius()) / jet->EffectiveAreaCharged(); // in FillFF() scaleFac = 1/norm; 

      GetTracksTiltedwrpJetAxis(TMath::Pi()/2., inputtracklist,tracklistperp1,jet,TMath::Abs(GetFFRadius()),sumPtPerp1);
      if(type==kBckgPerp2 || type==kBckgPerp2Area) GetTracksTiltedwrpJetAxis(-1*TMath::Pi()/2., inputtracklist,tracklistperp2,jet,TMath::Abs(GetFFRadius()),sumPtPerp2);

      tracklistperp->AddAll(tracklistperp1);
      tracklistperp->AddAll(tracklistperp2);

      if(tracklistperp->GetSize() != tracklistperp1->GetSize() + tracklistperp2->GetSize()){
	cout<<" ERROR: tracklistperp size "<<tracklistperp->GetSize()<<" perp1 "<<tracklistperp1->GetSize()<<" perp2 "<<tracklistperp2->GetSize()<<endl;
        exit(0); 
      }

      if(fh1Mult) fh1Mult->Fill(tracklistperp->GetSize());
      
      for(Int_t it=0; it<tracklistperp->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistperp->At(it));
	if(!trackVP)continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
       
	if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, norm );

	// Fill track QA for background
	if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	
	delete trackV;
      }
      // Increment jet pt with one entry in case #tracks outside jets = 0
      if(tracklistperp->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
      }


      if(fJSMode){
        // fill for tracklistperp1/2 separately, divide norm by 2
	if(type==kBckgPerp){ 
	  FillJetShape(jet, tracklistperp, fProNtracksLeadingJetBgrPerp2, fProDelRPtSumBgrPerp2, 0, TMath::Pi()/2., 0., kFALSE); 
        }
	if(type==kBckgPerp2){
          FillJetShape(jet, tracklistperp1, fProNtracksLeadingJetBgrPerp2, fProDelRPtSumBgrPerp2, 0,    TMath::Pi()/2., 0., kFALSE);
          FillJetShape(jet, tracklistperp2, fProNtracksLeadingJetBgrPerp2, fProDelRPtSumBgrPerp2, 0, -1*TMath::Pi()/2., 0., kFALSE);
        }
	if(type==kBckgPerp2Area){ // divide norm by 2: listperp1/2 filled separately
	  FillJetShape(jet, tracklistperp1, fProNtracksLeadingJetBgrPerp2, fProDelRPtSumBgrPerp2, 0,    TMath::Pi()/2., 0.5*norm, kFALSE);
	  FillJetShape(jet, tracklistperp2, fProNtracksLeadingJetBgrPerp2, fProDelRPtSumBgrPerp2, 0, -1*TMath::Pi()/2., 0.5*norm, kFALSE);
        }
      }
      
      delete tracklistperp;
      delete tracklistperp1;
      delete tracklistperp2;
    }

 if(type==kBckgASide)
    {
      Double_t sumPtASide = 0.;
      TList* tracklistaside = new TList();
      GetTracksTiltedwrpJetAxis(TMath::Pi(),inputtracklist,tracklistaside,jet,TMath::Abs(GetFFRadius()),sumPtASide);
      if(fh1Mult) fh1Mult->Fill(tracklistaside->GetSize());

      for(Int_t it=0; it<tracklistaside->GetSize(); ++it){
	
        AliVParticle*   trackVP = (AliVParticle*)(tracklistaside->At(it));
	if(!trackVP) continue;
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();

        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);

        delete trackV;
      }
      if(tracklistaside->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
      }

      delete tracklistaside;
    }

  if(type==kBckgASideWindow)
    {
      Double_t normFactorASide = 0.;
      Double_t sumPtASideW = 0.;
      TList* tracklistasidew = new TList();
      GetTracksTiltedwrpJetAxisWindow(TMath::Pi(),inputtracklist,tracklistasidew,jet,TMath::Abs(GetFFRadius()),sumPtASideW,normFactorASide);
      if(fh1Mult) fh1Mult->Fill(tracklistasidew->GetSize());

      for(Int_t it=0; it<tracklistasidew->GetSize(); ++it){

        AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistasidew->At(it));
	if(!trackVP) continue;
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();
        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorASide);

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt, kFALSE, normFactorASide);

        delete trackV;
      }
      if(tracklistasidew->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorASide);
      }

      delete tracklistasidew;
    }

  if(type==kBckgPerpWindow)
    {
      Double_t normFactorPerp = 0.;
      Double_t sumPtPerpW = 0.;
      TList* tracklistperpw = new TList();
      GetTracksTiltedwrpJetAxisWindow(TMath::Pi()/2.,inputtracklist,tracklistperpw,jet,TMath::Abs(GetFFRadius()),sumPtPerpW,normFactorPerp);
      if(fh1Mult) fh1Mult->Fill(tracklistperpw->GetSize());

      for(Int_t it=0; it<tracklistperpw->GetSize(); ++it){
	
        AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistperpw->At(it));
	if(!trackVP) continue;
        TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

        Float_t jetPt   = jet->Pt();
        Float_t trackPt = trackV->Pt();
        Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

        if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorPerp);

        // Fill track QA for background
        if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt, kFALSE, normFactorPerp);

        delete trackV;
      }
      if(tracklistperpw->GetSize()==0) {
         Float_t jetPt = jet->Pt();
         Bool_t incrementJetPt = kTRUE;
         if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactorPerp);
      }

      delete tracklistperpw;
    }


  if(type==kBckgOut2J || type==kBckgOutAJ)
    {
      if(type==kBckgOut2J && fh1Mult) fh1Mult->Fill(tracklistout2jets->GetSize());
      for(Int_t it=0; it<tracklistout2jets->GetSize(); ++it){

	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jets->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOut2J)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );	    
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==2 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
	    
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
          }
        // All cases included
        if(nRecJetsCuts==2 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
          }
      }
    }
  
  if(type==kBckgOut2JStat || type==kBckgOutAJStat)
    {
      for(Int_t it=0; it<tracklistout2jetsStat->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jetsStat->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(type==kBckgOut2JStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	    
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt ); // OB added bgr QA
	  }

	// All cases included
	if(nRecJetsCuts==2 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	     
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
          }
        // All cases included
        if(nRecJetsCuts==2 && type==kBckgOutAJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor2Jets);
          }
      }
      
    }

  if(type==kBckgOut3J || type==kBckgOutAJ)
    {
      if(type==kBckgOut3J && fh1Mult) fh1Mult->Fill(tracklistout3jets->GetSize());
      
      for(Int_t it=0; it<tracklistout3jets->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jets->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(type==kBckgOut3J)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
    
	    qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==3 && type==kBckgOutAJ)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt );
    
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
          }
        // All cases included
        if(nRecJetsCuts==3 && type==kBckgOutAJ)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt );
          }
      }
    }

  if(type==kBckgOut3JStat || type==kBckgOutAJStat)
    {
      for(Int_t it=0; it<tracklistout3jetsStat->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jetsStat->At(it));
	if(!trackVP) continue;
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	
	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(type==kBckgOut3JStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor3Jets);
	    	    
	    //if(fQAMode&1) 	qabckghistocuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	  }

	// All cases included
	if(nRecJetsCuts==3 && type==kBckgOutAJStat)
	  {
	    if(fFFMode) ffbckghistocuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor3Jets );
	    
	    if(fQAMode&1) qabckghistocuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt );

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
          }
        // All cases included
        if(nRecJetsCuts==3 && type==kBckgOutAJStat)
          {
            if(fFFMode) ffbckghistocuts->FillFF( -1., jetPt, incrementJetPt, normFactor3Jets);
          }
      }

    }

  if(type==kBckgClustersOutLeading){ // clusters bgr: all tracks in clusters out of leading jet
    
    TList* tracklistClustersOutLeading = new TList();
    Double_t normFactorClusters = 0;
    Float_t jetPt   = jet->Pt();
    
    GetClusterTracksOutOf1Jet(jet, tracklistClustersOutLeading, normFactorClusters);
    if(fh1Mult) fh1Mult->Fill(tracklistClustersOutLeading->GetSize());
    
    for(Int_t it=0; it<tracklistClustersOutLeading->GetSize(); ++it){
      
      AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistClustersOutLeading->At(it));
      if(!trackVP) continue;
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
    if(fh1Mult) fh1Mult->Fill(tracklistClustersMedian->GetSize());

    for(Int_t it=0; it<tracklistClustersMedian->GetSize(); ++it){
      
      AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistClustersMedian->At(it));
      if(!trackVP) continue;
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

// -----------------------------------------------------------------

Double_t AliAnalysisTaskFragmentationFunction::GetMCStrangenessFactor(const Double_t pt)
{
  // factor strangeness data/MC as function of pt from UE analysis (Sara Vallero)

  Double_t alpha = 1;

  if(0.150<pt && pt<0.200) alpha = 3.639;
  if(0.200<pt && pt<0.250) alpha = 2.097;
  if(0.250<pt && pt<0.300) alpha = 1.930;
  if(0.300<pt && pt<0.350) alpha = 1.932;
  if(0.350<pt && pt<0.400) alpha = 1.943;
  if(0.400<pt && pt<0.450) alpha = 1.993;
  if(0.450<pt && pt<0.500) alpha = 1.989;
  if(0.500<pt && pt<0.600) alpha = 1.963;
  if(0.600<pt && pt<0.700) alpha = 1.917;
  if(0.700<pt && pt<0.800) alpha = 1.861;
  if(0.800<pt && pt<0.900) alpha = 1.820;
  if(0.900<pt && pt<1.000) alpha = 1.741;
  if(1.000<pt && pt<1.500) alpha = 0.878;

  return alpha;
}

// ---------------------------------------------------------------------------------------------------------------------------------
void  AliAnalysisTaskFragmentationFunction::FillJetShape(AliAODJet* jet, TList* list,  
							 TProfile* hProNtracksLeadingJet, TProfile** hProDelRPtSum, TProfile* hProDelR80pcPt, 
							 Double_t dPhiUE, Double_t normUE, Bool_t scaleStrangeness){
  
  const Int_t   kNbinsR    = 50; 
  const Float_t kBinWidthR = 0.02;
  
  Int_t nJetTracks = list->GetEntries();
  
  Float_t PtSumA[kNbinsR]     = {0.0};
  Float_t PtWeightsA[kNbinsR] = {0.0};
  Float_t nTracksA[kNbinsR]   = {0.0};
  
  Float_t *delRA     = new Float_t[nJetTracks];
  Float_t *trackPtA  = new Float_t[nJetTracks];
  Int_t   *index     = new Int_t[nJetTracks];
  
  for(Int_t i=0; i<nJetTracks; i++){
    delRA[i]    = 0;
    trackPtA[i] = 0;
    index[i]    = 0;
  }
  
  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);
  
  if(TMath::Abs(dPhiUE)>0){
    Double_t phiTilted = jet3mom.Phi();
    phiTilted += dPhiUE;
    phiTilted = TVector2::Phi_0_2pi(phiTilted);
    jet3mom.SetPhi(phiTilted);
  }
  
  Double_t jetPt = jet->Pt();
  Double_t sumWeights = 0;
  
  for (Int_t j =0; j<nJetTracks; j++){
  
    AliVParticle* track = dynamic_cast<AliVParticle*>(list->At(j));
    if(!track)continue;
    
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);
    
    Double_t dR = jet3mom.DeltaR(track3mom);

    delRA[j]    = dR;
    trackPtA[j] = track->Pt();
    
    Double_t weight = GetMCStrangenessFactor(track->Pt()); // more correctly should be gen pt
    sumWeights += weight;

    for(Int_t ibin=1; ibin<=kNbinsR; ibin++){
      Float_t xlow = kBinWidthR*(ibin-1);
      Float_t xup  = kBinWidthR*ibin;
      if(xlow <= dR && dR < xup){
	PtSumA[ibin-1]     += track->Pt();
	PtWeightsA[ibin-1] += weight; 
	nTracksA[ibin-1]   += 1; 
      }
    }
  } // track loop
  
  Float_t jetPtMin=0;
  Float_t jetPtMax=0;
  
  for(Int_t ibin=0; ibin<kNbinsR; ibin++){
    Float_t fR =  kBinWidthR*(ibin+0.5);
    
    for(Int_t k=0; k<5; k++){
      if(k==0){jetPtMin=20.0;jetPtMax=30.0;}
      if(k==1){jetPtMin=30.0;jetPtMax=40.0;}
      if(k==2){jetPtMin=40.0;jetPtMax=60.0;}
      if(k==3){jetPtMin=60.0;jetPtMax=80.0;}
      if(k==4){jetPtMin=80.0;jetPtMax=100.0;}
      if(jetPt>jetPtMin && jetPt<jetPtMax){
	
	if(scaleStrangeness){
	  if(nTracksA[ibin]) hProDelRPtSum[k]->Fill(fR,PtSumA[ibin],PtWeightsA[ibin]/nTracksA[ibin]);
	  else               hProDelRPtSum[k]->Fill(fR,PtSumA[ibin],0);
	}
	else                 hProDelRPtSum[k]->Fill(fR,PtSumA[ibin]);
      }
    }
  }
  
  if(scaleStrangeness) hProNtracksLeadingJet->Fill(jetPt,sumWeights);
  else                 hProNtracksLeadingJet->Fill(jetPt,nJetTracks);
  
  if(normUE)           hProNtracksLeadingJet->Fill(jetPt,nJetTracks/normUE);
  
  if(hProDelR80pcPt){
    
    Float_t PtSum = 0;
    Float_t delRPtSum80pc = 0;
    
    TMath::Sort(nJetTracks,delRA,index,0);
    
    for(Int_t ii=0; ii<nJetTracks; ii++){
      
      if(scaleStrangeness){ 
        Double_t weight = GetMCStrangenessFactor(trackPtA[index[ii]]); // more correctly should be gen pt
	PtSum += weight*trackPtA[index[ii]];  
      }
      else PtSum += trackPtA[index[ii]];
      

      if(PtSum/jetPt >= 0.8000){
	delRPtSum80pc = delRA[index[ii]];
	break;
      }
    } 
    hProDelR80pcPt->Fill(jetPt,delRPtSum80pc);
  }
  
  delete[] delRA;
  delete[] trackPtA;
  delete[] index;
}
