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
#include "AliESDInputHandler.h"
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
#include "AliPIDResponse.h"
#include "AliIDFFUtils.h"

#include "AliAnalysisTaskIDFFTCF.h"
using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliAnalysisTaskIDFFTCF)

Bool_t AliAnalysisTaskIDFFTCF::fkDump = kFALSE;

//____________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliAnalysisTaskIDFFTCF()
   : AliAnalysisTaskSE()
   ,fESD(0)
   ,fAOD(0)
   ,fAODJets(0)  
   ,fAODExtension(0)
   ,fNonStdFile("")
   ,fBranchRecJets("jets")
   ,fBranchGenJets("")
   ,fTrackTypeGen(0)
   ,fJetTypeGen(0)
   ,fJetTypeRecEff(0)
   ,fUseAODInputJets(kTRUE)
   ,fFilterMask(0)
   ,fUsePhysicsSelection(kTRUE)
   ,fEvtSelectionMask(0)
   ,fEventClass(0)
   ,fMaxVertexZ(10)
   ,fLeadingJets(kFALSE)
   ,fTPCCutMode(kPIDNone) 
   ,fTOFCutMode(1)
   ,fStream(0x0)
   ,fTree(0x0)
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
   ,fFFRadius(0)
   ,fFFMinLTrackPt(-1)
   ,fFFMaxTrackPt(-1)
   ,fFFMinnTracks(0)
   ,fQAMode(0)
   ,fFFMode(0)
   ,fEffMode(0)
   ,fAvgTrials(0)
   ,fTracksRecCuts(0)
   ,fTracksGen(0)
   ,fTracksAODMCCharged(0)
   ,fTracksAODMCChargedSec(0)
   ,fTracksRecQualityCuts(0)
   ,fJetsRec(0)
   ,fJetsRecCuts(0)
   ,fJetsGen(0)
   ,fJetsRecEff(0)
   ,fQATrackHistosRecCuts(0)
   ,fQATrackHistosGen(0)
   ,fQAJetHistosRec(0)
   ,fQAJetHistosRecCuts(0)
   ,fQAJetHistosRecCutsLeading(0)
   ,fQAJetHistosGen(0)
   ,fQAJetHistosGenLeading(0)
   ,fQAJetHistosRecEffLeading(0)
   ,fFFHistosRecCutsInc(0)
   ,fFFHistosRecCutsIncPi(0)   
   ,fFFHistosRecCutsIncPro(0)  
   ,fFFHistosRecCutsIncK(0)    
   ,fFFHistosRecCutsIncEl(0)   
   ,fFFHistosRecCutsIncMu(0)   
   ,fFFHistosRecLeadingTrack(0)
   ,fFFHistosGenInc(0)
   ,fFFHistosGenIncPi(0)   
   ,fFFHistosGenIncPro(0)  
   ,fFFHistosGenIncK(0)    
   ,fFFHistosGenIncEl(0)   
   ,fFFHistosGenIncMu(0)   
   ,fFFHistosGenLeadingTrack(0)
   ,fQATrackHighPtThreshold(0)
   ,fTHnIDFF(0x0)
   ,fTHnIncl(0x0)
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
   ,fh2PtRecVsGenPrim(0)
   ,fh2PtRecVsGenSec(0)
   ,fQATrackHistosRecEffGen(0)  
   ,fQATrackHistosRecEffRec(0)
   ,fQATrackHistosSecRec(0)   
   ,fQATrackHistosRecEffGenPi(0)  
   ,fQATrackHistosRecEffGenPro(0) 
   ,fQATrackHistosRecEffGenK(0)   
   ,fQATrackHistosRecEffGenEl(0)  
   ,fQATrackHistosRecEffGenMu(0)  
   ,fQATrackHistosRecEffRecPi(0)  
   ,fQATrackHistosRecEffRecPro(0) 
   ,fQATrackHistosRecEffRecK(0)   
   ,fQATrackHistosRecEffRecEl(0)  
   ,fQATrackHistosRecEffRecMu(0)
   ,fQATrackHistosRecEffRecProGFL(0)
   ,fQATrackHistosRecEffRecKGFL(0)
   ,fQATrackHistosSecRecPi(0)  
   ,fQATrackHistosSecRecPro(0) 
   ,fQATrackHistosSecRecK(0)   
   ,fQATrackHistosSecRecEl(0)  
   ,fQATrackHistosSecRecMu(0)  
   ,fQATrackHistosSecRecProGFL(0) 
   ,fQATrackHistosSecRecKGFL(0)
   ,fFFHistosRecEffRec(0)
   ,fFFHistosSecRec(0)
   ,fFFHistosRecEffRecPi(0) 
   ,fFFHistosRecEffRecPro(0)
   ,fFFHistosRecEffRecK(0)  
   ,fFFHistosRecEffRecEl(0) 
   ,fFFHistosRecEffRecMu(0)
   ,fFFHistosRecEffRecProGFL(0)
   ,fFFHistosRecEffRecKGFL(0)
   ,fFFHistosSecRecPi(0)    
   ,fFFHistosSecRecPro(0)   
   ,fFFHistosSecRecK(0)     
   ,fFFHistosSecRecEl(0)    
   ,fFFHistosSecRecMu(0) 
   ,fFFHistosSecRecProGFL(0)   
   ,fFFHistosSecRecKGFL(0)
   ,fRandom(0)
{
   // default constructor
}

//_______________________________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliAnalysisTaskIDFFTCF(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fAODJets(0)  
  ,fAODExtension(0)
  ,fNonStdFile("")
  ,fBranchRecJets("jets")
  ,fBranchGenJets("")
  ,fTrackTypeGen(0)
  ,fJetTypeGen(0)
  ,fJetTypeRecEff(0)
  ,fUseAODInputJets(kTRUE)
  ,fFilterMask(0)
  ,fUsePhysicsSelection(kTRUE)
  ,fEvtSelectionMask(0)
  ,fEventClass(0)
  ,fMaxVertexZ(10)
  ,fLeadingJets(kFALSE)
  ,fTPCCutMode(kPIDNone)  
  ,fTOFCutMode(1)
  ,fStream(0x0)
  ,fTree(0x0)
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
  ,fFFRadius(0)
  ,fFFMinLTrackPt(-1)
  ,fFFMaxTrackPt(-1)
  ,fFFMinnTracks(0)
  ,fQAMode(0)
  ,fFFMode(0)
  ,fEffMode(0)
  ,fAvgTrials(0)
  ,fTracksRecCuts(0)
  ,fTracksGen(0)
  ,fTracksAODMCCharged(0)
  ,fTracksAODMCChargedSec(0)
  ,fTracksRecQualityCuts(0)
  ,fJetsRec(0)
  ,fJetsRecCuts(0)
  ,fJetsGen(0)
  ,fJetsRecEff(0)
  ,fQATrackHistosRecCuts(0)
  ,fQATrackHistosGen(0)
  ,fQAJetHistosRec(0)
  ,fQAJetHistosRecCuts(0)
  ,fQAJetHistosRecCutsLeading(0)
  ,fQAJetHistosGen(0)
  ,fQAJetHistosGenLeading(0)
  ,fQAJetHistosRecEffLeading(0)
  ,fFFHistosRecCutsInc(0)
  ,fFFHistosRecCutsIncPi(0)   
  ,fFFHistosRecCutsIncPro(0)  
  ,fFFHistosRecCutsIncK(0)    
  ,fFFHistosRecCutsIncEl(0)   
  ,fFFHistosRecCutsIncMu(0)   
  ,fFFHistosRecLeadingTrack(0)
  ,fFFHistosGenInc(0)
  ,fFFHistosGenIncPi(0)   
  ,fFFHistosGenIncPro(0)  
  ,fFFHistosGenIncK(0)    
  ,fFFHistosGenIncEl(0)   
  ,fFFHistosGenIncMu(0)   
  ,fFFHistosGenLeadingTrack(0)
  ,fQATrackHighPtThreshold(0) 
  ,fTHnIDFF(0x0)
  ,fTHnIncl(0x0)
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
  ,fh2PtRecVsGenPrim(0)
  ,fh2PtRecVsGenSec(0)
  ,fQATrackHistosRecEffGen(0)  
  ,fQATrackHistosRecEffRec(0)
  ,fQATrackHistosSecRec(0) 
  ,fQATrackHistosRecEffGenPi(0)  
  ,fQATrackHistosRecEffGenPro(0) 
  ,fQATrackHistosRecEffGenK(0)   
  ,fQATrackHistosRecEffGenEl(0)  
  ,fQATrackHistosRecEffGenMu(0)  
  ,fQATrackHistosRecEffRecPi(0)  
  ,fQATrackHistosRecEffRecPro(0) 
  ,fQATrackHistosRecEffRecK(0)   
  ,fQATrackHistosRecEffRecEl(0)  
  ,fQATrackHistosRecEffRecMu(0)
  ,fQATrackHistosRecEffRecProGFL(0)
  ,fQATrackHistosRecEffRecKGFL(0)
  ,fQATrackHistosSecRecPi(0)  
  ,fQATrackHistosSecRecPro(0) 
  ,fQATrackHistosSecRecK(0)   
  ,fQATrackHistosSecRecEl(0)  
  ,fQATrackHistosSecRecMu(0)
  ,fQATrackHistosSecRecProGFL(0) 
  ,fQATrackHistosSecRecKGFL(0)
  ,fFFHistosRecEffRec(0)
  ,fFFHistosSecRec(0)
  ,fFFHistosRecEffRecPi(0) 
  ,fFFHistosRecEffRecPro(0)
  ,fFFHistosRecEffRecK(0)  
  ,fFFHistosRecEffRecEl(0) 
  ,fFFHistosRecEffRecMu(0)
  ,fFFHistosRecEffRecProGFL(0)
  ,fFFHistosRecEffRecKGFL(0)
  ,fFFHistosSecRecPi(0)    
  ,fFFHistosSecRecPro(0)   
  ,fFFHistosSecRecK(0)     
  ,fFFHistosSecRecEl(0)    
  ,fFFHistosSecRecMu(0)
  ,fFFHistosSecRecProGFL(0)   
  ,fFFHistosSecRecKGFL(0)  
  ,fRandom(0)
{
  // constructor
   
  DefineOutput(1,TList::Class());

  if(fkDump){
    fStream = new TTreeStream("tree");
    DefineOutput(2, TTree::Class());
    fTree= fStream->GetTree();
  }
}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliAnalysisTaskIDFFTCF(const  AliAnalysisTaskIDFFTCF &copy)
  : AliAnalysisTaskSE()
  ,fESD(copy.fESD)
  ,fAOD(copy.fAOD)
  ,fAODJets(copy.fAODJets)  
  ,fAODExtension(copy.fAODExtension)
  ,fNonStdFile(copy.fNonStdFile)
  ,fBranchRecJets(copy.fBranchRecJets)
  ,fBranchGenJets(copy.fBranchGenJets)
  ,fTrackTypeGen(copy.fTrackTypeGen)
  ,fJetTypeGen(copy.fJetTypeGen)
  ,fJetTypeRecEff(copy.fJetTypeRecEff)
  ,fUseAODInputJets(copy.fUseAODInputJets)
  ,fFilterMask(copy.fFilterMask)
  ,fUsePhysicsSelection(copy.fUsePhysicsSelection)
  ,fEvtSelectionMask(copy.fEvtSelectionMask)
  ,fEventClass(copy.fEventClass)
  ,fMaxVertexZ(copy.fMaxVertexZ)
  ,fLeadingJets(copy.fLeadingJets)
  ,fTPCCutMode(copy.fTPCCutMode)
  ,fTOFCutMode(copy.fTOFCutMode)
  ,fStream(copy.fStream)
  ,fTree(copy.fTree)
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
  ,fFFRadius(copy.fFFRadius)
  ,fFFMinLTrackPt(copy.fFFMinLTrackPt)
  ,fFFMaxTrackPt(copy.fFFMaxTrackPt)
  ,fFFMinnTracks(copy.fFFMinnTracks)
  ,fQAMode(copy.fQAMode)
  ,fFFMode(copy.fFFMode)
  ,fEffMode(copy.fEffMode)
  ,fAvgTrials(copy.fAvgTrials)
  ,fTracksRecCuts(copy.fTracksRecCuts)
  ,fTracksGen(copy.fTracksGen)
  ,fTracksAODMCCharged(copy.fTracksAODMCCharged)
  ,fTracksAODMCChargedSec(copy.fTracksAODMCChargedSec)
  ,fTracksRecQualityCuts(copy.fTracksRecQualityCuts)
  ,fJetsRec(copy.fJetsRec)
  ,fJetsRecCuts(copy.fJetsRecCuts)
  ,fJetsGen(copy.fJetsGen)
  ,fJetsRecEff(copy.fJetsRecEff)
  ,fQATrackHistosRecCuts(copy.fQATrackHistosRecCuts)
  ,fQATrackHistosGen(copy.fQATrackHistosGen)
  ,fQAJetHistosRec(copy.fQAJetHistosRec)
  ,fQAJetHistosRecCuts(copy.fQAJetHistosRecCuts)
  ,fQAJetHistosRecCutsLeading(copy.fQAJetHistosRecCutsLeading)
  ,fQAJetHistosGen(copy.fQAJetHistosGen)
  ,fQAJetHistosGenLeading(copy.fQAJetHistosGenLeading)
  ,fQAJetHistosRecEffLeading(copy.fQAJetHistosRecEffLeading)
  ,fFFHistosRecCutsInc(copy.fFFHistosRecCutsInc)
  ,fFFHistosRecCutsIncPi(copy.fFFHistosRecCutsIncPi)
  ,fFFHistosRecCutsIncPro(copy.fFFHistosRecCutsIncPro) 
  ,fFFHistosRecCutsIncK(copy.fFFHistosRecCutsIncK)
  ,fFFHistosRecCutsIncEl(copy.fFFHistosRecCutsIncEl) 
  ,fFFHistosRecCutsIncMu(copy.fFFHistosRecCutsIncMu) 
  ,fFFHistosRecLeadingTrack(copy.fFFHistosRecLeadingTrack)
  ,fFFHistosGenInc(copy.fFFHistosGenInc)
  ,fFFHistosGenIncPi(copy.fFFHistosGenIncPi)
  ,fFFHistosGenIncPro(copy.fFFHistosGenIncPro) 
  ,fFFHistosGenIncK(copy.fFFHistosGenIncK)
  ,fFFHistosGenIncEl(copy.fFFHistosGenIncEl) 
  ,fFFHistosGenIncMu(copy.fFFHistosGenIncMu) 
  ,fFFHistosGenLeadingTrack(copy.fFFHistosGenLeadingTrack)
  ,fQATrackHighPtThreshold(copy.fQATrackHighPtThreshold) 
  ,fTHnIDFF(copy.fTHnIDFF)
  ,fTHnIncl(copy.fTHnIncl)
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
  ,fh2PtRecVsGenPrim(copy.fh2PtRecVsGenPrim)
  ,fh2PtRecVsGenSec(copy.fh2PtRecVsGenSec)
  ,fQATrackHistosRecEffGen(copy.fQATrackHistosRecEffGen)  
  ,fQATrackHistosRecEffRec(copy.fQATrackHistosRecEffRec)  
  ,fQATrackHistosSecRec(copy.fQATrackHistosSecRec)  
  ,fQATrackHistosRecEffGenPi(copy.fQATrackHistosRecEffGenPi)  
  ,fQATrackHistosRecEffGenPro(copy.fQATrackHistosRecEffGenPro) 
  ,fQATrackHistosRecEffGenK(copy.fQATrackHistosRecEffGenK)
  ,fQATrackHistosRecEffGenEl(copy.fQATrackHistosRecEffGenEl) 
  ,fQATrackHistosRecEffGenMu(copy.fQATrackHistosRecEffGenMu) 
  ,fQATrackHistosRecEffRecPi(copy.fQATrackHistosRecEffRecPi) 
  ,fQATrackHistosRecEffRecPro(copy.fQATrackHistosRecEffRecPro) 
  ,fQATrackHistosRecEffRecK(copy.fQATrackHistosRecEffRecK)
  ,fQATrackHistosRecEffRecEl(copy.fQATrackHistosRecEffRecEl) 
  ,fQATrackHistosRecEffRecMu(copy.fQATrackHistosRecEffRecMu)   
  ,fQATrackHistosRecEffRecProGFL(copy.fQATrackHistosRecEffRecProGFL)
  ,fQATrackHistosRecEffRecKGFL(copy.fQATrackHistosRecEffRecKGFL)  
  ,fQATrackHistosSecRecPi(copy.fQATrackHistosSecRecPi) 
  ,fQATrackHistosSecRecPro(copy.fQATrackHistosSecRecPro) 
  ,fQATrackHistosSecRecK(copy.fQATrackHistosSecRecK)
  ,fQATrackHistosSecRecEl(copy.fQATrackHistosSecRecEl) 
  ,fQATrackHistosSecRecMu(copy.fQATrackHistosSecRecMu) 
  ,fQATrackHistosSecRecProGFL(copy.fQATrackHistosSecRecProGFL) 
  ,fQATrackHistosSecRecKGFL(copy.fQATrackHistosSecRecKGFL)
  ,fFFHistosRecEffRec(copy.fFFHistosRecEffRec)  
  ,fFFHistosSecRec(copy.fFFHistosSecRec)   
  ,fFFHistosRecEffRecPi(copy.fFFHistosRecEffRecPi)  
  ,fFFHistosRecEffRecPro(copy.fFFHistosRecEffRecPro)
  ,fFFHistosRecEffRecK(copy.fFFHistosRecEffRecK)  
  ,fFFHistosRecEffRecEl(copy.fFFHistosRecEffRecEl) 
  ,fFFHistosRecEffRecMu(copy.fFFHistosRecEffRecMu) 
  ,fFFHistosRecEffRecProGFL(copy.fFFHistosRecEffRecProGFL)
  ,fFFHistosRecEffRecKGFL(copy.fFFHistosRecEffRecKGFL)  
  ,fFFHistosSecRecPi(copy.fFFHistosSecRecPi)    
  ,fFFHistosSecRecPro(copy.fFFHistosSecRecPro)   
  ,fFFHistosSecRecK(copy.fFFHistosSecRecK)     
  ,fFFHistosSecRecEl(copy.fFFHistosSecRecEl)    
  ,fFFHistosSecRecMu(copy.fFFHistosSecRecMu)
  ,fFFHistosSecRecProGFL(copy.fFFHistosSecRecProGFL)   
  ,fFFHistosSecRecKGFL(copy.fFFHistosSecRecKGFL)  
  ,fRandom(copy.fRandom)
{
  // copy constructor
}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskIDFFTCF& AliAnalysisTaskIDFFTCF::operator=(const AliAnalysisTaskIDFFTCF& o)
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
    fBranchGenJets                 = o.fBranchGenJets;
    fTrackTypeGen                  = o.fTrackTypeGen;
    fJetTypeGen                    = o.fJetTypeGen;
    fJetTypeRecEff                 = o.fJetTypeRecEff;
    fUseAODInputJets               = o.fUseAODInputJets;
    fFilterMask                    = o.fFilterMask;
    fUsePhysicsSelection           = o.fUsePhysicsSelection;
    fEvtSelectionMask              = o.fEvtSelectionMask;
    fEventClass                    = o.fEventClass;
    fMaxVertexZ                    = o.fMaxVertexZ;
    fLeadingJets                   = o.fLeadingJets;
    fTPCCutMode                    = o.fTPCCutMode;
    fTOFCutMode                    = o.fTOFCutMode;
    fStream                        = o.fStream;
    fTree                          = o.fTree;
    fTrackPtCut                    = o.fTrackPtCut;
    fTrackEtaMin                   = o.fTrackEtaMin;
    fTrackEtaMax                   = o.fTrackEtaMax;
    fTrackPhiMin                   = o.fTrackPhiMin;
    fTrackPhiMax                   = o.fTrackPhiMax;
    fJetPtCut                      = o.fJetPtCut;
    fJetEtaMin                     = o.fJetEtaMin;
    fJetEtaMax                     = o.fJetEtaMax;
    fJetPhiMin                     = o.fJetPhiMin;
    fJetPhiMax                     = o.fJetPhiMin;
    fFFRadius                      = o.fFFRadius;
    fFFMinLTrackPt                 = o.fFFMinLTrackPt;
    fFFMaxTrackPt                  = o.fFFMaxTrackPt;
    fFFMinnTracks                  = o.fFFMinnTracks;
    fQAMode                        = o.fQAMode;
    fFFMode                        = o.fFFMode;
    fEffMode                       = o.fEffMode;
    fAvgTrials                     = o.fAvgTrials;
    fTracksRecCuts                 = o.fTracksRecCuts;
    fTracksGen                     = o.fTracksGen;
    fTracksAODMCCharged            = o.fTracksAODMCCharged;
    fTracksAODMCChargedSec         = o.fTracksAODMCChargedSec;
    fTracksRecQualityCuts          = o.fTracksRecQualityCuts;
    fJetsRec                       = o.fJetsRec;
    fJetsRecCuts                   = o.fJetsRecCuts;
    fJetsGen                       = o.fJetsGen;
    fJetsRecEff                    = o.fJetsRecEff;
    fQATrackHistosRecCuts          = o.fQATrackHistosRecCuts;
    fQATrackHistosGen              = o.fQATrackHistosGen;
    fQAJetHistosRec                = o.fQAJetHistosRec;
    fQAJetHistosRecCuts            = o.fQAJetHistosRecCuts;
    fQAJetHistosRecCutsLeading     = o.fQAJetHistosRecCutsLeading;
    fQAJetHistosGen                = o.fQAJetHistosGen;
    fQAJetHistosGenLeading         = o.fQAJetHistosGenLeading;
    fQAJetHistosRecEffLeading      = o.fQAJetHistosRecEffLeading;
    fFFHistosRecCutsInc            = o.fFFHistosRecCutsInc;
    fFFHistosRecCutsIncPi          = o.fFFHistosRecCutsIncPi;
    fFFHistosRecCutsIncPro         = o.fFFHistosRecCutsIncPro; 
    fFFHistosRecCutsIncK           = o.fFFHistosRecCutsIncK;
    fFFHistosRecCutsIncEl          = o.fFFHistosRecCutsIncEl; 
    fFFHistosRecCutsIncMu          = o.fFFHistosRecCutsIncMu; 
    fFFHistosRecLeadingTrack       = o.fFFHistosRecLeadingTrack;
    fFFHistosGenInc                = o.fFFHistosGenInc;
    fFFHistosGenIncPi              = o.fFFHistosGenIncPi;
    fFFHistosGenIncPro             = o.fFFHistosGenIncPro; 
    fFFHistosGenIncK               = o.fFFHistosGenIncK;
    fFFHistosGenIncEl              = o.fFFHistosGenIncEl; 
    fFFHistosGenIncMu              = o.fFFHistosGenIncMu; 
    fFFHistosGenLeadingTrack       = o.fFFHistosGenLeadingTrack;
    fQATrackHighPtThreshold        = o.fQATrackHighPtThreshold; 
    fTHnIDFF                       = o.fTHnIDFF;
    fTHnIncl                       = o.fTHnIncl;
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
    fh2PtRecVsGenPrim              = o.fh2PtRecVsGenPrim;
    fh2PtRecVsGenSec               = o.fh2PtRecVsGenSec; 
    fQATrackHistosRecEffGen        = o.fQATrackHistosRecEffGen;  
    fQATrackHistosRecEffRec        = o.fQATrackHistosRecEffRec;  
    fQATrackHistosSecRec           = o.fQATrackHistosSecRec;  
    fQATrackHistosRecEffGenPi      = o.fQATrackHistosRecEffGenPi;  
    fQATrackHistosRecEffGenPro     = o.fQATrackHistosRecEffGenPro; 
    fQATrackHistosRecEffGenK       = o.fQATrackHistosRecEffGenK;
    fQATrackHistosRecEffGenEl      = o.fQATrackHistosRecEffGenEl; 
    fQATrackHistosRecEffGenMu      = o.fQATrackHistosRecEffGenMu; 
    fQATrackHistosRecEffRecPi      = o.fQATrackHistosRecEffRecPi; 
    fQATrackHistosRecEffRecPro     = o.fQATrackHistosRecEffRecPro; 
    fQATrackHistosRecEffRecK       = o.fQATrackHistosRecEffRecK;
    fQATrackHistosRecEffRecEl      = o.fQATrackHistosRecEffRecEl; 
    fQATrackHistosRecEffRecMu      = o.fQATrackHistosRecEffRecMu; 
    fQATrackHistosRecEffRecProGFL  = o.fQATrackHistosRecEffRecProGFL;
    fQATrackHistosRecEffRecKGFL    = o.fQATrackHistosRecEffRecKGFL;
    fQATrackHistosSecRecPi         = o.fQATrackHistosSecRecPi;  
    fQATrackHistosSecRecPro        = o.fQATrackHistosSecRecPro; 
    fQATrackHistosSecRecK          = o.fQATrackHistosSecRecK;
    fQATrackHistosSecRecEl         = o.fQATrackHistosSecRecEl; 
    fQATrackHistosSecRecMu         = o.fQATrackHistosSecRecMu; 
    fQATrackHistosSecRecProGFL     = o.fQATrackHistosSecRecProGFL;
    fQATrackHistosSecRecKGFL       = o.fQATrackHistosSecRecKGFL;
    fFFHistosRecEffRec             = o.fFFHistosRecEffRec;  
    fFFHistosSecRec                = o.fFFHistosSecRec;   
    fFFHistosRecEffRecPi           = o.fFFHistosRecEffRecPi;  
    fFFHistosRecEffRecPro          = o.fFFHistosRecEffRecPro;
    fFFHistosRecEffRecK            = o.fFFHistosRecEffRecK;  
    fFFHistosRecEffRecEl           = o.fFFHistosRecEffRecEl; 
    fFFHistosRecEffRecMu           = o.fFFHistosRecEffRecMu; 
    fFFHistosRecEffRecProGFL       = o.fFFHistosRecEffRecProGFL;
    fFFHistosRecEffRecKGFL         = o.fFFHistosRecEffRecKGFL;
    fFFHistosSecRecPi              = o.fFFHistosSecRecPi;    
    fFFHistosSecRecPro             = o.fFFHistosSecRecPro;   
    fFFHistosSecRecK               = o.fFFHistosSecRecK;     
    fFFHistosSecRecEl              = o.fFFHistosSecRecEl;    
    fFFHistosSecRecMu              = o.fFFHistosSecRecMu;    
    fFFHistosSecRecProGFL          = o.fFFHistosSecRecProGFL;   
    fFFHistosSecRecKGFL            = o.fFFHistosSecRecKGFL;  
    fRandom                        = o.fRandom;
  }
  
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskIDFFTCF::~AliAnalysisTaskIDFFTCF()
{
  // destructor
  
  if(fTracksRecCuts)           delete fTracksRecCuts;
  if(fTracksGen)               delete fTracksGen;
  if(fTracksAODMCCharged)      delete fTracksAODMCCharged;  
  if(fTracksAODMCChargedSec)   delete fTracksAODMCChargedSec;  
  if(fTracksRecQualityCuts)    delete fTracksRecQualityCuts; 
  if(fJetsRec)                 delete fJetsRec;
  if(fJetsRecCuts)             delete fJetsRecCuts;
  if(fJetsGen)                 delete fJetsGen;
  if(fJetsRecEff)              delete fJetsRecEff;

  if(fRandom)               delete fRandom;
}

//______________________________________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliFragFuncHistos::AliFragFuncHistos(const char* name, 
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
AliAnalysisTaskIDFFTCF::AliFragFuncHistos::AliFragFuncHistos(const AliFragFuncHistos& copy)
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
AliAnalysisTaskIDFFTCF::AliFragFuncHistos& AliAnalysisTaskIDFFTCF::AliFragFuncHistos::operator=(const AliAnalysisTaskIDFFTCF::AliFragFuncHistos& o)
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
AliAnalysisTaskIDFFTCF::AliFragFuncHistos::~AliFragFuncHistos()
{
  // destructor 

  if(fh1JetPt)   delete fh1JetPt;
  if(fh2TrackPt) delete fh2TrackPt;
  if(fh2Xi)      delete fh2Xi;
  if(fh2Z)       delete fh2Z;
}

//_________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncHistos::DefineHistos()
{
  // book FF histos

  fh1JetPt   = new TH1F(Form("fh1FFJetPt%s", fNameFF.Data()),"",fNBinsJetPt,fJetPtMin,fJetPtMax);
  fh2TrackPt = new TH2F(Form("fh2FFTrackPt%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsPt, fPtMin, fPtMax);
  fh2Z       = new TH2F(Form("fh2FFZ%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsZ, fZMin, fZMax);
  fh2Xi      = new TH2F(Form("fh2FFXi%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsXi, fXiMin, fXiMax);

  AliAnalysisTaskIDFFTCF::SetProperties(fh1JetPt, "p_{T} [GeV/c]", "entries"); 
  AliAnalysisTaskIDFFTCF::SetProperties(fh2TrackPt,"jet p_{T} [GeV/c]","p_{T} [GeV/c]","entries");
  AliAnalysisTaskIDFFTCF::SetProperties(fh2Xi,"jet p_{T} [GeV/c]","#xi", "entries");
  AliAnalysisTaskIDFFTCF::SetProperties(fh2Z,"jet p_{T} [GeV/c]","z","entries");
}

//_______________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncHistos::FillFF(Float_t trackPt, Float_t jetPt, Bool_t incrementJetPt, Float_t norm, 
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
void AliAnalysisTaskIDFFTCF::AliFragFuncHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  
  list->Add(fh2TrackPt);
  list->Add(fh2Xi);
  list->Add(fh2Z);
}

//_________________________________________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const char* name,
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
AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::AliFragFuncQAJetHistos(const AliFragFuncQAJetHistos& copy)
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
AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos& AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::operator=(const AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos& o)
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
AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::~AliFragFuncQAJetHistos()
{
  // destructor 
  
  if(fh2EtaPhi) delete fh2EtaPhi;
  if(fh1Pt)     delete fh1Pt;
}

//____________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::DefineHistos()
{
  // book jet QA histos

  fh2EtaPhi  = new TH2F(Form("fh2JetQAEtaPhi%s", fNameQAJ.Data()), Form("%s: #eta - #phi distribution", fNameQAJ.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt      = new TH1F(Form("fh1JetQAPt%s", fNameQAJ.Data()), Form("%s: p_{T} distribution", fNameQAJ.Data()), fNBinsPt, fPtMin, fPtMax);
	
  AliAnalysisTaskIDFFTCF::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskIDFFTCF::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
}

//____________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::FillJetQA(Float_t eta, Float_t phi, Float_t pt)
{
  // fill jet QA histos 

  fh2EtaPhi->Fill( eta, phi);
  fh1Pt->Fill( pt );
}

//____________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncQAJetHistos::AddToOutput(TList* list) const 
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh1Pt);
}

//___________________________________________________________________________________________________________
AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const char* name,
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
AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::AliFragFuncQATrackHistos(const AliFragFuncQATrackHistos& copy)
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
AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos& AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::operator=(const AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos& o)
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
AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::~AliFragFuncQATrackHistos()
{
  // destructor 
  
  if(fh2EtaPhi)       delete fh2EtaPhi;
  if(fh2HighPtEtaPhi) delete fh2HighPtEtaPhi;
  if(fh1Pt)           delete fh1Pt;
  if(fh2PhiPt)        delete fh2PhiPt;
}

//______________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::DefineHistos()
{
  // book track QA histos

  fh2EtaPhi       = new TH2F(Form("fh2TrackQAEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh2HighPtEtaPhi = new TH2F(Form("fh2TrackQAHighPtEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution for high-p_{T}", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt           = new TH1F(Form("fh1TrackQAPt%s", fNameQAT.Data()), Form("%s: p_{T} distribution", fNameQAT.Data()), fNBinsPt, fPtMin, fPtMax);
    fh2PhiPt        = new TH2F(Form("fh2TrackQAPhiPt%s", fNameQAT.Data()), Form("%s: #phi - p_{T} distribution", fNameQAT.Data()), fNBinsPhi, fPhiMin, fPhiMax, fNBinsPt, fPtMin, fPtMax);

  AliAnalysisTaskIDFFTCF::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskIDFFTCF::SetProperties(fh2HighPtEtaPhi, "#eta", "#phi");
  AliAnalysisTaskIDFFTCF::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
  AliAnalysisTaskIDFFTCF::SetProperties(fh2PhiPt, "#phi", "p_{T} [GeV/c]"); 
}



//________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::FillTrackQA(Float_t eta, Float_t phi, Float_t pt, Bool_t weightPt, Float_t norm, 
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
void AliAnalysisTaskIDFFTCF::AliFragFuncQATrackHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh2HighPtEtaPhi);
  list->Add(fh1Pt);
  list->Add(fh2PhiPt);
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskIDFFTCF::Notify()
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
void AliAnalysisTaskIDFFTCF::UserCreateOutputObjects()
{
  // create output objects

  if(fDebug > 1) Printf("AliAnalysisTaskIDFFTCF::UserCreateOutputObjects()");
 
  // create list of tracks and jets 

  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE);  

  fTracksGen = new TList();
  fTracksGen->SetOwner(kFALSE);

  fTracksAODMCCharged = new TList();
  fTracksAODMCCharged->SetOwner(kFALSE);
    
  fTracksAODMCChargedSec = new TList();
  fTracksAODMCChargedSec->SetOwner(kFALSE);

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

  fh2PtRecVsGenPrim          = new TH2F("fh2PtRecVsGenPrim","rec vs gen pt",fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax,fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax);
  fh2PtRecVsGenSec           = new TH2F("fh2PtRecVsGenSec","rec vs gen pt",fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax,fQATrackNBinsPt,fQATrackPtMin,fQATrackPtMax);
  

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
    


    fTHnIDFF  = AliIDFFUtils::GetTHn("THnIDFF");
    fTHnIncl  = AliIDFFUtils::GetTHn("THnIncl");
  } // end: FF
  
  // efficiency



    
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
    fCommonHistList->Add(fTHnIDFF);
    fCommonHistList->Add(fTHnIncl);
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
  fCommonHistList->Add(fh1Xsec);
  fCommonHistList->Add(fh1Trials);
  fCommonHistList->Add(fh1PtHard);
  fCommonHistList->Add(fh1PtHardTrials);
 
  if(genJets) fCommonHistList->Add(fh1nGenJets);
 
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

  if(fFFMode){
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsInc,"RecCutsInc",&fFFHistosGenInc,"GenInc");
    BookFFHistos(fCommonHistList,&fFFHistosRecLeadingTrack,"RecLeadingTrack",&fFFHistosGenLeadingTrack,"GenLeadingTrack");
  }
    

  if(fFFMode && genTracks){ 
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsIncPi,"RecCutsInc_piPlusPiMinus",&fFFHistosGenIncPi,"GenInc_piPlusPiMinus");
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsIncPro,"RecCutsInc_ppbar",&fFFHistosGenIncPro,"GenInc_ppbar");
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsIncK,"RecCutsInc_kPlusKMinus",&fFFHistosGenIncK,"GenInc_kPlusKMinus");
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsIncEl,"RecCutsInc_ePlusEMinus",&fFFHistosGenIncEl,"GenInc_ePlusEMinus");
    BookFFHistos(fCommonHistList,&fFFHistosRecCutsIncMu,"RecCutsInc_muPlusMuMinus",&fFFHistosGenIncMu,"GenInc_muPlusMuMinus");
  }


  if(fEffMode && recJetsEff && genTracks){
    if(fQAMode&1){
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRec,"RecEffRec",&fQATrackHistosRecEffGen,"RecEffGen",&fQATrackHistosSecRec,"SecRec");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecPi,"RecEffRec_piPlusPiMinus",&fQATrackHistosRecEffGenPi,"RecEffGen_piPlusPiMinus",&fQATrackHistosSecRecPi,"SecRec_piPlusPiMinus");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecPro,"RecEffRec_ppbar",&fQATrackHistosRecEffGenPro,"RecEffGen_ppbar",&fQATrackHistosSecRecPro,"SecRec_ppbar");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecK,"RecEffRec_kPlusKMinus",&fQATrackHistosRecEffGenK,"RecEffGen_kPlusKMinus",&fQATrackHistosSecRecK,"SecRec_kPlusKMinus");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecEl,"RecEffRec_ePlusEMinus",&fQATrackHistosRecEffGenEl,"RecEffGen_ePlusEMinus",&fQATrackHistosSecRecEl,"SecRec_ePlusEMinus");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecMu,"RecEffRec_muPlusMuMinus",&fQATrackHistosRecEffGenMu,"RecEffGen_muPlusMuMinus",&fQATrackHistosSecRecMu,"SecRec_muPlusMuMinus");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecProGFL,"RecEffRec_ppbarGFL",0,"",&fQATrackHistosSecRecProGFL,"SecRec_ppbarGFL");
      BookQAHistos(fCommonHistList,&fQATrackHistosRecEffRecKGFL,"RecEffRec_kPlusKMinusGFL",0,"",&fQATrackHistosSecRecKGFL,"SecRec_kPlusKMinusGFL");
    }
    if(fFFMode){
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRec,"RecEffRec",0x0,"",&fFFHistosSecRec,"SecRec");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecPi,"RecEffRec_piPlusPiMinus",0x0,"",&fFFHistosSecRecPi,"SecRec_piPlusPiMinus");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecPro,"RecEffRec_ppbar",0x0,"",&fFFHistosSecRecPro,"SecRec_ppbar");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecK,"RecEffRec_kPlusKMinus",0x0,"",&fFFHistosSecRecK,"SecRec_kPlusKMinus");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecEl,"RecEffRec_ePlusEMinus",0x0,"",&fFFHistosSecRecEl,"SecRec_ePlusEMinus");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecMu,"RecEffRec_muPlusMuMinus",0x0,"",&fFFHistosSecRecMu,"SecRec_muPlusMuMinus");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecProGFL,"RecEffRec_ppbarGFL",0x0,"",&fFFHistosSecRecProGFL,"SecRec_ppbarGFL");
      BookFFHistos(fCommonHistList,&fFFHistosRecEffRecKGFL,"RecEffRec_kPlusKMinusGFL",0x0,"",&fFFHistosSecRecKGFL,"SecRec_kPlusKMinusGFL");
    }

    fCommonHistList->Add(fh1nRecEffJets);
    fCommonHistList->Add(fh2PtRecVsGenPrim); 
    fCommonHistList->Add(fh2PtRecVsGenSec); 
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

  if(fStream){
    PostData(2, fTree);
  }
}

//_______________________________________________
void AliAnalysisTaskIDFFTCF::Init()
{
  // Initialization
  if(fDebug > 1) Printf("AliAnalysisTaskIDFFTCF::Init()");

}

//_____________________________________________________________
void AliAnalysisTaskIDFFTCF::UserExec(Option_t *) 
{
  AliIDFFUtils::fPid = 0x0;

  // Main loop
  // Called for each event
  if(fDebug > 1) Printf("AliAnalysisTaskIDFFTCF::UserExec()");
  
  
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
    AliIDFFUtils::fPid = ((AliAODInputHandler*)handler)->GetPIDResponse();
    if(fUseAODInputJets) fAODJets = fAOD;
    if (fDebug > 1)  Printf("%s:%d AOD event from input", (char*)__FILE__,__LINE__);
  }
  else {
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      
      AliInputEventHandler* esdinputHandler = 
	(AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      AliIDFFUtils::fPid = esdinputHandler->GetPIDResponse();
     
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
  if(!AliIDFFUtils::fPid){
    Printf("%s:%d PIDresponse not found", (char*)__FILE__,__LINE__);
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
  
      
  //____ fetch particles __________________________________________________________
  
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
      for(Int_t it=0; it<nRecPartCuts; ++it){
	AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(it));
	if(part)fQATrackHistosRecCuts->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt() );

        // fill inclusive tracks XXX, they have the same track cuts!
        AliAODTrack * inclusiveaod  = dynamic_cast<AliAODTrack*>(fTracksRecCuts->At(it));
        if(inclusiveaod){
          AliIDFFUtils::FillTHn(fTHnIncl, -999, inclusiveaod, fAOD, fTOFCutMode);
        }

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
      Float_t jetPt  = jet->Pt();
      
      if(fQAMode&2) fQAJetHistosRecCuts->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jetPt);
    
      if(fLeadingJets && ij>0) continue;  // leading/all jets
  	
      if(fQAMode&2 && (ij==0) ) fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jetPt );
      
      // get tracks in jet
      TList* jettracklist = new TList();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;

      if(GetFFRadius()<=0)
	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
      else 
	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, 
			     GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
      
      if(GetFFMinNTracks()>0 && jettracklist->GetSize()<=GetFFMinNTracks()) isBadJet = kTRUE;
      
      if(isBadJet) continue; 
      
      for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	if(!trackVP)continue;
	
	AliAODTrack * aodtrack  = dynamic_cast<AliAODTrack*>(jettracklist->At(it));
	if(!aodtrack) continue;
		
	Float_t trackPt = trackVP->Pt();
	
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(fFFMode){
	  
	  fFFHistosRecCutsInc->FillFF(trackPt, jetPt, incrementJetPt); 
	    
	  AliIDFFUtils::FillTHn(fTHnIDFF, jetPt, aodtrack, fAOD, fTOFCutMode);

	  if(it==0){ // leading track, all jets 
	    fFFHistosRecLeadingTrack->FillFF( trackPt, jetPt, kTRUE);
	  }
	} 

      }
      
      delete jettracklist;	
      
    } // end: rec. jets after cuts

    // loop again over rec jets: 
    // PID histos (only for MC) 

    if(fBranchGenJets.Length()>0){ // check if we're running over MC 
      
      for(Int_t ij=0; ij<nRecJetsCuts; ++ij){ // rec jets loop
	
	if(fLeadingJets && ij>0) continue;  // leading/all jets

	TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if(!tca) continue;
	
	AliAODJet* jet = (AliAODJet*)(fJetsRecCuts->At(ij));
	Float_t jetPt  = jet->Pt();
		
	// get tracks in jet
	TList* jettracklist = new TList();
	Double_t sumPt      = 0.;
	Bool_t isBadJet     = kFALSE;
	
	if(GetFFRadius()<=0)
	  GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	else 
	  GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, 
			       GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
      	
	if(GetFFMinNTracks()>0 && jettracklist->GetSize()<=GetFFMinNTracks()) isBadJet = kTRUE;
	
	if(isBadJet) continue; 
	
	Bool_t incrementJetPt = kTRUE; // there could be 0 tracks in jet: first fill jet pt histo once 
	fFFHistosRecCutsIncPi->FillFF(-1, jetPt, incrementJetPt);
	fFFHistosRecCutsIncPro->FillFF(-1, jetPt, incrementJetPt); 
	fFFHistosRecCutsIncK->FillFF(-1, jetPt, incrementJetPt);   
	fFFHistosRecCutsIncEl->FillFF(-1, jetPt, incrementJetPt);  
	fFFHistosRecCutsIncMu->FillFF(-1, jetPt, incrementJetPt);   
      
	incrementJetPt = kFALSE; 

	for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	  
	  AliAODTrack * rectrack  = dynamic_cast<AliAODTrack*>(jettracklist->At(it));
	  if(!rectrack) continue;

	  Int_t label   = TMath::Abs(rectrack->GetLabel());
	  Float_t recPt = rectrack->Pt();
	  
	  // find MC track in our list
	  AliAODMCParticle* gentrack = dynamic_cast<AliAODMCParticle*> (tca->At(label));
	  
	  // Float_t genPt = gentrack->Pt();
	  	    
	  if(gentrack){
	    Int_t mcpdg = TMath::Abs(gentrack->GetPdgCode());
	    
	    if(mcpdg == 211)  fFFHistosRecCutsIncPi->FillFF(recPt, jetPt, incrementJetPt);
	    if(mcpdg == 2212) fFFHistosRecCutsIncPro->FillFF(recPt, jetPt, incrementJetPt);
	    if(mcpdg == 321)  fFFHistosRecCutsIncK->FillFF(recPt, jetPt, incrementJetPt);
	    if(mcpdg == 11)   fFFHistosRecCutsIncEl->FillFF(recPt, jetPt, incrementJetPt);
	    if(mcpdg == 13)   fFFHistosRecCutsIncMu->FillFF(recPt, jetPt, incrementJetPt); 
	  }
	}
	delete jettracklist;	
        
      } // end: rec. jets after cuts
    } // MC 


    // generated jets
    for(Int_t ij=0; ij<nGenJets; ++ij){
      
      AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsGen->At(ij));
      if(!jet)continue;

      Float_t jetPt  = jet->Pt();

      if(fQAMode&2) fQAJetHistosGen->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jetPt);

      if(fLeadingJets && ij>0) continue;  // leading/all jets

      if(fQAMode&2 && (ij==0)) fQAJetHistosGenLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jetPt);

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

      Bool_t incrementJetPt = kTRUE; // first fill jet pt histo once 
      fFFHistosGenInc->FillFF(-1, jetPt, incrementJetPt);
      fFFHistosGenIncPi->FillFF(-1, jetPt, incrementJetPt);
      fFFHistosGenIncPro->FillFF(-1, jetPt, incrementJetPt); 
      fFFHistosGenIncK->FillFF(-1, jetPt, incrementJetPt);   
      fFFHistosGenIncEl->FillFF(-1, jetPt, incrementJetPt);  
      fFFHistosGenIncMu->FillFF(-1, jetPt, incrementJetPt);   

      incrementJetPt = kFALSE;

      for(Int_t it=0; it<jettracklist->GetSize(); ++it){
	
	AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	if(!trackVP)continue;
	
	Float_t trackPt = trackVP->Pt();
	
	if(fFFMode){

	  fFFHistosGenInc->FillFF( trackPt, jetPt, incrementJetPt);

	  AliAODMCParticle* gentrack = dynamic_cast<AliAODMCParticle*>(trackVP);
	  
	  if(gentrack){
	    
	    Int_t mcpdg = TMath::Abs(gentrack->GetPdgCode());
	    
	    if(mcpdg == 211)  fFFHistosGenIncPi->FillFF(trackPt, jetPt, incrementJetPt);
	    if(mcpdg == 2212) fFFHistosGenIncPro->FillFF(trackPt, jetPt, incrementJetPt);
	    if(mcpdg == 321)  fFFHistosGenIncK->FillFF(trackPt, jetPt, incrementJetPt);
	    if(mcpdg == 11)   fFFHistosGenIncEl->FillFF(trackPt, jetPt, incrementJetPt);
	    if(mcpdg == 13)   fFFHistosGenIncMu->FillFF(trackPt, jetPt, incrementJetPt);
	  }
	}
	
	if(it==0){ // leading track, all jets 
	  fFFHistosGenLeadingTrack->FillFF( trackPt, jetPt, kTRUE );	  
	}
      	
      } // track loop 

      delete jettracklist;

    } // gen jets loop     
  } // end: QA, FF and intra-jet
  
  
  // ____ efficiency _______________________________

  if(fEffMode && (fJetTypeRecEff != kJetsUndef)){

    // arrays holding for each generated particle the reconstructed AOD track index & isPrimary flag, are initialized in AssociateGenRec(...) function
    TArrayI indexAODTr; 
    TArrayS isGenPrim; 

    // array holding for each reconstructed AOD track generated particle index, initialized in AssociateGenRec(...) function
    TArrayI indexMCTr; 

    // ... and another set for secondaries from strange/non strange mothers (secondary MC tracks are stored in different lists)
    TArrayI indexAODTrSec; 
    TArrayS isGenSec; 
    TArrayI indexMCTrSec; 
   
    Int_t  nTracksAODMCCharged = GetListOfTracks(fTracksAODMCCharged, kTrackAODMCCharged);
    if(fDebug>2)Printf("%s:%d selected AODMC tracks: %d ",(char*)__FILE__,__LINE__,nTracksAODMCCharged);
  
    Int_t  nTracksAODMCChargedSec = GetListOfTracks(fTracksAODMCChargedSec, kTrackAODMCChargedSec);
    if(fDebug>2)Printf("%s:%d selected AODMC secondary tracks NS: %d ",(char*)__FILE__,__LINE__,nTracksAODMCChargedSec);
  

    Int_t  nTracksRecQualityCuts = GetListOfTracks(fTracksRecQualityCuts, kTrackAODQualityCuts);
    if(fDebug>2)Printf("%s:%d selected rec tracks quality after cuts, full acceptance/pt : %d ",(char*)__FILE__,__LINE__,nTracksRecQualityCuts);
  
    // associate gen and rec tracks, store indices in TArrays 
    AssociateGenRec(fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,indexMCTr,isGenPrim,fh2PtRecVsGenPrim); 
    AssociateGenRec(fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,indexMCTrSec,isGenSec,fh2PtRecVsGenSec);
  
    // single track eff 
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGen,fQATrackHistosRecEffRec,fTracksAODMCCharged,indexAODTr,isGenPrim);
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGenPi,fQATrackHistosRecEffRecPi,fTracksAODMCCharged,indexAODTr,isGenPrim, 211);
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGenPro,fQATrackHistosRecEffRecPro,fTracksAODMCCharged,indexAODTr,isGenPrim, 2212);
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGenK,fQATrackHistosRecEffRecK,fTracksAODMCCharged,indexAODTr,isGenPrim, 321);
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGenEl,fQATrackHistosRecEffRecEl,fTracksAODMCCharged,indexAODTr,isGenPrim, 11);
    if(fQAMode&1) FillSingleTrackHistosRecGen(fQATrackHistosRecEffGenMu,fQATrackHistosRecEffRecMu,fTracksAODMCCharged,indexAODTr,isGenPrim, 13);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0,fQATrackHistosRecEffRecProGFL,fTracksAODMCCharged,indexAODTr,isGenPrim, 2212,kTRUE);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0,fQATrackHistosRecEffRecKGFL,fTracksAODMCCharged,indexAODTr,isGenPrim, 321,kTRUE);

    // secondaries
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRec,fTracksAODMCChargedSec,indexAODTrSec,isGenSec);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecPi,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,211);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecPro,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,2212);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecK,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,321);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecEl,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,11);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecMu,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,13);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecProGFL,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,2212,kTRUE);
    if(fQAMode&1) FillSingleTrackHistosRecGen(0x0,fQATrackHistosSecRecKGFL,fTracksAODMCChargedSec,indexAODTrSec,isGenSec,321,kTRUE);
 

    // jet track eff    
    Double_t sumPtGenLeadingJetRecEff = 0;
    Double_t sumPtGenLeadingJetSec    = 0;
    Double_t sumPtRecLeadingJetRecEff = 0;
    
    for(Int_t ij=0; ij<nRecEffJets; ++ij){ // jet loop 
  
	if(fLeadingJets && ij>0) continue;  // leading/all jets
 
	AliAODJet* jet = (AliAODJet*)(fJetsRecEff->At(ij));
	
	Bool_t isBadJetGenPrim = kFALSE;
	Bool_t isBadJetGenSec  = kFALSE;
	Bool_t isBadJetRec     = kFALSE;
  	
	// for efficiency: gen tracks from pointing with gen/rec jet
	TList* jettracklistGenPrim = new TList();
	
	// if radius<0 -> trackRefs: collect gen tracks in wide radius + fill FF recEff rec histos with tracks contained in track refs
        // note : FF recEff gen histos will be somewhat useless in this approach
	
	if(GetFFRadius() >0)
	  GetJetTracksPointing(fTracksAODMCCharged, jettracklistGenPrim, jet, GetFFRadius(), sumPtGenLeadingJetRecEff, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetGenPrim); 
	else
	  GetJetTracksPointing(fTracksAODMCCharged, jettracklistGenPrim, jet, TMath::Abs(GetFFRadius())+0.2, sumPtGenLeadingJetRecEff, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetGenPrim); 
	
	TList* jettracklistGenSec = new TList();
	if(GetFFRadius() >0)
	  GetJetTracksPointing(fTracksAODMCChargedSec, jettracklistGenSec, jet, GetFFRadius(), sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 
        else
	  GetJetTracksPointing(fTracksAODMCChargedSec, jettracklistGenSec, jet, TMath::Abs(GetFFRadius())+0.2, sumPtGenLeadingJetSec, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetGenSec); 
	
	
	// bin efficiency in jet pt bins using rec tracks  
	TList* jettracklistRec = new TList();
	if(GetFFRadius() >0) GetJetTracksPointing(fTracksRecCuts,jettracklistRec, jet, GetFFRadius(), sumPtRecLeadingJetRecEff, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetRec); 
	else                 GetJetTracksTrackrefs(jettracklistRec, jet, GetFFMinLTrackPt() , GetFFMaxTrackPt(), isBadJetRec); 
	
	
	Double_t jetEta   = jet->Eta();
	Double_t jetPhi   = TVector2::Phi_0_2pi(jet->Phi());
	
	if(GetFFMinNTracks()>0 && jettracklistGenPrim->GetSize()<=GetFFMinNTracks())   isBadJetGenPrim = kTRUE;
	if(GetFFMinNTracks()>0 && jettracklistGenSec->GetSize()<=GetFFMinNTracks())  isBadJetGenSec  = kTRUE;
	if(GetFFMinNTracks()>0 && jettracklistRec->GetSize()<=GetFFMinNTracks())       isBadJetRec     = kTRUE;
	
	if(isBadJetRec) continue;
	
	if(fQAMode&2) fQAJetHistosRecEffLeading->FillJetQA( jetEta, jetPhi, sumPtGenLeadingJetRecEff ); 
	
	if(fFFMode){
	  
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRec,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRec,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec);
	  
 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecPi,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,211); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecPi,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,211);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecPro,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,2212); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecPro,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,2212);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecK,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,321); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecK,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,321);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecEl,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,11); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecEl,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,11);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecMu,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,13); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecMu,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,13);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecProGFL,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,2212,kTRUE); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecProGFL,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,2212,kTRUE);

 	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosRecEffRecKGFL,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,0,321,kTRUE); 
          else                FillJetTrackHistosRec(fFFHistosRecEffRecKGFL,jet,jettracklistGenPrim,fTracksAODMCCharged,fTracksRecQualityCuts,indexAODTr,isGenPrim,jettracklistRec,321,kTRUE);



	  // secondaries: use jet pt from primaries 
	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRec,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0);
	  else                FillJetTrackHistosRec(fFFHistosSecRec,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecPi,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,211);
	  else                FillJetTrackHistosRec(fFFHistosSecRecPi,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,211);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecPro,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,2212); 
	  else                FillJetTrackHistosRec(fFFHistosSecRecPro,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,2212);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecK,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,321); 
	  else                FillJetTrackHistosRec(fFFHistosSecRecK,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,321);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecEl,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,11); 
	  else                FillJetTrackHistosRec(fFFHistosSecRecEl,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,11);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecMu,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,13);
	  else                FillJetTrackHistosRec(fFFHistosSecRecMu,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,13);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecProGFL,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,2212,kTRUE); 
	  else                FillJetTrackHistosRec(fFFHistosSecRecProGFL,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,2212,kTRUE);

	  if(GetFFRadius()>0) FillJetTrackHistosRec(fFFHistosSecRecKGFL,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts, indexAODTrSec,isGenSec,0,321,kTRUE); 
	  else                FillJetTrackHistosRec(fFFHistosSecRecKGFL,jet,jettracklistGenSec,fTracksAODMCChargedSec,fTracksRecQualityCuts,indexAODTrSec,isGenSec,jettracklistRec,321,kTRUE); 
	}
	
	delete jettracklistGenPrim;
	delete jettracklistGenSec;
	delete jettracklistRec;	  
	
    } // jet loop 
  } // eff mode

  
  //___________________
  
  fTracksRecCuts->Clear();
  fTracksGen->Clear();
  fTracksAODMCCharged->Clear();
  fTracksAODMCChargedSec->Clear();
  fTracksRecQualityCuts->Clear();

  fJetsRec->Clear();
  fJetsRecCuts->Clear();
  fJetsGen->Clear();
  fJetsRecEff->Clear();
  
  //Post output data.
  PostData(1, fCommonHistList);
}

//______________________________________________________________
void AliAnalysisTaskIDFFTCF::Terminate(Option_t *) 
{
  // terminated

  if(fDebug > 1) printf("AliAnalysisTaskIDFFTCF::Terminate() \n");
}  

//_________________________________________________________________________________
Int_t AliAnalysisTaskIDFFTCF::GetListOfTracks(TList *list, Int_t type)
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
      
      if(type == kTrackAODCuts || type==kTrackAODQualityCuts){

	if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))   continue;
 
        //new cut on TPC
        if(fTPCCutMode==kPIDN && !AliIDFFUtils::TPCCutPIDN(tr)){
          continue;
        }
        if(fTPCCutMode==kMIGeo && !AliIDFFUtils::TPCCutMIGeo(tr, fAOD, fStream)){
          continue;
        }
	
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
  else if (type==kTrackAODMCCharged || type==kTrackAODMCAll || type==kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSec)  {
    // MC particles (from AOD), physical primaries, all or rather charged or rather charged within acceptance
    if(!fAOD) return -1;
    
    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    
    for(int it=0; it<tca->GetEntriesFast(); ++it){
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
      if(!part)continue;
      if(type != kTrackAODMCChargedSec  && !part->IsPhysicalPrimary())continue;
      if((type == kTrackAODMCChargedSec) && part->IsPhysicalPrimary())continue;

      if (type==kTrackAODMCCharged || type==kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSec){
	if(part->Charge()==0) continue;

	if(type==kTrackAODMCChargedSec){

	  Int_t iMother = part->GetMother();
          if(iMother < 0) continue; // throw out PYTHIA stack partons + incoming protons
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
Int_t AliAnalysisTaskIDFFTCF::GetListOfJets(TList *list, Int_t type)
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
  else{
    if(fDebug>0)Printf("%s:%d no such type %d",(char*)__FILE__,__LINE__,type);
    return 0;
  }
}

// _________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::SetProperties(THnSparse* h, Int_t dim, const char** labels)
{
  // Set properties of THnSparse 

  for(Int_t i=0; i<dim; i++){
    h->GetAxis(i)->SetTitle(labels[i]);
    h->GetAxis(i)->SetTitleColor(1);
  }
}

// __________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::SetProperties(TH1* h,const char* x, const char* y)
{
  //Set properties of histos (x and y title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
}

// _________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::SetProperties(TH1* h,const char* x, const char* y, const char* z)
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
void AliAnalysisTaskIDFFTCF::GetJetTracksPointing(TList* inputlist, TList* outputlist, const AliAODJet* jet, 
						  Double_t radius, Double_t& sumPt, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt)
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
void AliAnalysisTaskIDFFTCF::GetJetTracksTrackrefs(TList* list, const AliAODJet* jet, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt)
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
void  AliAnalysisTaskIDFFTCF::AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,
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
void AliAnalysisTaskIDFFTCF::FillSingleTrackHistosRecGen(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, 
							 const TArrayI& indexAODTr, const TArrayS& isRefGen, Int_t pdg, Bool_t scaleGFL){

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
    if(pdg && TMath::Abs(gentrack->GetPdgCode()) != pdg) continue;  

    if(trackQAGen) trackQAGen->FillTrackQA(etaGen, phiGen, ptGen);

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

    if(iRec>=0 && trackQARec){

      if(scaleGFL){ 
	
	Double_t GFLcorr = 1; 
  	if(gentrack->GetPdgCode() == -2212)      GFLcorr = TrackingPtGeantFlukaCorrectionPrMinus(ptGen);  
  	else if(gentrack->GetPdgCode() == -321)  GFLcorr = TrackingPtGeantFlukaCorrectionKaMinus(ptGen);  

	Double_t weight = 1;
	if(GFLcorr > 0) weight = 1/GFLcorr;
	
	trackQARec->FillTrackQA(etaGen, phiGen, ptGen, kFALSE, 0, kTRUE, weight);
      }
      else trackQARec->FillTrackQA(etaGen, phiGen, ptGen);
    }
  }
}

// ______________________________________________________________________________________________________________________________________________________

void  AliAnalysisTaskIDFFTCF::FillJetTrackHistosRec(AliFragFuncHistos* ffhistRec, AliAODJet* jet, 
						    TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr,
						    const TArrayS& isRefGen, TList* jetTrackListTR, Int_t pdg, Bool_t scaleGFL)
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
    if(pdg && TMath::Abs(gentrack->GetPdgCode()) != pdg) continue;  


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
	
	if(scaleGFL){
	
	  //Double_t weight = 1; 
	  //if(gentrack->GetPdgCode() == -2212)      weight = TrackingPtGeantFlukaCorrectionPrMinus(ptGen);  
	  //else if(gentrack->GetPdgCode() == -321)  weight = TrackingPtGeantFlukaCorrectionKaMinus(ptGen);  
	  //else                                     weight = 1; 
	  
	  Double_t GFLcorr = 1; 
	  if(gentrack->GetPdgCode() == -2212)      GFLcorr = TrackingPtGeantFlukaCorrectionPrMinus(ptGen);  
	  else if(gentrack->GetPdgCode() == -321)  GFLcorr = TrackingPtGeantFlukaCorrectionKaMinus(ptGen);  

	  Double_t weight = 1;
	  if(GFLcorr > 0) weight = 1/GFLcorr;

	  ffhistRec->FillFF( trackPt, jetPtRec, incrementJetPt, 0, kTRUE, weight );
	}
	else  ffhistRec->FillFF( trackPt, jetPtRec, incrementJetPt );
	
	listRecTracks->Add(rectrack);	
      }
    }
  }

  delete listRecTracks;

}

// ______________________________________________________________________________________________________________________________________________________
Float_t AliAnalysisTaskIDFFTCF::CalcJetArea(const Float_t etaJet, const Float_t rc) const
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

//____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::BookQAHistos(TList* list, AliFragFuncQATrackHistos** rec, TString strTitRec, AliFragFuncQATrackHistos** gen, TString strTitGen, 
					  AliFragFuncQATrackHistos** sec, TString strTitSec){
  
  // book QA histos 

  if(strTitRec.Length()>0){
    
    *rec = new AliFragFuncQATrackHistos(strTitRec, fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
					fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
					fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
					fQATrackHighPtThreshold);

    (*rec)->DefineHistos(); 
    (*rec)->AddToOutput(list);  
  }

  if(strTitGen.Length()>0){

    *gen = new AliFragFuncQATrackHistos(strTitGen, fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
				       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
					fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
				       fQATrackHighPtThreshold);
    
    (*gen)->DefineHistos(); 
    (*gen)->AddToOutput(list); 
  }

  if(strTitSec.Length()>0){
    
    *sec = new AliFragFuncQATrackHistos(strTitSec, fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
				       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
				       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
				       fQATrackHighPtThreshold);

    (*sec)->DefineHistos(); 
    (*sec)->AddToOutput(list); 
  }    
  
}

//____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskIDFFTCF::BookFFHistos(TList* list, AliFragFuncHistos** rec, TString strTitRec, AliFragFuncHistos** gen, TString strTitGen, 
					  AliFragFuncHistos** sec, TString strTitSec){
  
  // book FF histos 

  if(strTitRec.Length()>0){

    *rec = new AliFragFuncHistos(strTitRec, fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
				 fFFNBinsPt, fFFPtMin, fFFPtMax, 
				 fFFNBinsXi, fFFXiMin, fFFXiMax,  
				 fFFNBinsZ , fFFZMin , fFFZMax);

    (*rec)->DefineHistos(); 
    (*rec)->AddToOutput(list);
  }
  
  if(strTitGen.Length()>0){
    
    *gen = new AliFragFuncHistos(strTitGen, fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
				 fFFNBinsPt, fFFPtMin, fFFPtMax, 
				 fFFNBinsXi, fFFXiMin, fFFXiMax,  
				 fFFNBinsZ , fFFZMin , fFFZMax);

    (*gen)->DefineHistos(); 
    (*gen)->AddToOutput(list);
  }

  if(strTitSec.Length()>0){
    
    *sec = new AliFragFuncHistos(strTitSec, fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
				 fFFNBinsPt, fFFPtMin, fFFPtMax, 
				 fFFNBinsXi, fFFXiMin, fFFXiMax,  
				 fFFNBinsZ , fFFZMin , fFFZMax);

    (*sec)->DefineHistos(); 
    (*sec)->AddToOutput(list); 
  }
}

//____________________________________________________________________________________
Double_t  AliAnalysisTaskIDFFTCF::TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  // GEANT-FLUKA correction for pbar from Roberto via Benjamin

  Double_t corr = 1. - 0.129758 * TMath::Exp(-pTmc * 0.679612);
  return corr;
}

//____________________________________________________________________________________
Double_t  AliAnalysisTaskIDFFTCF::TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{

  // GEANT-FLUKA correction for K- from Roberto via Benjamin

  Double_t corr = TMath::Min((0.972865 + 0.0117093 * pTmc), 1.);
  return corr; 
}
