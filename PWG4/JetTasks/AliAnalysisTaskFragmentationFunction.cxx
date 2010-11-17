/*************************************************************************
 *                                                                       *
 * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
 *                                                                       *
 *************************************************************************/


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
#include "TString.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom.h"

#include "AliAODInputHandler.h" 
#include "AliAODHandler.h" 
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODJet.h"
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
   ,fBranchGenJets("")
   ,fTrackTypeGen(0)
   ,fJetTypeGen(0)
   ,fJetTypeRecEff(0)
   ,fFilterMask(0)
   ,fUsePhysicsSelection(kTRUE)
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
   ,fUseRecEffRecJetPtBins(1)
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
   ,fCommonHistList(0)
   ,fh1EvtSelection(0)
   ,fh1VertexNContributors(0)
   ,fh1VertexZ(0)
   ,fh1EvtMult(0)
   ,fh1Xsec(0)
   ,fh1Trials(0)
   ,fh1PtHard(0)
   ,fh1PtHardTrials(0)
   ,fh1nRecJetsCuts(0)
   ,fh1nGenJets(0)
   ,fh1nRecEffJets(0)
   ,fh2PtRecVsGenPrim(0)
   ,fQATrackHistosRecEffGen(0)  
   ,fQATrackHistosRecEffRec(0)  
   ,fFFHistosRecEffGen(0)    
   ,fFFHistosRecEffRec(0)   
   // Background
   ,fh1OutLeadingMult(0)
   ,fh1PerpMult(0)
   ,fh1Out2JetsMult(0)
   ,fh1Out3JetsMult(0)
   ,fQABckgNoJetTrackHistosRec(0)      
   ,fQABckgNoJetTrackHistosRecCuts(0)  
   ,fQABckgNoJetTrackHistosGen(0)      
   ,fQABckgLeadingTrackHistosRec(0)      
   ,fQABckgLeadingTrackHistosRecCuts(0)  
   ,fQABckgLeadingTrackHistosGen(0)      
   ,fQABckg2JetsTrackHistosRec(0)      
   ,fQABckg2JetsTrackHistosRecCuts(0)  
   ,fQABckg2JetsTrackHistosGen(0)
   ,fQABckg3JetsTrackHistosRec(0)      
   ,fQABckg3JetsTrackHistosRecCuts(0)  
   ,fQABckg3JetsTrackHistosGen(0)      
   ,fQABckgPerpTrackHistosRec(0)      
   ,fQABckgPerpTrackHistosRecCuts(0)  
   ,fQABckgPerpTrackHistosGen(0)  
   ,fFFBckgNoJetHistosRecCuts(0)
   ,fFFBckgNoJetHistosRecLeading(0)
   ,fFFBckgNoJetHistosGen(0)       
   ,fFFBckgNoJetHistosGenLeading(0)
   ,fFFBckgHistosRecCuts(0)
   ,fFFBckgHistosRecLeading(0)
   ,fFFBckgHistosGen(0)       
   ,fFFBckgHistosGenLeading(0)
   ,fFFBckgLeadingHistosRecCuts(0)
   ,fFFBckgLeadingHistosRecLeading(0)
   ,fFFBckgLeadingHistosGen(0)       
   ,fFFBckgLeadingHistosGenLeading(0)
   ,fFFBckg2JetsHistosRecCuts(0)     
   ,fFFBckg2JetsHistosRecLeading(0)  
   ,fFFBckg2JetsHistosGen(0)         
   ,fFFBckg2JetsHistosGenLeading(0)  
   ,fFFBckg3JetsHistosRecCuts(0)     
   ,fFFBckg3JetsHistosRecLeading(0)  
   ,fFFBckg3JetsHistosGen(0)         
   ,fFFBckg3JetsHistosGenLeading(0)  
   ,fFFBckgPerpHistosRecCuts(0)
   ,fFFBckgPerpHistosRecLeading(0)
   ,fFFBckgPerpHistosGen(0)       
   ,fFFBckgPerpHistosGenLeading(0)
   ,fFFBckgHistosStatRecCuts(0)
   ,fFFBckgHistosStatRecLeading(0)
   ,fFFBckgHistosStatGen(0)       
   ,fFFBckgHistosStatGenLeading(0)
   ,fFFBckgLeadingHistosStatRecCuts(0)
   ,fFFBckgLeadingHistosStatRecLeading(0)
   ,fFFBckgLeadingHistosStatGen(0)       
   ,fFFBckgLeadingHistosStatGenLeading(0)
   ,fFFBckg2JetsHistosStatRecCuts(0)     
   ,fFFBckg2JetsHistosStatRecLeading(0)  
   ,fFFBckg2JetsHistosStatGen(0)         
   ,fFFBckg2JetsHistosStatGenLeading(0)  
   ,fFFBckg3JetsHistosStatRecCuts(0)     
   ,fFFBckg3JetsHistosStatRecLeading(0)  
   ,fFFBckg3JetsHistosStatGen(0)         
   ,fFFBckg3JetsHistosStatGenLeading(0)  
   ,fIJBckgNoJetHistosRecCuts(0)   
   ,fIJBckgNoJetHistosRecLeading(0)
   ,fIJBckgNoJetHistosGen(0)       
   ,fIJBckgNoJetHistosGenLeading(0)
   ,fIJBckgHistosRecCuts(0)   
   ,fIJBckgHistosRecLeading(0)
   ,fIJBckgHistosGen(0)       
   ,fIJBckgHistosGenLeading(0)
   ,fIJBckgLeadingHistosRecCuts(0)   
   ,fIJBckgLeadingHistosRecLeading(0)
   ,fIJBckgLeadingHistosGen(0)       
   ,fIJBckgLeadingHistosGenLeading(0)
   ,fIJBckg2JetsHistosRecCuts(0)     
   ,fIJBckg2JetsHistosRecLeading(0)  
   ,fIJBckg2JetsHistosGen(0)         
   ,fIJBckg2JetsHistosGenLeading(0)  
   ,fIJBckg3JetsHistosRecCuts(0)     
   ,fIJBckg3JetsHistosRecLeading(0)  
   ,fIJBckg3JetsHistosGen(0)         
   ,fIJBckg3JetsHistosGenLeading(0)  
   ,fIJBckgPerpHistosRecCuts(0)   
   ,fIJBckgPerpHistosRecLeading(0)
   ,fIJBckgPerpHistosGen(0)       
   ,fIJBckgPerpHistosGenLeading(0)
   ,fIJBckgHistosStatRecCuts(0)   
   ,fIJBckgHistosStatRecLeading(0)
   ,fIJBckgHistosStatGen(0)       
   ,fIJBckgHistosStatGenLeading(0)
   ,fIJBckgLeadingHistosStatRecCuts(0)   
   ,fIJBckgLeadingHistosStatRecLeading(0)
   ,fIJBckgLeadingHistosStatGen(0)       
   ,fIJBckgLeadingHistosStatGenLeading(0)
   ,fIJBckg2JetsHistosStatRecCuts(0)     
   ,fIJBckg2JetsHistosStatRecLeading(0)  
   ,fIJBckg2JetsHistosStatGen(0)         
   ,fIJBckg2JetsHistosStatGenLeading(0)  
   ,fIJBckg3JetsHistosStatRecCuts(0)     
   ,fIJBckg3JetsHistosStatRecLeading(0)  
   ,fIJBckg3JetsHistosStatGen(0)         
   ,fIJBckg3JetsHistosStatGenLeading(0) 
{
   // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fAOD(0)
  ,fBranchRecJets("jets")
  ,fBranchGenJets("")
  ,fTrackTypeGen(0)
  ,fJetTypeGen(0)
  ,fJetTypeRecEff(0)
  ,fFilterMask(0)
  ,fUsePhysicsSelection(kTRUE)
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
  ,fUseRecEffRecJetPtBins(1)
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
  ,fCommonHistList(0)
  ,fh1EvtSelection(0)
  ,fh1VertexNContributors(0)
  ,fh1VertexZ(0)
  ,fh1EvtMult(0)
  ,fh1Xsec(0)
  ,fh1Trials(0)
  ,fh1PtHard(0)
  ,fh1PtHardTrials(0)
  ,fh1nRecJetsCuts(0)
  ,fh1nGenJets(0)
  ,fh1nRecEffJets(0)
  ,fh2PtRecVsGenPrim(0)
  ,fQATrackHistosRecEffGen(0)  
  ,fQATrackHistosRecEffRec(0)  
  ,fFFHistosRecEffGen(0)    
  ,fFFHistosRecEffRec(0)   
   // Background
   ,fh1OutLeadingMult(0)
   ,fh1PerpMult(0)
   ,fh1Out2JetsMult(0)
   ,fh1Out3JetsMult(0)
   ,fQABckgNoJetTrackHistosRec(0)      
   ,fQABckgNoJetTrackHistosRecCuts(0)  
   ,fQABckgNoJetTrackHistosGen(0)      
   ,fQABckgLeadingTrackHistosRec(0)      
   ,fQABckgLeadingTrackHistosRecCuts(0)  
   ,fQABckgLeadingTrackHistosGen(0)      
   ,fQABckg2JetsTrackHistosRec(0)      
   ,fQABckg2JetsTrackHistosRecCuts(0)  
   ,fQABckg2JetsTrackHistosGen(0)
   ,fQABckg3JetsTrackHistosRec(0)      
   ,fQABckg3JetsTrackHistosRecCuts(0)  
   ,fQABckg3JetsTrackHistosGen(0)      
   ,fQABckgPerpTrackHistosRec(0)      
   ,fQABckgPerpTrackHistosRecCuts(0)  
   ,fQABckgPerpTrackHistosGen(0)  
   ,fFFBckgNoJetHistosRecCuts(0)
   ,fFFBckgNoJetHistosRecLeading(0)
   ,fFFBckgNoJetHistosGen(0)       
   ,fFFBckgNoJetHistosGenLeading(0)
   ,fFFBckgHistosRecCuts(0)
   ,fFFBckgHistosRecLeading(0)
   ,fFFBckgHistosGen(0)       
   ,fFFBckgHistosGenLeading(0)
   ,fFFBckgLeadingHistosRecCuts(0)
   ,fFFBckgLeadingHistosRecLeading(0)
   ,fFFBckgLeadingHistosGen(0)       
   ,fFFBckgLeadingHistosGenLeading(0)
   ,fFFBckg2JetsHistosRecCuts(0)     
   ,fFFBckg2JetsHistosRecLeading(0)  
   ,fFFBckg2JetsHistosGen(0)         
   ,fFFBckg2JetsHistosGenLeading(0)  
   ,fFFBckg3JetsHistosRecCuts(0)     
   ,fFFBckg3JetsHistosRecLeading(0)  
   ,fFFBckg3JetsHistosGen(0)         
   ,fFFBckg3JetsHistosGenLeading(0)  
   ,fFFBckgPerpHistosRecCuts(0)
   ,fFFBckgPerpHistosRecLeading(0)
   ,fFFBckgPerpHistosGen(0)       
   ,fFFBckgPerpHistosGenLeading(0)
   ,fFFBckgHistosStatRecCuts(0)
   ,fFFBckgHistosStatRecLeading(0)
   ,fFFBckgHistosStatGen(0)       
   ,fFFBckgHistosStatGenLeading(0)
   ,fFFBckgLeadingHistosStatRecCuts(0)
   ,fFFBckgLeadingHistosStatRecLeading(0)
   ,fFFBckgLeadingHistosStatGen(0)       
   ,fFFBckgLeadingHistosStatGenLeading(0)
   ,fFFBckg2JetsHistosStatRecCuts(0)     
   ,fFFBckg2JetsHistosStatRecLeading(0)  
   ,fFFBckg2JetsHistosStatGen(0)         
   ,fFFBckg2JetsHistosStatGenLeading(0)  
   ,fFFBckg3JetsHistosStatRecCuts(0)     
   ,fFFBckg3JetsHistosStatRecLeading(0)  
   ,fFFBckg3JetsHistosStatGen(0)         
   ,fFFBckg3JetsHistosStatGenLeading(0)  
   ,fIJBckgNoJetHistosRecCuts(0)   
   ,fIJBckgNoJetHistosRecLeading(0)
   ,fIJBckgNoJetHistosGen(0)       
   ,fIJBckgNoJetHistosGenLeading(0)
   ,fIJBckgHistosRecCuts(0)   
   ,fIJBckgHistosRecLeading(0)
   ,fIJBckgHistosGen(0)       
   ,fIJBckgHistosGenLeading(0)
   ,fIJBckgLeadingHistosRecCuts(0)   
   ,fIJBckgLeadingHistosRecLeading(0)
   ,fIJBckgLeadingHistosGen(0)       
   ,fIJBckgLeadingHistosGenLeading(0)
   ,fIJBckg2JetsHistosRecCuts(0)     
   ,fIJBckg2JetsHistosRecLeading(0)  
   ,fIJBckg2JetsHistosGen(0)         
   ,fIJBckg2JetsHistosGenLeading(0)  
   ,fIJBckg3JetsHistosRecCuts(0)     
   ,fIJBckg3JetsHistosRecLeading(0)  
   ,fIJBckg3JetsHistosGen(0)         
   ,fIJBckg3JetsHistosGenLeading(0)  
   ,fIJBckgPerpHistosRecCuts(0)   
   ,fIJBckgPerpHistosRecLeading(0)
   ,fIJBckgPerpHistosGen(0)       
   ,fIJBckgPerpHistosGenLeading(0)
   ,fIJBckgHistosStatRecCuts(0)   
   ,fIJBckgHistosStatRecLeading(0)
   ,fIJBckgHistosStatGen(0)       
   ,fIJBckgHistosStatGenLeading(0)
   ,fIJBckgLeadingHistosStatRecCuts(0)   
   ,fIJBckgLeadingHistosStatRecLeading(0)
   ,fIJBckgLeadingHistosStatGen(0)       
   ,fIJBckgLeadingHistosStatGenLeading(0)
   ,fIJBckg2JetsHistosStatRecCuts(0)     
   ,fIJBckg2JetsHistosStatRecLeading(0)  
   ,fIJBckg2JetsHistosStatGen(0)         
   ,fIJBckg2JetsHistosStatGenLeading(0)  
   ,fIJBckg3JetsHistosStatRecCuts(0)     
   ,fIJBckg3JetsHistosStatRecLeading(0)  
   ,fIJBckg3JetsHistosStatGen(0)         
   ,fIJBckg3JetsHistosStatGenLeading(0) 
{
  // constructor
  
  DefineOutput(1,TList::Class());
  

}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const  AliAnalysisTaskFragmentationFunction &copy)
  : AliAnalysisTaskSE()
  ,fESD(copy.fESD)
  ,fAOD(copy.fAOD)
  ,fBranchRecJets(copy.fBranchRecJets)
  ,fBranchGenJets(copy.fBranchGenJets)
  ,fTrackTypeGen(copy.fTrackTypeGen)
  ,fJetTypeGen(copy.fJetTypeGen)
  ,fJetTypeRecEff(copy.fJetTypeRecEff)
  ,fFilterMask(copy.fFilterMask)
  ,fUsePhysicsSelection(copy.fUsePhysicsSelection)
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
  ,fUseRecEffRecJetPtBins(copy.fUseRecEffRecJetPtBins)
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
  ,fCommonHistList(copy.fCommonHistList)
  ,fh1EvtSelection(copy.fh1EvtSelection)
  ,fh1VertexNContributors(copy.fh1VertexNContributors)
  ,fh1VertexZ(copy.fh1VertexZ)
  ,fh1EvtMult(copy.fh1EvtMult)
  ,fh1Xsec(copy.fh1Xsec)
  ,fh1Trials(copy.fh1Trials)
  ,fh1PtHard(copy.fh1PtHard)  
  ,fh1PtHardTrials(copy.fh1PtHardTrials)  
  ,fh1nRecJetsCuts(copy.fh1nRecJetsCuts)
  ,fh1nGenJets(copy.fh1nGenJets)
  ,fh1nRecEffJets(copy.fh1nRecEffJets)
  ,fh2PtRecVsGenPrim(copy.fh2PtRecVsGenPrim)
  ,fQATrackHistosRecEffGen(copy.fQATrackHistosRecEffGen)  
  ,fQATrackHistosRecEffRec(copy.fQATrackHistosRecEffRec)  
  ,fFFHistosRecEffGen(copy.fFFHistosRecEffGen)    
  ,fFFHistosRecEffRec(copy.fFFHistosRecEffRec)   
  // Background
  ,fh1OutLeadingMult(copy.fh1OutLeadingMult)
  ,fh1PerpMult(copy.fh1PerpMult)
  ,fh1Out2JetsMult(copy.fh1Out2JetsMult)
  ,fh1Out3JetsMult(copy.fh1Out3JetsMult)
  ,fQABckgNoJetTrackHistosRec(copy.fQABckgNoJetTrackHistosRec)      
  ,fQABckgNoJetTrackHistosRecCuts(copy.fQABckgNoJetTrackHistosRecCuts)  
  ,fQABckgNoJetTrackHistosGen(copy.fQABckgNoJetTrackHistosGen)      
  ,fQABckgLeadingTrackHistosRec(copy.fQABckgLeadingTrackHistosRec)      
  ,fQABckgLeadingTrackHistosRecCuts(copy.fQABckgLeadingTrackHistosRecCuts)  
  ,fQABckgLeadingTrackHistosGen(copy.fQABckgLeadingTrackHistosGen)      
  ,fQABckg2JetsTrackHistosRec(copy.fQABckg2JetsTrackHistosRec)      
  ,fQABckg2JetsTrackHistosRecCuts(copy.fQABckg2JetsTrackHistosRecCuts)  
  ,fQABckg2JetsTrackHistosGen(copy.fQABckg2JetsTrackHistosGen)
  ,fQABckg3JetsTrackHistosRec(copy.fQABckg3JetsTrackHistosRec)      
  ,fQABckg3JetsTrackHistosRecCuts(copy.fQABckg3JetsTrackHistosRecCuts)  
  ,fQABckg3JetsTrackHistosGen(copy.fQABckg3JetsTrackHistosGen)      
  ,fQABckgPerpTrackHistosRec(copy.fQABckgPerpTrackHistosRec)      
  ,fQABckgPerpTrackHistosRecCuts(copy.fQABckgPerpTrackHistosRecCuts)  
  ,fQABckgPerpTrackHistosGen(copy.fQABckgPerpTrackHistosGen)      
   ,fFFBckgNoJetHistosRecCuts(copy.fFFBckgNoJetHistosRecCuts)
  ,fFFBckgNoJetHistosRecLeading(copy.fFFBckgNoJetHistosRecLeading)
  ,fFFBckgNoJetHistosGen(copy.fFFBckgNoJetHistosGen)       
  ,fFFBckgNoJetHistosGenLeading(copy.fFFBckgNoJetHistosGenLeading)
  ,fFFBckgHistosRecCuts(copy.fFFBckgHistosRecCuts)
  ,fFFBckgHistosRecLeading(copy.fFFBckgHistosRecLeading)
  ,fFFBckgHistosGen(copy.fFFBckgHistosGen)       
  ,fFFBckgHistosGenLeading(copy.fFFBckgHistosGenLeading)
  ,fFFBckgLeadingHistosRecCuts(copy.fFFBckgLeadingHistosRecCuts)
  ,fFFBckgLeadingHistosRecLeading(copy.fFFBckgLeadingHistosRecLeading)
  ,fFFBckgLeadingHistosGen(copy.fFFBckgLeadingHistosGen)       
  ,fFFBckgLeadingHistosGenLeading(copy.fFFBckgLeadingHistosGenLeading)
  ,fFFBckg2JetsHistosRecCuts(copy.fFFBckg2JetsHistosRecCuts)     
  ,fFFBckg2JetsHistosRecLeading(copy.fFFBckg2JetsHistosRecLeading)  
  ,fFFBckg2JetsHistosGen(copy.fFFBckg2JetsHistosGen)         
  ,fFFBckg2JetsHistosGenLeading(copy.fFFBckg2JetsHistosGenLeading)  
  ,fFFBckg3JetsHistosRecCuts(copy.fFFBckg3JetsHistosRecCuts)     
  ,fFFBckg3JetsHistosRecLeading(copy.fFFBckg3JetsHistosRecLeading)  
  ,fFFBckg3JetsHistosGen(copy.fFFBckg3JetsHistosGen)         
  ,fFFBckg3JetsHistosGenLeading(copy.fFFBckg3JetsHistosGenLeading)  
  ,fFFBckgPerpHistosRecCuts(copy.fFFBckgPerpHistosRecCuts)     
  ,fFFBckgPerpHistosRecLeading(copy.fFFBckgPerpHistosRecLeading)  
  ,fFFBckgPerpHistosGen(copy.fFFBckgPerpHistosGen)         
  ,fFFBckgPerpHistosGenLeading(copy.fFFBckgPerpHistosGenLeading) 
  ,fFFBckgHistosStatRecCuts(copy.fFFBckgHistosStatRecCuts)
  ,fFFBckgHistosStatRecLeading(copy.fFFBckgHistosStatRecLeading)
  ,fFFBckgHistosStatGen(copy.fFFBckgHistosStatGen)       
  ,fFFBckgHistosStatGenLeading(copy.fFFBckgHistosStatGenLeading)
  ,fFFBckgLeadingHistosStatRecCuts(copy.fFFBckgLeadingHistosStatRecCuts)
  ,fFFBckgLeadingHistosStatRecLeading(copy.fFFBckgLeadingHistosStatRecLeading)
  ,fFFBckgLeadingHistosStatGen(copy.fFFBckgLeadingHistosStatGen)       
  ,fFFBckgLeadingHistosStatGenLeading(copy.fFFBckgLeadingHistosStatGenLeading)
  ,fFFBckg2JetsHistosStatRecCuts(copy.fFFBckg2JetsHistosStatRecCuts)     
  ,fFFBckg2JetsHistosStatRecLeading(copy.fFFBckg2JetsHistosStatRecLeading)  
  ,fFFBckg2JetsHistosStatGen(copy.fFFBckg2JetsHistosStatGen)         
  ,fFFBckg2JetsHistosStatGenLeading(copy.fFFBckg2JetsHistosStatGenLeading)  
  ,fFFBckg3JetsHistosStatRecCuts(copy.fFFBckg3JetsHistosStatRecCuts)     
  ,fFFBckg3JetsHistosStatRecLeading(copy.fFFBckg3JetsHistosStatRecLeading)  
  ,fFFBckg3JetsHistosStatGen(copy.fFFBckg3JetsHistosStatGen)         
  ,fFFBckg3JetsHistosStatGenLeading(copy.fFFBckg3JetsHistosStatGenLeading)  
  ,fIJBckgNoJetHistosRecCuts(copy.fIJBckgNoJetHistosRecCuts)   
  ,fIJBckgNoJetHistosRecLeading(copy.fIJBckgNoJetHistosRecLeading)
  ,fIJBckgNoJetHistosGen(copy.fIJBckgNoJetHistosGen)       
  ,fIJBckgNoJetHistosGenLeading(copy.fIJBckgNoJetHistosGenLeading)
  ,fIJBckgHistosRecCuts(copy.fIJBckgHistosRecCuts)   
  ,fIJBckgHistosRecLeading(copy.fIJBckgHistosRecLeading)
  ,fIJBckgHistosGen(copy.fIJBckgHistosGen)       
  ,fIJBckgHistosGenLeading(copy.fIJBckgHistosGenLeading)
  ,fIJBckgLeadingHistosRecCuts(copy.fIJBckgLeadingHistosRecCuts)   
  ,fIJBckgLeadingHistosRecLeading(copy.fIJBckgLeadingHistosRecLeading)
  ,fIJBckgLeadingHistosGen(copy.fIJBckgLeadingHistosGen)       
  ,fIJBckgLeadingHistosGenLeading(copy.fIJBckgLeadingHistosGenLeading)
  ,fIJBckg2JetsHistosRecCuts(copy.fIJBckg2JetsHistosRecCuts)     
  ,fIJBckg2JetsHistosRecLeading(copy.fIJBckg2JetsHistosRecLeading)  
  ,fIJBckg2JetsHistosGen(copy.fIJBckg2JetsHistosGen)         
  ,fIJBckg2JetsHistosGenLeading(copy.fIJBckg2JetsHistosGenLeading)  
  ,fIJBckg3JetsHistosRecCuts(copy.fIJBckg3JetsHistosRecCuts)     
  ,fIJBckg3JetsHistosRecLeading(copy.fIJBckg3JetsHistosRecLeading)  
  ,fIJBckg3JetsHistosGen(copy.fIJBckg3JetsHistosGen)         
  ,fIJBckg3JetsHistosGenLeading(copy.fIJBckg3JetsHistosGenLeading) 
  ,fIJBckgPerpHistosRecCuts(copy.fIJBckgPerpHistosRecCuts)   
  ,fIJBckgPerpHistosRecLeading(copy.fIJBckgPerpHistosRecLeading)
  ,fIJBckgPerpHistosGen(copy.fIJBckgPerpHistosGen)       
  ,fIJBckgPerpHistosGenLeading(copy.fIJBckgPerpHistosGenLeading)
  ,fIJBckgHistosStatRecCuts(copy.fIJBckgHistosStatRecCuts)   
  ,fIJBckgHistosStatRecLeading(copy.fIJBckgHistosStatRecLeading)
  ,fIJBckgHistosStatGen(copy.fIJBckgHistosStatGen)       
  ,fIJBckgHistosStatGenLeading(copy.fIJBckgHistosStatGenLeading)
  ,fIJBckgLeadingHistosStatRecCuts(copy.fIJBckgLeadingHistosStatRecCuts)   
  ,fIJBckgLeadingHistosStatRecLeading(copy.fIJBckgLeadingHistosStatRecLeading)
  ,fIJBckgLeadingHistosStatGen(copy.fIJBckgLeadingHistosStatGen)       
  ,fIJBckgLeadingHistosStatGenLeading(copy.fIJBckgLeadingHistosStatGenLeading)
  ,fIJBckg2JetsHistosStatRecCuts(copy.fIJBckg2JetsHistosStatRecCuts)     
  ,fIJBckg2JetsHistosStatRecLeading(copy.fIJBckg2JetsHistosStatRecLeading)  
  ,fIJBckg2JetsHistosStatGen(copy.fIJBckg2JetsHistosStatGen)         
  ,fIJBckg2JetsHistosStatGenLeading(copy.fIJBckg2JetsHistosStatGenLeading)  
  ,fIJBckg3JetsHistosStatRecCuts(copy.fIJBckg3JetsHistosStatRecCuts)     
  ,fIJBckg3JetsHistosStatRecLeading(copy.fIJBckg3JetsHistosStatRecLeading)  
  ,fIJBckg3JetsHistosStatGen(copy.fIJBckg3JetsHistosStatGen)         
  ,fIJBckg3JetsHistosStatGenLeading(copy.fIJBckg3JetsHistosStatGenLeading) 
{
  // copy constructor

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
    fBranchGenJets                = o.fBranchGenJets;
    fTrackTypeGen                 = o.fTrackTypeGen;
    fJetTypeGen                   = o.fJetTypeGen;
    fJetTypeRecEff                = o.fJetTypeRecEff;
    fFilterMask                   = o.fFilterMask;
    fUsePhysicsSelection          = o.fUsePhysicsSelection;
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
    fUseRecEffRecJetPtBins        = o.fUseRecEffRecJetPtBins;
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
    fCommonHistList               = o.fCommonHistList;
    fh1EvtSelection               = o.fh1EvtSelection;
    fh1VertexNContributors        = o.fh1VertexNContributors;
    fh1VertexZ                    = o.fh1VertexZ;
    fh1EvtMult                    = o.fh1EvtMult;
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
    // Background
    fh1OutLeadingMult                  = o.fh1OutLeadingMult;
    fh1PerpMult                        = o.fh1PerpMult;
    fh1Out2JetsMult                    = o.fh1Out2JetsMult;
    fh1Out3JetsMult                    = o.fh1Out3JetsMult;
    fQABckgNoJetTrackHistosRec         = o.fQABckgNoJetTrackHistosRec;      
    fQABckgNoJetTrackHistosRecCuts     = o.fQABckgNoJetTrackHistosRecCuts;  
    fQABckgNoJetTrackHistosGen         = o.fQABckgNoJetTrackHistosGen;      
    fQABckgLeadingTrackHistosRec       = o.fQABckgLeadingTrackHistosRec;      
    fQABckgLeadingTrackHistosRecCuts   = o.fQABckgLeadingTrackHistosRecCuts;  
    fQABckgLeadingTrackHistosGen       = o.fQABckgLeadingTrackHistosGen;      
    fQABckg2JetsTrackHistosRec         = o.fQABckg2JetsTrackHistosRec;      
    fQABckg2JetsTrackHistosRecCuts     = o.fQABckg2JetsTrackHistosRecCuts;  
    fQABckg2JetsTrackHistosGen         = o.fQABckg2JetsTrackHistosGen;
    fQABckg3JetsTrackHistosRec         = o.fQABckg3JetsTrackHistosRec;      
    fQABckg3JetsTrackHistosRecCuts     = o.fQABckg3JetsTrackHistosRecCuts;  
    fQABckg3JetsTrackHistosGen         = o.fQABckg3JetsTrackHistosGen;      
    fQABckgPerpTrackHistosRec          = o.fQABckgPerpTrackHistosRec;      
    fQABckgPerpTrackHistosRecCuts      = o.fQABckgPerpTrackHistosRecCuts;  
    fQABckgPerpTrackHistosGen          = o.fQABckgPerpTrackHistosGen;      
    fFFBckgNoJetHistosRecCuts          = o.fFFBckgNoJetHistosRecCuts;
    fFFBckgNoJetHistosRecLeading       = o.fFFBckgNoJetHistosRecLeading;
    fFFBckgNoJetHistosGen              = o.fFFBckgNoJetHistosGen;       
    fFFBckgNoJetHistosGenLeading       = o.fFFBckgNoJetHistosGenLeading;
    fFFBckgHistosRecCuts               = o.fFFBckgHistosRecCuts;
    fFFBckgHistosRecLeading            = o.fFFBckgHistosRecLeading;
    fFFBckgHistosGen                   = o.fFFBckgHistosGen;       
    fFFBckgHistosGenLeading            = o.fFFBckgHistosGenLeading;
    fFFBckgLeadingHistosRecCuts        = o.fFFBckgLeadingHistosRecCuts;
    fFFBckgLeadingHistosRecLeading     = o.fFFBckgLeadingHistosRecLeading;
    fFFBckgLeadingHistosGen            = o.fFFBckgLeadingHistosGen;       
    fFFBckgLeadingHistosGenLeading     = o.fFFBckgLeadingHistosGenLeading;
    fFFBckg2JetsHistosRecCuts          = o.fFFBckg2JetsHistosRecCuts;     
    fFFBckg2JetsHistosRecLeading       = o.fFFBckg2JetsHistosRecLeading;  
    fFFBckg2JetsHistosGen              = o.fFFBckg2JetsHistosGen;         
    fFFBckg2JetsHistosGenLeading       = o.fFFBckg2JetsHistosGenLeading;  
    fFFBckg3JetsHistosRecCuts          = o.fFFBckg3JetsHistosRecCuts;     
    fFFBckg3JetsHistosRecLeading       = o.fFFBckg3JetsHistosRecLeading;  
    fFFBckg3JetsHistosGen              = o.fFFBckg3JetsHistosGen;         
    fFFBckg3JetsHistosGenLeading       = o.fFFBckg3JetsHistosGenLeading;  
    fFFBckgPerpHistosRecCuts           = o.fFFBckgPerpHistosRecCuts;     
    fFFBckgPerpHistosRecLeading        = o.fFFBckgPerpHistosRecLeading;  
    fFFBckgPerpHistosGen               = o.fFFBckgPerpHistosGen;         
    fFFBckgPerpHistosGenLeading        = o.fFFBckgPerpHistosGenLeading;
    fFFBckgHistosStatRecCuts           = o.fFFBckgHistosStatRecCuts;
    fFFBckgHistosStatRecLeading        = o.fFFBckgHistosStatRecLeading;
    fFFBckgHistosStatGen               = o.fFFBckgHistosStatGen;       
    fFFBckgHistosStatGenLeading        = o.fFFBckgHistosStatGenLeading;
    fFFBckgLeadingHistosStatRecCuts    = o.fFFBckgLeadingHistosStatRecCuts;
    fFFBckgLeadingHistosStatRecLeading = o.fFFBckgLeadingHistosStatRecLeading;
    fFFBckgLeadingHistosStatGen        = o.fFFBckgLeadingHistosStatGen;       
    fFFBckgLeadingHistosStatGenLeading = o.fFFBckgLeadingHistosStatGenLeading;
    fFFBckg2JetsHistosStatRecCuts      = o.fFFBckg2JetsHistosStatRecCuts;     
    fFFBckg2JetsHistosStatRecLeading   = o.fFFBckg2JetsHistosStatRecLeading;  
    fFFBckg2JetsHistosStatGen          = o.fFFBckg2JetsHistosStatGen;         
    fFFBckg2JetsHistosStatGenLeading   = o.fFFBckg2JetsHistosStatGenLeading;  
    fFFBckg3JetsHistosStatRecCuts      = o.fFFBckg3JetsHistosStatRecCuts;     
    fFFBckg3JetsHistosStatRecLeading   = o.fFFBckg3JetsHistosStatRecLeading;  
    fFFBckg3JetsHistosStatGen          = o.fFFBckg3JetsHistosStatGen;         
    fFFBckg3JetsHistosStatGenLeading   = o.fFFBckg3JetsHistosStatGenLeading;  
    fIJBckgNoJetHistosRecCuts          = o.fIJBckgNoJetHistosRecCuts;   
    fIJBckgNoJetHistosRecLeading       = o.fIJBckgNoJetHistosRecLeading;
    fIJBckgNoJetHistosGen              = o.fIJBckgNoJetHistosGen;       
    fIJBckgNoJetHistosGenLeading       = o.fIJBckgNoJetHistosGenLeading;
    fIJBckgHistosRecCuts               = o.fIJBckgHistosRecCuts;   
    fIJBckgHistosRecLeading            = o.fIJBckgHistosRecLeading;
    fIJBckgHistosGen                   = o.fIJBckgHistosGen;       
    fIJBckgHistosGenLeading            = o.fIJBckgHistosGenLeading;
    fIJBckgLeadingHistosRecCuts        = o.fIJBckgLeadingHistosRecCuts;   
    fIJBckgLeadingHistosRecLeading     = o.fIJBckgLeadingHistosRecLeading;
    fIJBckgLeadingHistosGen            = o.fIJBckgLeadingHistosGen;       
    fIJBckgLeadingHistosGenLeading     = o.fIJBckgLeadingHistosGenLeading;
    fIJBckg2JetsHistosRecCuts          = o.fIJBckg2JetsHistosRecCuts;     
    fIJBckg2JetsHistosRecLeading       = o.fIJBckg2JetsHistosRecLeading;  
    fIJBckg2JetsHistosGen              = o.fIJBckg2JetsHistosGen;         
    fIJBckg2JetsHistosGenLeading       = o.fIJBckg2JetsHistosGenLeading;  
    fIJBckg3JetsHistosRecCuts          = o.fIJBckg3JetsHistosRecCuts;     
    fIJBckg3JetsHistosRecLeading       = o.fIJBckg3JetsHistosRecLeading;  
    fIJBckg3JetsHistosGen              = o.fIJBckg3JetsHistosGen;         
    fIJBckg3JetsHistosGenLeading       = o.fIJBckg3JetsHistosGenLeading;  
    fIJBckgPerpHistosRecCuts           = o.fIJBckgPerpHistosRecCuts;     
    fIJBckgPerpHistosRecLeading        = o.fIJBckgPerpHistosRecLeading;  
    fIJBckgPerpHistosGen               = o.fIJBckgPerpHistosGen;         
    fIJBckgPerpHistosGenLeading        = o.fIJBckgPerpHistosGenLeading;  
    fIJBckgHistosStatRecCuts           = o.fIJBckgHistosStatRecCuts;   
    fIJBckgHistosStatRecLeading        = o.fIJBckgHistosStatRecLeading;
    fIJBckgHistosStatGen               = o.fIJBckgHistosStatGen;       
    fIJBckgHistosStatGenLeading        = o.fIJBckgHistosStatGenLeading;
    fIJBckgLeadingHistosStatRecCuts    = o.fIJBckgLeadingHistosStatRecCuts;   
    fIJBckgLeadingHistosStatRecLeading = o.fIJBckgLeadingHistosStatRecLeading;
    fIJBckgLeadingHistosStatGen        = o.fIJBckgLeadingHistosStatGen;       
    fIJBckgLeadingHistosStatGenLeading = o.fIJBckgLeadingHistosStatGenLeading;
    fIJBckg2JetsHistosStatRecCuts      = o.fIJBckg2JetsHistosStatRecCuts;     
    fIJBckg2JetsHistosStatRecLeading   = o.fIJBckg2JetsHistosStatRecLeading;  
    fIJBckg2JetsHistosStatGen          = o.fIJBckg2JetsHistosStatGen;         
    fIJBckg2JetsHistosStatGenLeading   = o.fIJBckg2JetsHistosStatGenLeading;  
    fIJBckg3JetsHistosStatRecCuts      = o.fIJBckg3JetsHistosStatRecCuts;     
    fIJBckg3JetsHistosStatRecLeading   = o.fIJBckg3JetsHistosStatRecLeading;  
    fIJBckg3JetsHistosStatGen          = o.fIJBckg3JetsHistosStatGen;         
    fIJBckg3JetsHistosStatGenLeading   = o.fIJBckg3JetsHistosStatGenLeading;  
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
  if(incrementJetPt) fh1JetPt->Fill(jetPt);
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
}

//______________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::DefineHistos()
{
  // book track QA histos

  fh2EtaPhi       = new TH2F(Form("fh2TrackQAEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh2HighPtEtaPhi = new TH2F(Form("fh2TrackQAHighPtEtaPhi%s", fNameQAT.Data()), Form("%s: #eta - #phi distribution for high-p_{T}", fNameQAT.Data()), fNBinsEta, fEtaMin, fEtaMax, fNBinsPhi, fPhiMin, fPhiMax);
  fh1Pt           = new TH1F(Form("fh1TrackQAPt%s", fNameQAT.Data()), Form("%s: p_{T} distribution", fNameQAT.Data()), fNBinsPt, fPtMin, fPtMax);
  
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2EtaPhi, "#eta", "#phi"); 
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2HighPtEtaPhi, "#eta", "#phi");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh1Pt, "p_{T} [GeV/c]", "entries");
}

//________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::FillTrackQA(Float_t eta, Float_t phi, Float_t pt)
{
  // fill track QA histos
    
  fh2EtaPhi->Fill( eta, phi);
  if(pt > fHighPtThreshold) fh2HighPtEtaPhi->Fill( eta, phi);
  fh1Pt->Fill( pt );	
}

//______________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQATrackHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2EtaPhi);
  list->Add(fh2HighPtEtaPhi);
  list->Add(fh1Pt);
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
  ,fh2Theta(0)
  ,fh2CosTheta(0)
  ,fh2Jt(0)
  ,fh2PtvsZ(0)
  ,fhnIntraJet(0)
  ,fnDim(6)
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
  ,fh2Theta(copy.fh2Theta)
  ,fh2CosTheta(copy.fh2CosTheta)
  ,fh2Jt(copy.fh2Jt)
  ,fh2PtvsZ(copy.fh2PtvsZ)
  ,fhnIntraJet(copy.fhnIntraJet)
  ,fnDim(copy.fnDim)
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
    fh2Theta          = o.fh2Theta;
    fh2CosTheta       = o.fh2CosTheta;
    fh2Jt             = o.fh2Jt;
    fh2PtvsZ          = o.fh2PtvsZ;
    fhnIntraJet       = o.fhnIntraJet;
    fnDim             = o.fnDim;
    fNameIJ           = o.fNameIJ;
  }
    
  return *this;
}

//_________________________________________________________
AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::~AliFragFuncIntraJetHistos()
{
  // destructor 


  if(fh2Theta)          delete fh2Theta;
  if(fh2CosTheta)       delete fh2CosTheta;
  if(fh2Jt)             delete fh2Jt;
  if(fh2PtvsZ)          delete fh2PtvsZ;
  if(fhnIntraJet)       delete fhnIntraJet;

}

//_________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::DefineHistos()
{
  // book FF histos

  fh2Theta    = new TH2F(Form("fh2IJTheta%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsTheta, fThetaMin, fThetaMax);
  fh2CosTheta = new TH2F(Form("fh2IJcosTheta%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax,fNBinsCosTheta, fCosThetaMin, fCosThetaMax);
  fh2Jt       = new TH2F(Form("fh2IJJt%s",fNameIJ.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsJt, fJtMin, fJtMax);

  // Create 3D histograms
  Int_t    *iBin = new Int_t[fnDim];
  Double_t *min  = new Double_t[fnDim];
  Double_t *max  = new Double_t[fnDim];

  iBin[0] = fNBinsJetPt; iBin[1] = fNBinsTheta; iBin[2] = fNBinsCosTheta; iBin[3] = fNBinsJt; iBin[4] = fNBinsZ; iBin[5] = fNBinsPt;
  min[0]  = fJetPtMin; min[1] = fThetaMin; min[2] = fCosThetaMin; min[3] = fJtMin; min[4] = fZMin; min[5] = fPtMin; 
  max[0]  = fJetPtMax; max[1] = fThetaMax; max[2] = fCosThetaMax; max[3] = fJtMax; max[4] = fZMax; max[5] = fPtMax;

  const char* title = Form("fhnIntraJetPart%s",fNameIJ.Data());
  const char* comment = "THnSparseF p_{T} jet [GeV/c] : #Theta : cos(#Theta) : j_{T} : Z : p_{T} part [GeV/c]";
  fhnIntraJet = new THnSparseF(title,comment,fnDim,iBin,min,max);

  const char** axisTitle = new const char*[fnDim];
  axisTitle[0] = "p_{T}^{jet} [GeV/c]";
  axisTitle[1] = "#Theta";
  axisTitle[2] = "Cos(#Theta)";
  axisTitle[3] = "j_{T} [GeV]";
  axisTitle[4] = "z = p_{T}^{had}/p_{T}^{jet}";
  axisTitle[5] = "p_{T}^{had} [GeV/c]";
  
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Theta,"jet p_{T} [GeV/c]","#Theta","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2CosTheta,"jet p_{T} [GeV/c]","cos(#Theta)", "entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2Jt,"jet p_{T} [GeV/c]","j_{T}","entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fhnIntraJet,fnDim,axisTitle);
  delete[] iBin;
  delete[] min;
  delete[] max;
  delete[] axisTitle;
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

  // Fill histos and THnSparse
  if(norm){
    fh2CosTheta->Fill(ptJ,cosTheta,1/norm);
    fh2Theta->Fill(ptJ,theta,1/norm);
    fh2Jt->Fill(ptJ,jt,1/norm);
  }
  else {
    fh2CosTheta->Fill(ptJ,cosTheta);
    fh2Theta->Fill(ptJ,theta);
    fh2Jt->Fill(ptJ,jt);
  }

  // Fill THnSparse
  Double_t *content = new Double_t[fnDim];
  content[0]= ptJ; content[1] = theta; content[2] = cosTheta; content[3] = jt; content[4] = z; content[5] = ptT; 

  if(norm)fhnIntraJet->Fill(content,1/norm);
  else fhnIntraJet->Fill(content);

  delete [] content;

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
      if(z!=0) xi = TMath::Log(1/z);
      
      fh2Xi->Fill(jetBin, xi);
      fh2Z->Fill(jetBin, z);
    }
  if(jetType == 1)
    {
      if(incrementJetPt) fh1Jet1Pt->Fill(jetPt);
      
      fh2TrackPtJet1->Fill(jetBin, trackPt);
      
      Double_t z = trackPt / jetPt;
      Double_t xi = 0;
      if(z!=0) xi = TMath::Log(1/z);
      
      fh2Xi1->Fill(jetBin, xi);
      fh2Z1->Fill(jetBin, z);
    }
  if(jetType == 2)
    {
      if(incrementJetPt) fh1Jet2Pt->Fill(jetPt);
      
      fh2TrackPtJet2->Fill(jetBin, trackPt);
      
      Double_t z = trackPt / jetPt;
      Double_t xi = 0;
      if(z!=0) xi = TMath::Log(1/z);
      
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
                                          Int_t nDeltaPt, Float_t deltaPtMin, Float_t deltaPtMax)
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
  ,fh2InvMass(0)
  ,fh2DeltaPhi(0)
  ,fh2DeltaEta(0)
  ,fh2DeltaPt(0)
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
  ,fh2InvMass(copy.fh2InvMass)
  ,fh2DeltaPhi(copy.fh2DeltaPhi)
  ,fh2DeltaEta(copy.fh2DeltaEta)
  ,fh2DeltaPt(copy.fh2DeltaPt)
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
    fh2InvMass        = o.fh2InvMass;
    fh2DeltaPhi       = o.fh2DeltaPhi;
    fh2DeltaEta       = o.fh2DeltaEta;
    fh2DeltaPt        = o.fh2DeltaPt;
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
  
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2InvMass, xaxis, "Invariant Mass", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaPhi, xaxis, "#Delta #phi", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaEta, xaxis, "#Delta #eta", "Entries");
  AliAnalysisTaskFragmentationFunction::SetProperties(fh2DeltaPt, xaxis, "#Delta p_{T}", "Entries");

}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::FillDiJetQA(Double_t invMass, Double_t deltaPhi, Double_t deltaEta,Double_t deltaPt, Double_t jetBin)
{
  // fill dijet QA

  fh2InvMass->Fill(jetBin, invMass);
  fh2DeltaPhi->Fill(jetBin, deltaPhi);
  fh2DeltaEta->Fill(jetBin, deltaEta);
  fh2DeltaPt->Fill(jetBin, deltaPt);
}

//________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncQADiJetHistos::AddToOutput(TList* list)const
{
  // add histos to list

  list->Add(fh2InvMass);
  list->Add(fh2DeltaPhi);
  list->Add(fh2DeltaEta);
  list->Add(fh2DeltaPt);
}

//_________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::AliFragFuncIntraJetHistos::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh2CosTheta);
  list->Add(fh2Theta);
  list->Add(fh2Jt);

  list->Add(fhnIntraJet);

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

  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);

  fJetsGen = new TList();
  fJetsGen->SetOwner(kFALSE);

  fJetsRecEff = new TList();
  fJetsRecEff->SetOwner(kFALSE);

  // fJetsKine = new TList();
  // fJetsKine->SetOwner(kTRUE); // delete AOD jets using mom from Kine Tree via TList::Clear()


  //
  // Create histograms / output container
  //

  OpenFile(1);
  fCommonHistList = new TList();
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  
  // Histograms	
  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 11,-.5, 10.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1EvtMult 	             = new TH1F("fh1EvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",120,0.,120.);

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

  fh1OutLeadingMult          = new TH1F("fh1OutLeadingMult","Background multiplicity - Cone outside leading jet",120,0,120.);
  fh1PerpMult                = new TH1F("fh1PerpMult","Background multiplicity - Cone perpendicular to leading jet axis",120,0.,120.);
  fh1Out2JetsMult            = new TH1F("fh1Out2JetsMult","Background multiplicity - Cone outside 2 jets",120,0.,120.);
  fh1Out3JetsMult            = new TH1F("fh1Out3JetsMult","Background multiplicity - Cone outside 3 jets",120,0.,120.);


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
  fQAJetHistosRecEffLeading  = new AliFragFuncQAJetHistos("RecEffLeading", fQAJetNBinsPt, fQAJetPtMin, fQAJetPtMax, 
							  fQAJetNBinsEta, fQAJetEtaMin, fQAJetEtaMax,fQAJetNBinsPhi, fQAJetPhiMin, fQAJetPhiMax);


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

  fQADiJetHistosRecCuts = new AliFragFuncQADiJetHistos("RecCuts", fDiJetKindBins, 
						       fQADiJetNBinsInvMass, fQADiJetInvMassMin, fQADiJetInvMassMax,
						       fQADiJetNBinsJetPt, fQADiJetJetPtMin, fQADiJetJetPtMax,
						       fQADiJetNBinsDeltaPhi, fQADiJetDeltaPhiMin, fQADiJetDeltaPhiMax , 
						       fQADiJetNBinsDeltaEta, fQADiJetDeltaEtaMin, fQADiJetDeltaEtaMax , 
						       fQADiJetNBinsDeltaPt, fQADiJetDeltaPtMin, fQADiJetDeltaPtMax);
  fQADiJetHistosGen     = new AliFragFuncQADiJetHistos("Gen", fDiJetKindBins, 
						       fQADiJetNBinsInvMass, fQADiJetInvMassMin,  fQADiJetInvMassMax, 
						       fDiJetNBinsJetPt, fDiJetJetPtMin, fDiJetJetPtMax,
						       fQADiJetNBinsDeltaPhi, fQADiJetDeltaPhiMin, fQADiJetDeltaPhiMax,
						       fQADiJetNBinsDeltaEta, fQADiJetDeltaEtaMin, fQADiJetDeltaEtaMax,
						       fQADiJetNBinsDeltaPt, fQADiJetDeltaPtMin, fQADiJetDeltaPtMax);

  // efficiency

  fQATrackHistosRecEffGen = new AliFragFuncQATrackHistos("RecEffGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
						      fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
						      fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
						      fQATrackHighPtThreshold);

  fQATrackHistosRecEffRec = new AliFragFuncQATrackHistos("RecEffRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
						      fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
						      fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
						      fQATrackHighPtThreshold);
  
  fFFHistosRecEffGen      = new AliFragFuncHistos("RecEffGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
					       fFFNBinsPt, fFFPtMin, fFFPtMax, 
					       fFFNBinsXi, fFFXiMin, fFFXiMax,  
					       fFFNBinsZ , fFFZMin , fFFZMax);
  
  fFFHistosRecEffRec      = new AliFragFuncHistos("RecEffRec", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
					       fFFNBinsPt, fFFPtMin, fFFPtMax, 
					       fFFNBinsXi, fFFXiMin, fFFXiMax,  
					       fFFNBinsZ , fFFZMin , fFFZMax);

  // Background
  if(fBckgMode){
    fQABckgNoJetTrackHistosRec          = new AliFragFuncQATrackHistos("BckgNoJetRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    fQABckgNoJetTrackHistosRecCuts      = new AliFragFuncQATrackHistos("BckgNoJetRecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    fQABckgNoJetTrackHistosGen          = new AliFragFuncQATrackHistos("BckgNoJetGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    
    fQABckgLeadingTrackHistosRec          = new AliFragFuncQATrackHistos("BckgLeadingRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
									 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
									 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
									 fQATrackHighPtThreshold);
    fQABckgLeadingTrackHistosRecCuts      = new AliFragFuncQATrackHistos("BckgLeadingRecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
									 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
									 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
									 fQATrackHighPtThreshold);
    fQABckgLeadingTrackHistosGen          = new AliFragFuncQATrackHistos("BckgLeadingGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
									 fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
									 fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
									 fQATrackHighPtThreshold);
    
    fQABckg2JetsTrackHistosRec          = new AliFragFuncQATrackHistos("Bckg2JetsRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    fQABckg2JetsTrackHistosRecCuts      = new AliFragFuncQATrackHistos("Bckg2JetsRecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							    fQATrackHighPtThreshold);
    fQABckg2JetsTrackHistosGen          = new AliFragFuncQATrackHistos("Bckg2JetsGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    
    fQABckg3JetsTrackHistosRec          = new AliFragFuncQATrackHistos("Bckg3JetsRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    fQABckg3JetsTrackHistosRecCuts      = new AliFragFuncQATrackHistos("Bckg3JetsRecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							    fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    fQABckg3JetsTrackHistosGen          = new AliFragFuncQATrackHistos("Bckg3JetsGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								       fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								       fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								       fQATrackHighPtThreshold);
    
    
    fQABckgPerpTrackHistosRec          = new AliFragFuncQATrackHistos("BckgPerpRec", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								      fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							    fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								      fQATrackHighPtThreshold);
    fQABckgPerpTrackHistosRecCuts      = new AliFragFuncQATrackHistos("BckgPerpRecCuts", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								      fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								      fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								      fQATrackHighPtThreshold);
    fQABckgPerpTrackHistosGen          = new AliFragFuncQATrackHistos("BckgPerpGen", fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
								      fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
								      fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
								      fQATrackHighPtThreshold);

    // event with no jets
    fFFBckgNoJetHistosRecCuts    = new AliFragFuncHistos("NoJetBckgRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgNoJetHistosRecLeading = new AliFragFuncHistos("NoJetBckgRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgNoJetHistosGen   	 = new AliFragFuncHistos("NoJetBckgGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgNoJetHistosGenLeading = new AliFragFuncHistos("NoJetBckgGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    
    fIJBckgNoJetHistosRecCuts   	     = new AliFragFuncIntraJetHistos("NoJetBckgRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgNoJetHistosRecLeading        = new AliFragFuncIntraJetHistos("NoJetBckgRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgNoJetHistosGen   	     = new AliFragFuncIntraJetHistos("NoJetBckgGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgNoJetHistosGenLeading        = new AliFragFuncIntraJetHistos("NoJetBckgGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);

    
    // outside leading jet or 2 jets or more
    fFFBckgHistosRecCuts    = new AliFragFuncHistos("BckgRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						    fFFNBinsPt, fFFPtMin, fFFPtMax, 
						    fFFNBinsXi, fFFXiMin, fFFXiMax,  
						    fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosRecLeading = new AliFragFuncHistos("BckgRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						    fFFNBinsPt, fFFPtMin, fFFPtMax, 
						    fFFNBinsXi, fFFXiMin, fFFXiMax,  
						    fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosGen   	 = new AliFragFuncHistos("BckgGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						 fFFNBinsPt, fFFPtMin, fFFPtMax, 
						 fFFNBinsXi, fFFXiMin, fFFXiMax,  
						 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosGenLeading = new AliFragFuncHistos("BckgGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						    fFFNBinsPt, fFFPtMin, fFFPtMax, 
						    fFFNBinsXi, fFFXiMin, fFFXiMax,  
						    fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside leading jet
    fFFBckgLeadingHistosRecCuts    = new AliFragFuncHistos("BckgLeadingRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							   fFFNBinsPt, fFFPtMin, fFFPtMax, 
							   fFFNBinsXi, fFFXiMin, fFFXiMax,  
							   fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosRecLeading = new AliFragFuncHistos("BckgLeadingRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							   fFFNBinsPt, fFFPtMin, fFFPtMax, 
							   fFFNBinsXi, fFFXiMin, fFFXiMax,  
							   fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosGen   	 = new AliFragFuncHistos("BckgLeadingGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosGenLeading = new AliFragFuncHistos("BckgLeadingGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							   fFFNBinsPt, fFFPtMin, fFFPtMax, 
							   fFFNBinsXi, fFFXiMin, fFFXiMax,  
							   fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside leading jet or 2 jets or more
    fIJBckgHistosRecCuts   	     = new AliFragFuncIntraJetHistos("BckgRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosRecLeading        = new AliFragFuncIntraJetHistos("BckgRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								   fIJNBinsPt, fIJPtMin, fIJPtMax, 
								   fIJNBinsZ, fIJZMin, fIJZMax,  
								   fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								   fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								   fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosGen   	     = new AliFragFuncIntraJetHistos("BckgGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
							     fIJNBinsPt, fIJPtMin, fIJPtMax, 
							     fIJNBinsZ, fIJZMin, fIJZMax,  
							     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
							     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
							     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosGenLeading        = new AliFragFuncIntraJetHistos("BckgGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								   fIJNBinsPt, fIJPtMin, fIJPtMax, 
								   fIJNBinsZ, fIJZMin, fIJZMax,  
								   fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								   fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								   fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    // outside leading jet
    fIJBckgLeadingHistosRecCuts   	     = new AliFragFuncIntraJetHistos("BckgLeadingRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosRecLeading        = new AliFragFuncIntraJetHistos("BckgLeadingRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									  fIJNBinsPt, fIJPtMin, fIJPtMax, 
									  fIJNBinsZ, fIJZMin, fIJZMax,  
									  fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									  fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									  fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosGen   	     = new AliFragFuncIntraJetHistos("BckgLeadingGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosGenLeading        = new AliFragFuncIntraJetHistos("BckgLeadingGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									  fIJNBinsPt, fIJPtMin, fIJPtMax, 
									  fIJNBinsZ, fIJZMin, fIJZMax,  
									  fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									  fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									  fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    // outside 2 jets
    fFFBckg2JetsHistosRecCuts    = new AliFragFuncHistos("Bckg2JetsRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosRecLeading = new AliFragFuncHistos("Bckg2JetsRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosGen   	 = new AliFragFuncHistos("Bckg2JetsGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosGenLeading = new AliFragFuncHistos("Bckg2JetsGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside 2 jets
    fIJBckg2JetsHistosRecCuts   	     = new AliFragFuncIntraJetHistos("Bckg2JetsRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosRecLeading        = new AliFragFuncIntraJetHistos("Bckg2JetsRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosGen   	     = new AliFragFuncIntraJetHistos("Bckg2JetsGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosGenLeading        = new AliFragFuncIntraJetHistos("Bckg2JetsGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    // outside 3 jets
    fFFBckg3JetsHistosRecCuts    = new AliFragFuncHistos("Bckg3JetsRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosRecLeading = new AliFragFuncHistos("Bckg3JetsRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosGen   	 = new AliFragFuncHistos("Bckg3JetsGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosGenLeading = new AliFragFuncHistos("Bckg3JetsGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);

    // outside 3 jets
    fIJBckg3JetsHistosRecCuts   	     = new AliFragFuncIntraJetHistos("Bckg3JetsRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosRecLeading        = new AliFragFuncIntraJetHistos("Bckg3JetsRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosGen   	     = new AliFragFuncIntraJetHistos("Bckg3JetsGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosGenLeading        = new AliFragFuncIntraJetHistos("Bckg3JetsGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									fIJNBinsPt, fIJPtMin, fIJPtMax, 
									fIJNBinsZ, fIJZMin, fIJZMax,  
									fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    
    // in the transverse jet axis direction
    fFFBckgPerpHistosRecCuts    = new AliFragFuncHistos("BckgPerpRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgPerpHistosRecLeading = new AliFragFuncHistos("BckgPerpRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgPerpHistosGen   	 = new AliFragFuncHistos("BckgPerpGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgPerpHistosGenLeading = new AliFragFuncHistos("BckgPerpGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);

    // In the transverse jet axis direction 
    fIJBckgPerpHistosRecCuts   	     = new AliFragFuncIntraJetHistos("BckgPerpRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgPerpHistosRecLeading        = new AliFragFuncIntraJetHistos("BckgPerpRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								       fIJNBinsZ, fIJZMin, fIJZMax,  
								       fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								       fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								       fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgPerpHistosGen   	     = new AliFragFuncIntraJetHistos("BckgPerpGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgPerpHistosGenLeading        = new AliFragFuncIntraJetHistos("BckgPerpGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								       fIJNBinsZ, fIJZMin, fIJZMax,  
								       fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								       fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								       fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    
    // Statistical study
    // outside leading jet or 2 jets or more
    fFFBckgHistosStatRecCuts    = new AliFragFuncHistos("StatBckgRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosStatRecLeading = new AliFragFuncHistos("StatBckgRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosStatGen   	 = new AliFragFuncHistos("StatBckgGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							 fFFNBinsPt, fFFPtMin, fFFPtMax, 
							 fFFNBinsXi, fFFXiMin, fFFXiMax,  
							 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgHistosStatGenLeading = new AliFragFuncHistos("StatBckgGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							fFFNBinsPt, fFFPtMin, fFFPtMax, 
							fFFNBinsXi, fFFXiMin, fFFXiMax,  
							fFFNBinsZ , fFFZMin , fFFZMax);
  
    // outside leading jet
    fFFBckgLeadingHistosStatRecCuts    = new AliFragFuncHistos("StatBckgLeadingRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							       fFFNBinsPt, fFFPtMin, fFFPtMax, 
							       fFFNBinsXi, fFFXiMin, fFFXiMax,  
							       fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosStatRecLeading = new AliFragFuncHistos("StatBckgLeadingRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							       fFFNBinsPt, fFFPtMin, fFFPtMax, 
							       fFFNBinsXi, fFFXiMin, fFFXiMax,  
							       fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosStatGen   	 = new AliFragFuncHistos("StatBckgLeadingGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
								 fFFNBinsPt, fFFPtMin, fFFPtMax, 
								 fFFNBinsXi, fFFXiMin, fFFXiMax,  
								 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckgLeadingHistosStatGenLeading = new AliFragFuncHistos("StatBckgLeadingGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							       fFFNBinsPt, fFFPtMin, fFFPtMax, 
							       fFFNBinsXi, fFFXiMin, fFFXiMax,  
							       fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside leading jet or 2 jets or more
    fIJBckgHistosStatRecCuts   	     = new AliFragFuncIntraJetHistos("StatBckgRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosStatRecLeading        = new AliFragFuncIntraJetHistos("StatBckgRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								       fIJNBinsZ, fIJZMin, fIJZMax,  
								       fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								       fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								       fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosStatGen   	     = new AliFragFuncIntraJetHistos("StatBckgGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								     fIJNBinsPt, fIJPtMin, fIJPtMax, 
								     fIJNBinsZ, fIJZMin, fIJZMax,  
								     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgHistosStatGenLeading        = new AliFragFuncIntraJetHistos("StatBckgGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
								       fIJNBinsPt, fIJPtMin, fIJPtMax, 
								       fIJNBinsZ, fIJZMin, fIJZMax,  
								       fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
								       fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
								       fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    // outside leading jet
    fIJBckgLeadingHistosStatRecCuts   	     = new AliFragFuncIntraJetHistos("StatBckgLeadingRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosStatRecLeading        = new AliFragFuncIntraJetHistos("StatBckgLeadingRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									      fIJNBinsPt, fIJPtMin, fIJPtMax, 
									      fIJNBinsZ, fIJZMin, fIJZMax,  
									      fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									      fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									      fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosStatGen   	     = new AliFragFuncIntraJetHistos("StatBckgLeadingGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckgLeadingHistosStatGenLeading        = new AliFragFuncIntraJetHistos("StatBckgLeadingGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									      fIJNBinsPt, fIJPtMin, fIJPtMax, 
									      fIJNBinsZ, fIJZMin, fIJZMax,  
									      fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									      fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									      fIJNBinsJt , fIJJtMin , fIJJtMax);

    // outside 2 jets
    fFFBckg2JetsHistosStatRecCuts    = new AliFragFuncHistos("StatBckg2JetsRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosStatRecLeading = new AliFragFuncHistos("StatBckg2JetsRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosStatGen   	 = new AliFragFuncHistos("StatBckg2JetsGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
								 fFFNBinsPt, fFFPtMin, fFFPtMax, 
								 fFFNBinsXi, fFFXiMin, fFFXiMax,  
								 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg2JetsHistosStatGenLeading = new AliFragFuncHistos("StatBckg2JetsGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside 2 jets
    fIJBckg2JetsHistosStatRecCuts   	     = new AliFragFuncIntraJetHistos("StatBckg2JetsRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosStatRecLeading        = new AliFragFuncIntraJetHistos("StatBckg2JetsRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									    fIJNBinsPt, fIJPtMin, fIJPtMax, 
									    fIJNBinsZ, fIJZMin, fIJZMax,  
									    fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									    fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									    fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosStatGen   	     = new AliFragFuncIntraJetHistos("StatBckg2JetsGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg2JetsHistosStatGenLeading        = new AliFragFuncIntraJetHistos("StatBckg2JetsGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									    fIJNBinsPt, fIJPtMin, fIJPtMax, 
									    fIJNBinsZ, fIJZMin, fIJZMax,  
									    fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									    fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									    fIJNBinsJt , fIJJtMin , fIJJtMax);
    
    // outside 3 jets
    fFFBckg3JetsHistosStatRecCuts    = new AliFragFuncHistos("StatBckg3JetsRecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosStatRecLeading = new AliFragFuncHistos("StatBckg3JetsRecLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosStatGen   	 = new AliFragFuncHistos("StatBckg3JetsGen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
								 fFFNBinsPt, fFFPtMin, fFFPtMax, 
								 fFFNBinsXi, fFFXiMin, fFFXiMax,  
								 fFFNBinsZ , fFFZMin , fFFZMax);
    fFFBckg3JetsHistosStatGenLeading = new AliFragFuncHistos("StatBckg3JetsGenLeading", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
							     fFFNBinsPt, fFFPtMin, fFFPtMax, 
							     fFFNBinsXi, fFFXiMin, fFFXiMax,  
							     fFFNBinsZ , fFFZMin , fFFZMax);
    
    // outside 3 jets
    fIJBckg3JetsHistosStatRecCuts   	     = new AliFragFuncIntraJetHistos("StatBckg3JetsRecCuts", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax,
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosStatRecLeading        = new AliFragFuncIntraJetHistos("StatBckg3JetsRecLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									    fIJNBinsPt, fIJPtMin, fIJPtMax, 
									    fIJNBinsZ, fIJZMin, fIJZMax,  
									    fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									    fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									    fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosStatGen   	     = new AliFragFuncIntraJetHistos("StatBckg3JetsGen", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									     fIJNBinsPt, fIJPtMin, fIJPtMax, 
									     fIJNBinsZ, fIJZMin, fIJZMax,  
									     fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									     fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									     fIJNBinsJt , fIJJtMin , fIJJtMax);
    fIJBckg3JetsHistosStatGenLeading        = new AliFragFuncIntraJetHistos("StatBckg3JetsGenLeading", fIJNBinsJetPt, fIJJetPtMin, fIJJetPtMax, 
									    fIJNBinsPt, fIJPtMin, fIJPtMax, 
									    fIJNBinsZ, fIJZMin, fIJZMax,  
									    fIJNBinsCosTheta , fIJCosThetaMin , fIJCosThetaMax, 
									    fIJNBinsTheta , fIJThetaMin , fIJThetaMax, 
									    fIJNBinsJt , fIJJtMin , fIJJtMax);
  }

  fQATrackHistosRec->DefineHistos();
  fQATrackHistosRecCuts->DefineHistos();
  fQATrackHistosGen->DefineHistos();

  fQAJetHistosRec->DefineHistos();
  fQAJetHistosRecCuts->DefineHistos();
  fQAJetHistosRecCutsLeading->DefineHistos();
  fQAJetHistosGen->DefineHistos();
  fQAJetHistosGenLeading->DefineHistos();
  fQAJetHistosRecEffLeading->DefineHistos();

  fFFHistosRecCuts->DefineHistos();
  fFFHistosRecLeading->DefineHistos();
  fFFHistosRecLeadingTrack->DefineHistos();
  fFFHistosGen->DefineHistos();
  fFFHistosGenLeading->DefineHistos();
  fFFHistosGenLeadingTrack->DefineHistos();

  fIJHistosRecCuts->DefineHistos();
  fIJHistosRecLeading->DefineHistos();
  fIJHistosRecLeadingTrack->DefineHistos();
  fIJHistosGen->DefineHistos();
  fIJHistosGenLeading->DefineHistos();
  fIJHistosGenLeadingTrack->DefineHistos();
  
  fFFDiJetHistosRecCuts->DefineDiJetHistos();
  fFFDiJetHistosRecLeading->DefineDiJetHistos();
  fFFDiJetHistosRecLeadingTrack->DefineDiJetHistos();
  fFFDiJetHistosGen->DefineDiJetHistos();
  fFFDiJetHistosGenLeading->DefineDiJetHistos();
  fFFDiJetHistosGenLeadingTrack->DefineDiJetHistos();
  fQADiJetHistosRecCuts->DefineQADiJetHistos();
  fQADiJetHistosGen->DefineQADiJetHistos();

  fQATrackHistosRecEffGen->DefineHistos();
  fQATrackHistosRecEffRec->DefineHistos(); 
  fFFHistosRecEffGen->DefineHistos();
  fFFHistosRecEffRec->DefineHistos();

  // Background
  if(fBckgMode){
    fFFBckgNoJetHistosRecCuts->DefineHistos();
    fFFBckgNoJetHistosRecLeading->DefineHistos();
    fFFBckgNoJetHistosGen->DefineHistos();
    fFFBckgNoJetHistosGenLeading->DefineHistos();
    
    fFFBckgHistosRecCuts->DefineHistos();
    fFFBckgHistosRecLeading->DefineHistos();
    fFFBckgHistosGen->DefineHistos();
    fFFBckgHistosGenLeading->DefineHistos();
    
    fFFBckgLeadingHistosRecCuts->DefineHistos();
    fFFBckgLeadingHistosRecLeading->DefineHistos();
    fFFBckgLeadingHistosGen->DefineHistos();
    fFFBckgLeadingHistosGenLeading->DefineHistos();
    
    fIJBckgNoJetHistosRecCuts->DefineHistos();
    fIJBckgNoJetHistosRecLeading->DefineHistos();
    fIJBckgNoJetHistosGen->DefineHistos();
    fIJBckgNoJetHistosGenLeading->DefineHistos();
    
    fIJBckgHistosRecCuts->DefineHistos();
    fIJBckgHistosRecLeading->DefineHistos();
    fIJBckgHistosGen->DefineHistos();
    fIJBckgHistosGenLeading->DefineHistos();
    
    fIJBckgLeadingHistosRecCuts->DefineHistos();
    fIJBckgLeadingHistosRecLeading->DefineHistos();
    fIJBckgLeadingHistosGen->DefineHistos();
    fIJBckgLeadingHistosGenLeading->DefineHistos();
    
    fFFBckg2JetsHistosRecCuts->DefineHistos();
    fFFBckg2JetsHistosRecLeading->DefineHistos();
    fFFBckg2JetsHistosGen->DefineHistos();
    fFFBckg2JetsHistosGenLeading->DefineHistos();
    
    fIJBckg2JetsHistosRecCuts->DefineHistos();
    fIJBckg2JetsHistosRecLeading->DefineHistos();
    fIJBckg2JetsHistosGen->DefineHistos();
    fIJBckg2JetsHistosGenLeading->DefineHistos();
    
    fFFBckg3JetsHistosRecCuts->DefineHistos();
    fFFBckg3JetsHistosRecLeading->DefineHistos();
    fFFBckg3JetsHistosGen->DefineHistos();
    fFFBckg3JetsHistosGenLeading->DefineHistos();
    
    fIJBckg3JetsHistosRecCuts->DefineHistos();
    fIJBckg3JetsHistosRecLeading->DefineHistos();
    fIJBckg3JetsHistosGen->DefineHistos();
    fIJBckg3JetsHistosGenLeading->DefineHistos();
    
    fFFBckgPerpHistosRecCuts->DefineHistos();
    fFFBckgPerpHistosRecLeading->DefineHistos();
    fFFBckgPerpHistosGen->DefineHistos();
    fFFBckgPerpHistosGenLeading->DefineHistos();
    
    fIJBckgPerpHistosRecCuts->DefineHistos();
    fIJBckgPerpHistosRecLeading->DefineHistos();
    fIJBckgPerpHistosGen->DefineHistos();
    fIJBckgPerpHistosGenLeading->DefineHistos();
    
    fFFBckgHistosStatRecCuts->DefineHistos();
    fFFBckgHistosStatRecLeading->DefineHistos();
    fFFBckgHistosStatGen->DefineHistos();
    fFFBckgHistosStatGenLeading->DefineHistos();
    
    fFFBckgLeadingHistosStatRecCuts->DefineHistos();
    fFFBckgLeadingHistosStatRecLeading->DefineHistos();
    fFFBckgLeadingHistosStatGen->DefineHistos();
    fFFBckgLeadingHistosStatGenLeading->DefineHistos();
    
    fIJBckgHistosStatRecCuts->DefineHistos();
    fIJBckgHistosStatRecLeading->DefineHistos();
    fIJBckgHistosStatGen->DefineHistos();
    fIJBckgHistosStatGenLeading->DefineHistos();
    
    fIJBckgLeadingHistosStatRecCuts->DefineHistos();
    fIJBckgLeadingHistosStatRecLeading->DefineHistos();
    fIJBckgLeadingHistosStatGen->DefineHistos();
    fIJBckgLeadingHistosStatGenLeading->DefineHistos();
    
    fFFBckg2JetsHistosStatRecCuts->DefineHistos();
    fFFBckg2JetsHistosStatRecLeading->DefineHistos();
    fFFBckg2JetsHistosStatGen->DefineHistos();
    fFFBckg2JetsHistosStatGenLeading->DefineHistos();

    fIJBckg2JetsHistosStatRecCuts->DefineHistos();
    fIJBckg2JetsHistosStatRecLeading->DefineHistos();
    fIJBckg2JetsHistosStatGen->DefineHistos();
    fIJBckg2JetsHistosStatGenLeading->DefineHistos();
    
    fFFBckg3JetsHistosStatRecCuts->DefineHistos();
    fFFBckg3JetsHistosStatRecLeading->DefineHistos();
    fFFBckg3JetsHistosStatGen->DefineHistos();
    fFFBckg3JetsHistosStatGenLeading->DefineHistos();
    
    fIJBckg3JetsHistosStatRecCuts->DefineHistos();
    fIJBckg3JetsHistosStatRecLeading->DefineHistos();
    fIJBckg3JetsHistosStatGen->DefineHistos();
    fIJBckg3JetsHistosStatGenLeading->DefineHistos();
    
    fQABckgNoJetTrackHistosRec->DefineHistos();
    fQABckgNoJetTrackHistosRecCuts->DefineHistos();
    fQABckgNoJetTrackHistosGen->DefineHistos();
    fQABckgLeadingTrackHistosRec->DefineHistos();
    fQABckgLeadingTrackHistosRecCuts->DefineHistos();
    fQABckgLeadingTrackHistosGen->DefineHistos();
    fQABckg2JetsTrackHistosRec->DefineHistos();
    fQABckg2JetsTrackHistosRecCuts->DefineHistos();
    fQABckg2JetsTrackHistosGen->DefineHistos();
    fQABckg3JetsTrackHistosRec->DefineHistos();
    fQABckg3JetsTrackHistosRecCuts->DefineHistos();
    fQABckg3JetsTrackHistosGen->DefineHistos();
    fQABckgPerpTrackHistosRec->DefineHistos();
    fQABckgPerpTrackHistosRecCuts->DefineHistos();
    fQABckgPerpTrackHistosGen->DefineHistos();
  }
  
  Bool_t genJets    = (fJetTypeGen != kJetsUndef) ? kTRUE : kFALSE;
  Bool_t genTracks  = (fTrackTypeGen != kTrackUndef) ? kTRUE : kFALSE;
  Bool_t recJetsEff = (fJetTypeRecEff != kJetsUndef) ? kTRUE : kFALSE;


  const Int_t saveLevel = 5;
  if(saveLevel>0){
    fCommonHistList->Add(fh1EvtSelection);
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecLeading->AddToOutput(fCommonHistList);
    fFFHistosRecLeadingTrack->AddToOutput(fCommonHistList);

    // Background
    if(fBckgMode){
      fFFBckgNoJetHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckgNoJetHistosRecLeading->AddToOutput(fCommonHistList);
      fFFBckgHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckgHistosRecLeading->AddToOutput(fCommonHistList);
      fFFBckgLeadingHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckgLeadingHistosRecLeading->AddToOutput(fCommonHistList);
      fFFBckg2JetsHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckg2JetsHistosRecLeading->AddToOutput(fCommonHistList);
      fFFBckg3JetsHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckg3JetsHistosRecLeading->AddToOutput(fCommonHistList);
      fFFBckgPerpHistosRecCuts->AddToOutput(fCommonHistList);
      fFFBckgPerpHistosRecLeading->AddToOutput(fCommonHistList);
      
      fFFBckgHistosStatRecCuts->AddToOutput(fCommonHistList);
      fFFBckgHistosStatRecLeading->AddToOutput(fCommonHistList);
      fFFBckgLeadingHistosStatRecCuts->AddToOutput(fCommonHistList);
      fFFBckgLeadingHistosStatRecLeading->AddToOutput(fCommonHistList);
      fFFBckg2JetsHistosStatRecCuts->AddToOutput(fCommonHistList);
      fFFBckg2JetsHistosStatRecLeading->AddToOutput(fCommonHistList);
      fFFBckg3JetsHistosStatRecCuts->AddToOutput(fCommonHistList);
      fFFBckg3JetsHistosStatRecLeading->AddToOutput(fCommonHistList);
      
      fQABckgNoJetTrackHistosRec->AddToOutput(fCommonHistList);
      fQABckgNoJetTrackHistosRecCuts->AddToOutput(fCommonHistList);
      fQABckgLeadingTrackHistosRec->AddToOutput(fCommonHistList);
      fQABckgLeadingTrackHistosRecCuts->AddToOutput(fCommonHistList);
      fQABckg2JetsTrackHistosRec->AddToOutput(fCommonHistList);
      fQABckg2JetsTrackHistosRecCuts->AddToOutput(fCommonHistList);
      fQABckg3JetsTrackHistosRec->AddToOutput(fCommonHistList);
      fQABckg3JetsTrackHistosRecCuts->AddToOutput(fCommonHistList);
      fQABckgPerpTrackHistosRec->AddToOutput(fCommonHistList);
      fQABckgPerpTrackHistosRecCuts->AddToOutput(fCommonHistList);
    }

    fCommonHistList->Add(fh1OutLeadingMult);
    fCommonHistList->Add(fh1PerpMult);
    fCommonHistList->Add(fh1Out2JetsMult);
    fCommonHistList->Add(fh1Out3JetsMult);


    if(genJets && genTracks){
       fFFHistosGen->AddToOutput(fCommonHistList);
       fFFHistosGenLeading->AddToOutput(fCommonHistList);
       fFFHistosGenLeadingTrack->AddToOutput(fCommonHistList);

       fCommonHistList->Add(fh1Xsec);
       fCommonHistList->Add(fh1Trials);
       fCommonHistList->Add(fh1PtHard);
       fCommonHistList->Add(fh1PtHardTrials);

       // Background
       if(fBckgMode){
	 fFFBckgNoJetHistosGen->AddToOutput(fCommonHistList);
	 fFFBckgNoJetHistosGenLeading->AddToOutput(fCommonHistList);
	 fFFBckgHistosGen->AddToOutput(fCommonHistList);
	 fFFBckgHistosGenLeading->AddToOutput(fCommonHistList);
	 fFFBckgLeadingHistosGen->AddToOutput(fCommonHistList);
	 fFFBckgLeadingHistosGenLeading->AddToOutput(fCommonHistList);
	 fFFBckg2JetsHistosGen->AddToOutput(fCommonHistList);
	 fFFBckg2JetsHistosGenLeading->AddToOutput(fCommonHistList);
	 fFFBckg3JetsHistosGen->AddToOutput(fCommonHistList);
	 fFFBckg3JetsHistosGenLeading->AddToOutput(fCommonHistList);
	 fFFBckgPerpHistosGen->AddToOutput(fCommonHistList);
	 fFFBckgPerpHistosGenLeading->AddToOutput(fCommonHistList);
	 
	 fFFBckgHistosStatGen->AddToOutput(fCommonHistList);
	 fFFBckgHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fFFBckgLeadingHistosStatGen->AddToOutput(fCommonHistList);
	 fFFBckgLeadingHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fFFBckg2JetsHistosStatGen->AddToOutput(fCommonHistList);
	 fFFBckg2JetsHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fFFBckg3JetsHistosStatGen->AddToOutput(fCommonHistList);
	 fFFBckg3JetsHistosStatGenLeading->AddToOutput(fCommonHistList);
	 
	 fQABckgNoJetTrackHistosGen->AddToOutput(fCommonHistList);
	 fQABckgLeadingTrackHistosGen->AddToOutput(fCommonHistList);
	 fQABckg2JetsTrackHistosGen->AddToOutput(fCommonHistList);
	 fQABckg3JetsTrackHistosGen->AddToOutput(fCommonHistList);
	 fQABckgPerpTrackHistosGen->AddToOutput(fCommonHistList);
       }
    }
  }
  if(saveLevel>1){
    fQATrackHistosRec->AddToOutput(fCommonHistList);
    fQATrackHistosRecCuts->AddToOutput(fCommonHistList);
    if(genTracks) fQATrackHistosGen->AddToOutput(fCommonHistList);
    
    fQAJetHistosRec->AddToOutput(fCommonHistList);
    fQAJetHistosRecCuts->AddToOutput(fCommonHistList);
    fQAJetHistosRecCutsLeading->AddToOutput(fCommonHistList);
    if(recJetsEff) fQAJetHistosRecEffLeading->AddToOutput(fCommonHistList); 

    if(genJets){
       fQAJetHistosGen->AddToOutput(fCommonHistList);
       fQAJetHistosGenLeading->AddToOutput(fCommonHistList);
    }

    fCommonHistList->Add(fh1EvtMult);
    fCommonHistList->Add(fh1nRecJetsCuts);
    if(genJets) fCommonHistList->Add(fh1nGenJets);
  }
  if(saveLevel>2){
    fCommonHistList->Add(fh1VertexNContributors);
    fCommonHistList->Add(fh1VertexZ);    
  }
  if(saveLevel>3){
    fIJHistosRecCuts->AddToOutput(fCommonHistList);
    fIJHistosRecLeading->AddToOutput(fCommonHistList);
    fIJHistosRecLeadingTrack->AddToOutput(fCommonHistList);

    // Background
    if(fBckgMode){
      fIJBckgNoJetHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckgNoJetHistosRecLeading->AddToOutput(fCommonHistList);
      fIJBckgHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckgHistosRecLeading->AddToOutput(fCommonHistList);
      fIJBckgLeadingHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckgLeadingHistosRecLeading->AddToOutput(fCommonHistList);
      fIJBckg2JetsHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckg2JetsHistosRecLeading->AddToOutput(fCommonHistList);
      fIJBckg3JetsHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckg3JetsHistosRecLeading->AddToOutput(fCommonHistList);
      fIJBckgPerpHistosRecCuts->AddToOutput(fCommonHistList);
      fIJBckgPerpHistosRecLeading->AddToOutput(fCommonHistList);
      
      fIJBckgHistosStatRecCuts->AddToOutput(fCommonHistList);
      fIJBckgHistosStatRecLeading->AddToOutput(fCommonHistList);
      fIJBckgLeadingHistosStatRecCuts->AddToOutput(fCommonHistList);
      fIJBckgLeadingHistosStatRecLeading->AddToOutput(fCommonHistList);
      fIJBckg2JetsHistosStatRecCuts->AddToOutput(fCommonHistList);
      fIJBckg2JetsHistosStatRecLeading->AddToOutput(fCommonHistList);
      fIJBckg3JetsHistosStatRecCuts->AddToOutput(fCommonHistList);
      fIJBckg3JetsHistosStatRecLeading->AddToOutput(fCommonHistList);
    }

    if(genJets && genTracks){
       fIJHistosGen->AddToOutput(fCommonHistList);
       fIJHistosGenLeading->AddToOutput(fCommonHistList);
       fIJHistosGenLeadingTrack->AddToOutput(fCommonHistList);

       // Background
       if(fBckgMode){
	 fIJBckgNoJetHistosGen->AddToOutput(fCommonHistList);
	 fIJBckgNoJetHistosGenLeading->AddToOutput(fCommonHistList);
	 fIJBckgHistosGen->AddToOutput(fCommonHistList);
	 fIJBckgHistosGenLeading->AddToOutput(fCommonHistList);
	 fIJBckgLeadingHistosGen->AddToOutput(fCommonHistList);
	 fIJBckgLeadingHistosGenLeading->AddToOutput(fCommonHistList);
	 fIJBckg2JetsHistosGen->AddToOutput(fCommonHistList);
	 fIJBckg2JetsHistosGenLeading->AddToOutput(fCommonHistList);
	 fIJBckg3JetsHistosGen->AddToOutput(fCommonHistList);
	 fIJBckg3JetsHistosGenLeading->AddToOutput(fCommonHistList);
	 fIJBckgPerpHistosGen->AddToOutput(fCommonHistList);
	 fIJBckgPerpHistosGenLeading->AddToOutput(fCommonHistList);
	 
	 fIJBckgHistosStatGen->AddToOutput(fCommonHistList);
	 fIJBckgHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fIJBckgLeadingHistosStatGen->AddToOutput(fCommonHistList);
	 fIJBckgLeadingHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fIJBckg2JetsHistosStatGen->AddToOutput(fCommonHistList);
	 fIJBckg2JetsHistosStatGenLeading->AddToOutput(fCommonHistList);
	 fIJBckg3JetsHistosStatGen->AddToOutput(fCommonHistList);
	 fIJBckg3JetsHistosStatGenLeading->AddToOutput(fCommonHistList);
       }
    }
  }
  if(saveLevel>4){
    fFFDiJetHistosRecCuts->AddToOutput(fCommonHistList);
    fFFDiJetHistosRecLeading->AddToOutput(fCommonHistList);
    fFFDiJetHistosRecLeadingTrack->AddToOutput(fCommonHistList);
    fQADiJetHistosRecCuts->AddToOutput(fCommonHistList);

    if(genJets && genTracks){
        fFFDiJetHistosGen->AddToOutput(fCommonHistList);
        fFFDiJetHistosGenLeading->AddToOutput(fCommonHistList);
        fFFDiJetHistosGenLeadingTrack->AddToOutput(fCommonHistList);
        fQADiJetHistosGen->AddToOutput(fCommonHistList);
    }
    
    if(recJetsEff && genTracks){
      fQATrackHistosRecEffGen->AddToOutput(fCommonHistList);
      fQATrackHistosRecEffRec->AddToOutput(fCommonHistList);
      fFFHistosRecEffGen->AddToOutput(fCommonHistList);
      fFFHistosRecEffRec->AddToOutput(fCommonHistList);
      fCommonHistList->Add(fh1nRecEffJets);
      fCommonHistList->Add(fh2PtRecVsGenPrim); 
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
  if(inputHandler->IsEventSelected() & AliVEvent::kMB){
    if(fDebug > 1)  Printf(" Trigger Selection: event ACCEPTED ... ");
    fh1EvtSelection->Fill(1.);
  } else {
    fh1EvtSelection->Fill(0.);
    if(inputHandler->InheritsFrom("AliESDInputHandler") && fUsePhysicsSelection){ // PhysicsSelection only with ESD input
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
  

  // event selection (vertex) *****************************************
  
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  if(!nTracksPrim){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: event REJECTED...",(char*)__FILE__,__LINE__); 
    fh1EvtSelection->Fill(2.);
    PostData(1, fCommonHistList);
    return;
  }

  fh1VertexZ->Fill(primVtx->GetZ());
  
  if(TMath::Abs(primVtx->GetZ())>10){
    if (fDebug > 1) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ()); 
    fh1EvtSelection->Fill(3.);
    PostData(1, fCommonHistList);
    return; 
  }

  TString primVtxName(primVtx->GetName());

  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(4.);
    PostData(1, fCommonHistList);
    return;
  }
  if (fDebug > 1) Printf("%s:%d primary vertex selection: event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  fh1EvtSelection->Fill(5.);


  //___ get MC information __________________________________________________________________

  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  
  if(fMCEvent){
     AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
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
  
  // loop over tracks

  for(Int_t it=0; it<nRecPart; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRec->At(it));
    fQATrackHistosRec->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  for(Int_t it=0; it<nRecPartCuts; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(it));
    fQATrackHistosRecCuts->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  for(Int_t it=0; it<nGenPart; ++it){
    AliVParticle *part = dynamic_cast<AliVParticle*>(fTracksGen->At(it));
    fQATrackHistosGen->FillTrackQA( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt());
  }
  
  // loop over jets

  for(Int_t ij=0; ij<nRecJets; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRec->At(ij));
    fQAJetHistosRec->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
  }
  
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

  if(fBckgMode && nRecJetsCuts>1) {
    GetOutNJetsTracks(2,fTracksRecCuts, tracklistout2jets, fJetsRecCuts, sumPtOut2Jets);
    GetOutNJetsTracksStat(2,fTracksRecCuts, tracklistout2jetsStat, fJetsRecCuts,sumPtOut2JetsStat, normFactor2Jets);
  }
  if(fBckgMode && nRecJetsCuts>2) {
    GetOutNJetsTracks(3,fTracksRecCuts, tracklistout3jets, fJetsRecCuts, sumPtOut3Jets);
    GetOutNJetsTracksStat(3,fTracksRecCuts, tracklistout3jetsStat, fJetsRecCuts, sumPtOut3JetsStat, normFactor3Jets);
  }

  for(Int_t ij=0; ij<nRecJetsCuts; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRecCuts->At(ij));
    fQAJetHistosRecCuts->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());

    if(ij==0){ // leading jet
      
      fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
      
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
	TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	Float_t jetPt   = jet->Pt();
	Float_t trackPt = trackV->Pt();

	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	fFFHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt);
	fIJHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	
	if(it==0){ // leading track 
	  leadTrackPt = trackPt;
	  leadTrackV->SetPxPyPzE(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  fFFHistosRecLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE);
	  fIJHistosRecLeadingTrack->FillIntraJet( leadTrackV, jet->MomentumVector() );
	}
	fFFHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt);
	fIJHistosRecLeading->FillIntraJet( trackV, leadTrackV );

	delete trackV;
      }

      // ff and ij for background study
      if(fBckgMode){
	TList* tracklistoutleading       = new TList();
	TList* tracklistoutleadingStat   = new TList();
	Double_t sumPtOutLeading     = 0.; 
	Double_t sumPtOutLeadingStat = 0.; 
	Double_t normFactorLeading   = 0.;
	GetOutNJetsTracks(1,fTracksRecCuts, tracklistoutleading, fJetsRecCuts, sumPtOutLeading);
	GetOutNJetsTracksStat(1,fTracksRecCuts, tracklistoutleadingStat, fJetsRecCuts, sumPtOutLeadingStat, normFactorLeading);
	fh1OutLeadingMult->Fill(tracklistoutleading->GetSize());

	for(Int_t it=0; it<tracklistoutleading->GetSize(); ++it){

	  AliVParticle* trackVP   = dynamic_cast<AliVParticle*>(tracklistoutleading->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  fFFBckgLeadingHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt);
	  fIJBckgLeadingHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	
	  fFFBckgLeadingHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt);
	  fIJBckgLeadingHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  
	  // Fill track QA for background
	  fQABckgLeadingTrackHistosRecCuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  
	  // All cases included
	  if(nRecJetsCuts==1){
	    fFFBckgHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	    fIJBckgHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    fFFBckgHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    fIJBckgHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  }
	  
	  delete trackV;
	}

	for(Int_t it=0; it<tracklistoutleadingStat->GetSize(); ++it){

	  AliVParticle* trackVP   = dynamic_cast<AliVParticle*>(tracklistoutleadingStat->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();

	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  // Stat plots
	  fFFBckgLeadingHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);
	  fIJBckgLeadingHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorLeading);
	  
	  fFFBckgLeadingHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorLeading);
	  fIJBckgLeadingHistosStatRecLeading->FillIntraJet( trackV, leadTrackV, normFactorLeading);
	
	  // Fill track QA for background
	  //fQABckgLeadingTrackHistosRecCuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	  
	  // All cases included
	  if(nRecJetsCuts==1){
	    fFFBckgHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt, normFactorLeading);
	    fIJBckgHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactorLeading );
	    
	    fFFBckgHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactorLeading);
	    fIJBckgHistosStatRecLeading->FillIntraJet( trackV, leadTrackV, normFactorLeading );
	  }
	  
	  delete trackV;
	}
	
	Double_t sumPtPerp = 0.;
	TList* tracklistperp = new TList();
	GetOutPerpJetTracks(fTracksRecCuts,tracklistperp,jet,GetFFRadius(),sumPtPerp);
	fh1PerpMult->Fill(tracklistperp->GetSize());
	
	for(Int_t it=0; it<tracklistperp->GetSize(); ++it){

	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistperp->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  fFFBckgPerpHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	  fIJBckgPerpHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	  
	  fFFBckgPerpHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	  fIJBckgPerpHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  
	  // Fill track QA for background
	  fQABckgPerpTrackHistosRecCuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  
	  delete trackV;
	}

	fh1Out2JetsMult->Fill(tracklistout2jets->GetSize());
	for(Int_t it=0; it<tracklistout2jets->GetSize(); ++it){

	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jets->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  fFFBckg2JetsHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	  fIJBckg2JetsHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	  
	  fFFBckg2JetsHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	  fIJBckg2JetsHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  
	  fQABckg2JetsTrackHistosRecCuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	
	  // All cases included
	  if(nRecJetsCuts==2){
	    fFFBckgHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	    fIJBckgHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    fFFBckgHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    fIJBckgHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  }
	  
	  delete trackV;
	}

	for(Int_t it=0; it<tracklistout2jetsStat->GetSize(); ++it){

	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout2jetsStat->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  fFFBckg2JetsHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	  fIJBckg2JetsHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor2Jets );
	  
	  fFFBckg2JetsHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor2Jets);
	  fIJBckg2JetsHistosStatRecLeading->FillIntraJet( trackV, leadTrackV, normFactor2Jets );
	  
	  //fQABckg2JetsTrackHistosRecCuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	
	  // All cases included
	  if(nRecJetsCuts==2){
	    fFFBckgHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor2Jets);
	    fIJBckgHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor2Jets );
	    
	    fFFBckgHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor2Jets);
	    fIJBckgHistosStatRecLeading->FillIntraJet( trackV, leadTrackV, normFactor2Jets );
	  }
	  
	  delete trackV;
	}
	
	fh1Out3JetsMult->Fill(tracklistout3jets->GetSize());
      
	for(Int_t it=0; it<tracklistout3jets->GetSize(); ++it){

	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jets->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	  fFFBckg3JetsHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	  fIJBckg3JetsHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	  
	  fFFBckg3JetsHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	  fIJBckg3JetsHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  
	  fQABckg3JetsTrackHistosRecCuts->FillTrackQA( trackV->Eta(), TVector2::Phi_0_2pi(trackV->Phi()), trackPt);
	  
	  // All cases included
	  if(nRecJetsCuts==3){
	    fFFBckgHistosRecCuts->FillFF( trackPt, jetPt, incrementJetPt );
	    fIJBckgHistosRecCuts->FillIntraJet( trackV, jet->MomentumVector() );
	    
	    fFFBckgHistosRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt );
	    fIJBckgHistosRecLeading->FillIntraJet( trackV, leadTrackV );
	  }
	  
	  delete trackV;
	  
	}
	
	for(Int_t it=0; it<tracklistout3jetsStat->GetSize(); ++it){

	  AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(tracklistout3jetsStat->At(it));
	  TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  Float_t jetPt   = jet->Pt();
	  Float_t trackPt = trackV->Pt();
	  
	  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  
	  fFFBckg3JetsHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt,normFactor3Jets);
	  fIJBckg3JetsHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor3Jets);
	  
	  fFFBckg3JetsHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt,normFactor3Jets);
	  fIJBckg3JetsHistosStatRecLeading->FillIntraJet( trackV, leadTrackV, normFactor3Jets);
	  
	  //	fQABckg3JetsTrackHistosRecCuts->FillTrackQA( trackEta, TVector2::Phi_0_2pi(trackPhi), trackPt);
	  
	  // All cases included
	  if(nRecJetsCuts==3){
	    fFFBckgHistosStatRecCuts->FillFF( trackPt, jetPt, incrementJetPt, normFactor3Jets );
	    fIJBckgHistosStatRecCuts->FillIntraJet( trackV, jet->MomentumVector(), normFactor3Jets);
	    
	    fFFBckgHistosStatRecLeading->FillFF( trackPt, leadTrackPt , incrementJetPt, normFactor3Jets );
	    fIJBckgHistosStatRecLeading->FillIntraJet( trackV, leadTrackV,normFactor3Jets );
	  }

	  delete trackV;
	}
	
	delete tracklistoutleading;
	delete tracklistoutleadingStat;
	delete tracklistperp;
      } // end if(fBckgMode)

      delete leadTrackV;
      delete jettracklist;
    }
  }

  delete tracklistout2jets;
  delete tracklistout3jets;
  delete tracklistout2jetsStat;
  delete tracklistout3jetsStat;

  // generated jets

  for(Int_t ij=0; ij<nGenJets; ++ij){

    AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsGen->At(ij));
    fQAJetHistosGen->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
    
    if(ij==0){ // leading jet

      fQAJetHistosGenLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt());
      
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
	
	fFFHistosGen->FillFF( trackPt, jetPt, incrementJetPt );
	fIJHistosGen->FillIntraJet( trackV, jet->MomentumVector() );
	
	if(it==0){ // leading track
	  leadTrackPt = trackPt;
	  leadTrackV->SetPxPyPzE(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());

	  fFFHistosGenLeadingTrack->FillFF( leadTrackPt, jetPt, kTRUE );	  
	  fIJHistosGenLeadingTrack->FillIntraJet( leadTrackV, jet->MomentumVector() );
	}
	fFFHistosGenLeading->FillFF( trackPt, leadTrackPt, incrementJetPt );
	fIJHistosGenLeading->FillIntraJet( trackV, leadTrackV );

	delete trackV;
      }
      
      delete leadTrackV;
      delete jettracklist;
    }
  }

  //_______ DiJet part _____________________________________________________

  if (nRecJetsCuts > 1) 
  {

    AliAODJet* jet1 = dynamic_cast<AliAODJet*>(fJetsRecCuts->At(0));
    AliAODJet* jet2 = dynamic_cast<AliAODJet*>(fJetsRecCuts->At(1));
    
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
	Double_t meanEt        = (Double_t)((et1+et2)/2.);
	Double_t invariantMass = (Double_t)InvMass(jet1,jet2);
	
	Double_t  jetBin = GetDiJetBin(invariantMass, jet1->Pt(), meanEt, fDiJetKindBins);

	if (jetBin > 0)
	  {
	    fQADiJetHistosRecCuts->FillDiJetQA(invariantMass, deltaPhi, deltaEta, deltaPt, jetBin);
	    
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
		    Float_t trackPt1 = (dynamic_cast<AliVParticle*> (jettracklist1->At(it)))->Pt();
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
		    Float_t trackPt2   = (dynamic_cast<AliVParticle*>(jettracklist2->At(it)))->Pt();
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

	  } // End if(jetBin > 0)
	else { Printf("Jet bins for di-jet studies not set !");}
      } // End if(go)
  } // End if(nRecJets > 1)

  if (nGenJets > 1)
  {
    AliAODJet* jet1 = dynamic_cast<AliAODJet*>(fJetsGen->At(0));
    AliAODJet* jet2 = dynamic_cast<AliAODJet*>(fJetsGen->At(1));

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
      Double_t meanEt        = (Double_t)((et1+et2)/2.);
      Double_t invariantMass = (Double_t)InvMass(jet1,jet2);

      Double_t jetBin = GetDiJetBin(invariantMass, jet1->Pt(), meanEt, fDiJetKindBins);

      if(jetBin > 0)
      {
        fQADiJetHistosGen->FillDiJetQA(invariantMass, deltaPhi, deltaEta, deltaPt, jetBin);

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
  	    Float_t trackPt1 = (dynamic_cast<AliVParticle*>(jettracklist1->At(it)))->Pt();
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

      } // End if(jetBin > 0)
      else { Printf("Jet bins for di-jet studies not set !");}
    } // End if (go)
  } // End if(nGenJets > 1)

  
  // ____ efficiency _______________________________

  if(fJetTypeRecEff != kJetsUndef){

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
    FillSingleTrackRecEffHisto(fQATrackHistosRecEffGen,fQATrackHistosRecEffRec,fTracksAODMCCharged,indexAODTr,isGenPrim);

    // jet track eff
    
    Double_t sumPtGenLeadingJetRecEff = 0;
    Double_t sumPtRecLeadingJetRecEff = 0;
    
    for(Int_t ij=0; ij<nRecEffJets; ++ij){ 
    
      AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsRecEff->At(ij));
      
      if(ij==0){ // leading jet
	
	TList* jettracklistGen = new TList();
	GetJetTracksPointing(fTracksGen, jettracklistGen, jet, GetFFRadius(), sumPtGenLeadingJetRecEff); // for efficiency: gen tracks from pointing with gen/rec jet
	
	TList* jettracklistRec = new TList();
	GetJetTracksPointing(fTracksRecCuts,jettracklistRec, jet, GetFFRadius(), sumPtRecLeadingJetRecEff); // bin efficiency in jet pt bins using rec tracks  
	
	Double_t jetEta   = jet->Eta();
	Double_t jetPhi   = TVector2::Phi_0_2pi(jet->Phi());
	
	fQAJetHistosRecEffLeading->FillJetQA( jetEta, jetPhi, sumPtGenLeadingJetRecEff ); 
	
	FillJetTrackRecEffHisto(fFFHistosRecEffGen,fFFHistosRecEffRec,sumPtGenLeadingJetRecEff,sumPtRecLeadingJetRecEff,
				jettracklistGen,fTracksAODMCCharged,indexAODTr,isGenPrim,fUseRecEffRecJetPtBins); 
      
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
	  
	  GetOutPerpJetTracks(fTracksGen, perpjettracklistGen, jet, GetFFBckgRadius() , sumPtGen); // for efficiency: gen tracks perp to gen/rec jet
	  
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
      
      GetOutNJetsTracks(nCases, fTracksGen, outjettracklistGen, fJetsRecEff, sumPtGen); // for efficiency: gen tracks outide n gen/rec jets 
      
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

  if(type==kTrackUndef) return 0;
  
  Int_t iCount = 0;
  if(type==kTrackAODCuts || type==kTrackAODQualityCuts || type==kTrackAOD){

    // all rec. tracks, esd filter mask, eta range
    if(!fAOD) return -1;
    
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

// _________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::SetProperties(THnSparse* h,const Int_t dim, const char** labels)
{
  //Set properties of THnSparse 

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
void AliAnalysisTaskFragmentationFunction::SetProperties(TH2* h,const char* x, const char* y, const char* z)
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

    Int_t pdg = gentrack->GetPdgCode();    

    // 211 - pi, 2212 - proton, 321 - Kaon, 11 - electron, 13 - muon
    if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || 
       TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){
      
      isGenPrim[iGen] = kTRUE;

      Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track 

      if(iRec>=0){
	Float_t genPt = gentrack->Pt();
	Float_t recPt = (dynamic_cast<AliAODTrack*>(tracksRec->At(iRec)))->Pt(); 
	
	fh2PtRecVsGenPrim->Fill(genPt,recPt);
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
    
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF 
    // could in principle also be done setting THNsparse axis limits before projecting, 
    // but then the binning needs to be fine grained enough 

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
								   TList* tracksGen, const TArrayI& indexAODTr, const TArrayS& isGenPrim, const Bool_t useRecJetPt)
{

  // fill objects for jet track reconstruction efficiency

  Int_t nTracksJet = jetTrackList->GetSize(); // list with AODMC tracks

  if(!nTracksJet) return; 

  Bool_t incrementJetPtGenFF = kTRUE; // needed in case we fill FFHistos
  Bool_t incrementJetPtRecFF = kTRUE; // needed in case we fill FFHistos

  for(Int_t iTr=0; iTr<nTracksJet; iTr++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (jetTrackList->At(iTr));

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
      
      AliFragFuncHistos* effFFhistGen = dynamic_cast<AliFragFuncHistos*>(histGen); 
      AliFragFuncHistos* effFFhistRec = dynamic_cast<AliFragFuncHistos*>(histRec); 

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

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetOutPerpJetTracks(TList* inputlist, TList* outputlist, AliAODJet* jet, Double_t radius,Double_t& sumPt)
{
  // List of tracks in cone perpendicular to the jet azimuthal direction

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);
  // Rotate phi and keep eta unchanged
  Double_t ptPerp = jet3mom.Pt();
  Double_t etaPerp = jet3mom.Eta();
  Double_t phiPerp = TVector2::Phi_0_2pi(jet3mom.Phi()) + TMath::Pi()/2;
  if(phiPerp > 2*TMath::Pi()) phiPerp = phiPerp - 2*TMath::Pi();
  TVector3 vPerp;
  vPerp.SetPtEtaPhi(ptPerp,etaPerp,phiPerp);

  // Take orthogonal vector to jet direction
  // TVector3 vPerp(jet3mom.Orthogonal());

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));

    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = vPerp.DeltaR(track3mom);
    if(dR<=radius){
      outputlist->Add(track);
      sumPt += track->Pt();
    }
  }

}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetOutNJetsTracks(Int_t nCases, TList* inputlist, TList* outputlist, TList* jetlist, Double_t& sumPt)
{
  // List of tracks outside cone around N jet axis  
  // Particles taken randomly

  sumPt = 0;
  //  Int_t   nj  = jetlist->GetSize();
  Float_t rc  = GetFFRadius();
  Float_t rcl = GetFFBckgRadius();

  // Estimate jet and background areas
  Float_t* areaJet = new Float_t[nCases];
  Float_t* areaJetLarge = new Float_t[nCases];
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

  Int_t nScaled = (Int_t) (nOut * areaJet[0] / areaOut + 0.5); 

  TRandom* rd = new TRandom();

  Int_t* t = new Int_t[nScaled];
  for(Int_t i = 0; i < nScaled; i++)
    t[i] = -1;

  // Take randomly nScaled particles out of nOut and fill the list
  for (Int_t it=0; it<nScaled; it++)
    {
      rd->SetSeed(0);
      Int_t temp = rd->Integer(nOut);

      Bool_t ok = 1;
      for(Int_t j=0; j<nScaled; j++)
        {
          if(temp==t[j])
            {
              ok = 0;
              break;
            }
        }

      if(ok) 
	{ 
	  t[it] = temp;
	  AliVParticle* track = dynamic_cast<AliVParticle*>(templist->At(temp));
	  outputlist->Add(track);
	  sumPt += track->Pt();
	}
  }

  outputlist->Sort();

  delete rd;
  delete [] t;
  delete vect3Jet;
  delete templist;
  delete [] areaJetLarge;
  delete [] areaJet;

}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskFragmentationFunction::GetOutNJetsTracksStat(Int_t nCases, TList* inputlist, TList* outputlist, TList* jetlist, Double_t& sumPt, Double_t &normFactor)
{
  // List of tracks outside cone around N jet axis  
  // All particles taken + final scaling factor 

  sumPt = 0;
  Float_t rc  = GetFFRadius();
  Float_t rcl = GetFFBckgRadius();

  // Estimate jet and background areas
  Float_t* areaJet = new Float_t[nCases];
  Float_t* areaJetLarge = new Float_t[nCases];
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
Float_t AliAnalysisTaskFragmentationFunction::CalcJetArea(Float_t etaJet, Float_t rc)
{

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
