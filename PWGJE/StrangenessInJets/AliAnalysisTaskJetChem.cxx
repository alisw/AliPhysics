 /*************************************************************************
 *                                                                       *
 *                                                                       *
 *      Task for Jet Chemistry Analysis in PWG4 Jet Task Force Train     *
 *    Analysis of K0s, Lambda and Antilambda with and without Jetevents  *
 *                                                                       *
 *************************************************************************/

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby grante    *
 *                                                                        *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
  

/* $Id: */

#include <iostream>
#include "TH2.h"
#include "TH3.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "THnSparse.h"
#include <algorithm>
#include <string> 
#include "AliAnalysisHelperJetTasks.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h" 
#include "AliAODInputHandler.h" 
#include "AliESDEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include "TLorentzVector.h"
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliPID.h" 
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliExternalTrackParam.h"
#include "AliAnalysisTaskJetChem.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliInputEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODPid.h"
#include "AliVEvent.h"
#include "AliAODMCParticle.h"
#include "TVector3.h"
#include "TRandom.h"

ClassImp(AliAnalysisTaskJetChem)

//____________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem()
   : AliAnalysisTaskFragmentationFunction()

   ,fAnalysisMC(0)
   ,fDeltaVertexZ(0)
   ,fCuttrackNegNcls(0)
   ,fCuttrackPosNcls(0)
   ,fCutPostrackRap(0)
   ,fCutNegtrackRap(0)
   ,fCutRap(0)
   ,fCutPostrackEta(0)
   ,fCutNegtrackEta(0)
   ,fCutEta(0)
   ,fCutV0cosPointAngle(0)
   ,fCutChi2PosDaughter(0)
   ,fCutChi2NegDaughter(0)
   ,fKinkDaughters(0)
   ,fRequireTPCRefit(0)
   ,fCutArmenteros(0)
   ,fCutV0DecayMin(0)
   ,fCutV0DecayMax(0)
   ,fCutV0totMom(0)
   ,fCutDcaV0Daughters(0)
   ,fCutDcaPosToPrimVertex(0)
   ,fCutDcaNegToPrimVertex(0)
   ,fCutV0RadiusMin(0)
   ,fCutV0RadiusMax(0)
   ,fCutBetheBloch(0)
   ,fCutRatio(0)
   ,fK0Type(0)  
   ,fFilterMaskK0(0)
   ,fListK0s(0)
   ,fPIDResponse(0)
   ,fV0QAK0(0)
   ,fFFHistosRecCutsK0Evt(0)      
   ,fFFHistosIMK0AllEvt(0)        
   ,fFFHistosIMK0Jet(0)           
   ,fFFHistosIMK0Cone(0)
   ,fFFHistosPhiCorrIMK0(0)
   ,fLaType(0) 
   ,fFilterMaskLa(0)
   ,fListLa(0)
   ,fFFHistosIMLaAllEvt(0)        
   ,fFFHistosIMLaJet(0)           
   ,fFFHistosIMLaCone(0)
   ,fFFHistosPhiCorrIMLa(0)
   ,fALaType(0) 
   ,fFilterMaskALa(0)
   ,fListALa(0)
   ,fListFeeddownLaCand(0)
   ,fListFeeddownALaCand(0)
   ,jetConeFDLalist(0)
   ,jetConeFDALalist(0)
   ,fListMCgenK0s(0)
   ,fListMCgenLa(0)
   ,fListMCgenALa(0)
   ,fListMCgenK0sCone(0)
   ,fListMCgenLaCone(0)
   ,fListMCgenALaCone(0)
   ,IsArmenterosSelected(0)
   ,fFFHistosIMALaAllEvt(0)        
   ,fFFHistosIMALaJet(0)           
   ,fFFHistosIMALaCone(0)
   ,fFFHistosPhiCorrIMALa(0)
   ,fFFIMNBinsJetPt(0)  
   ,fFFIMJetPtMin(0) 
   ,fFFIMJetPtMax(0)
   ,fFFIMNBinsInvM(0) 
   ,fFFIMInvMMin(0)   
   ,fFFIMInvMMax(0)   
   ,fFFIMNBinsPt(0)      
   ,fFFIMPtMin(0)        
   ,fFFIMPtMax(0)        
   ,fFFIMNBinsXi(0)      
   ,fFFIMXiMin(0)        
   ,fFFIMXiMax(0)        
   ,fFFIMNBinsZ(0)       
   ,fFFIMZMin(0)         
   ,fFFIMZMax(0)
   ,fFFIMLaNBinsJetPt(0)    
   ,fFFIMLaJetPtMin(0) 
   ,fFFIMLaJetPtMax(0)
   ,fFFIMLaNBinsInvM(0) 
   ,fFFIMLaInvMMin(0)   
   ,fFFIMLaInvMMax(0)   
   ,fFFIMLaNBinsPt(0)      
   ,fFFIMLaPtMin(0)        
   ,fFFIMLaPtMax(0)        
   ,fFFIMLaNBinsXi(0)      
   ,fFFIMLaXiMin(0)        
   ,fFFIMLaXiMax(0)        
   ,fFFIMLaNBinsZ(0)       
   ,fFFIMLaZMin(0)         
   ,fFFIMLaZMax(0)
   ,fPhiCorrIMNBinsPt(0)  
   ,fPhiCorrIMPtMin(0)
   ,fPhiCorrIMPtMax(0)
   ,fPhiCorrIMNBinsPhi(0)
   ,fPhiCorrIMPhiMin(0)
   ,fPhiCorrIMPhiMax(0)
   ,fPhiCorrIMNBinsInvM(0)
   ,fPhiCorrIMInvMMin(0)
   ,fPhiCorrIMInvMMax(0)
   ,fPhiCorrIMLaNBinsPt(0)  
   ,fPhiCorrIMLaPtMin(0)
   ,fPhiCorrIMLaPtMax(0)
   ,fPhiCorrIMLaNBinsPhi(0)
   ,fPhiCorrIMLaPhiMin(0)
   ,fPhiCorrIMLaPhiMax(0)
   ,fPhiCorrIMLaNBinsInvM(0)
   ,fPhiCorrIMLaInvMMin(0)
   ,fPhiCorrIMLaInvMMax(0)
   ,fh1EvtAllCent(0)
   ,fh1Evt(0)
   ,fh1K0Mult(0)
   ,fh1dPhiJetK0(0)
   ,fh1LaMult(0)
   ,fh1dPhiJetLa(0)
   ,fh1ALaMult(0)
   ,fh1dPhiJetALa(0)
   ,fh1JetEta(0)         
   ,fh1JetPhi(0)                 
   ,fh2JetEtaPhi(0)
   ,fh1V0JetPt(0)
   ,fh2FFJetTrackEta(0)
   ,fh1trackPosNCls(0)           
   ,fh1trackNegNCls(0)   
   ,fh1trackPosRap(0)            
   ,fh1trackNegRap(0)          
   ,fh1V0Rap(0)        
   ,fh1trackPosEta(0)            
   ,fh1trackNegEta(0)          
   ,fh1V0Eta(0)
   ,fh1V0totMom(0)           
   ,fh1CosPointAngle(0)        
   ,fh1Chi2Pos(0)                 
   ,fh1Chi2Neg(0)   
   ,fh1DecayLengthV0(0)    
   ,fh2ProperLifetimeK0sVsPtBeforeCut(0)    
   ,fh2ProperLifetimeK0sVsPtAfterCut(0)
   ,fh1ProperLifetimeV0BeforeCut(0) 
   ,fh1ProperLifetimeV0AfterCut(0) 
   ,fh1V0Radius(0)     
   ,fh1DcaV0Daughters(0)        
   ,fh1DcaPosToPrimVertex(0)   
   ,fh1DcaNegToPrimVertex(0)    
   ,fh2ArmenterosBeforeCuts(0)
   ,fh2ArmenterosAfterCuts(0)
   ,fh2BB3SigProton(0)
   ,fh2BBLaPos(0)
   ,fh2BBLaNeg(0)
   ,fh1CrossedRowsOverFindableNeg(0)
   ,fh1CrossedRowsOverFindablePos(0)
   ,fh1PosDaughterCharge(0)
   ,fh1NegDaughterCharge(0)
   ,fh1PtMCK0s(0)
   ,fh1PtMCLa(0)
   ,fh1PtMCALa(0)
   ,fh1EtaK0s(0)
   ,fh1EtaLa(0)
   ,fh1EtaALa(0)
   ,fh3InvMassEtaTrackPtK0s(0)
   ,fh3InvMassEtaTrackPtLa(0)
   ,fh3InvMassEtaTrackPtALa(0)
   ,fh1noAssociatedK0s(0)
   ,fh1TrackMultCone(0)
   ,fh2TrackMultCone(0) 
   ,fh2MCgenK0Cone(0)
   ,fh2MCgenLaCone(0)
   ,fh2MCgenALaCone(0) 
   ,fh2MCEtagenK0Cone(0)
   ,fh2MCEtagenLaCone(0)
   ,fh2MCEtagenALaCone(0)
   ,fh1FFIMK0ConeSmear(0)
   ,fh1FFIMLaConeSmear(0)
   ,fh1FFIMALaConeSmear(0)
   ,fh3MCrecK0Cone(0)   
   ,fh3MCrecLaCone(0)   
   ,fh3MCrecALaCone(0)
   ,fh3MCrecK0ConeSmear(0) 
   ,fh3MCrecLaConeSmear(0)   
   ,fh3MCrecALaConeSmear(0)
   ,fh3SecContinCone(0)
   ,fh3StrContinCone(0)
   ,fh3IMK0PerpCone(0)
   ,fh3IMLaPerpCone(0)
   ,fh3IMALaPerpCone(0)
   ,fh3IMK0MedianCone(0)
   ,fh3IMLaMedianCone(0)
   ,fh3IMALaMedianCone(0)
   ,fh1MCMultiplicityPrimary(0)
   ,fh1MCMultiplicityTracks(0)
   ,fh1MCmotherLa(0)
   ,fh1MCmotherALa(0)
   ,fh3FeedDownLa(0)
   ,fh3FeedDownALa(0)
   ,fh1MCProdRadiusK0s(0)
   ,fh1MCProdRadiusLambda(0)
   ,fh1MCProdRadiusAntiLambda(0)
   ,fh1MCPtV0s(0)
   ,fh1MCPtK0s(0) 
   ,fh1MCPtLambda(0) 
   ,fh1MCPtAntiLambda(0) 
   ,fh1MCXiPt(0)
   ,fh1MCXibarPt(0)
   ,fh2MCEtaVsPtK0s(0)
   ,fh2MCEtaVsPtLa(0)
   ,fh2MCEtaVsPtALa(0)
   ,fh1MCRapK0s(0) 
   ,fh1MCRapLambda(0)
   ,fh1MCRapAntiLambda(0)
   ,fh1MCEtaAllK0s(0) 
   ,fh1MCEtaK0s(0) 
   ,fh1MCEtaLambda(0)
   ,fh1MCEtaAntiLambda(0)

{
   // default constructor
}

//__________________________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem(const char *name) 
  : AliAnalysisTaskFragmentationFunction(name)

  ,fAnalysisMC(0)
  ,fDeltaVertexZ(0)
  ,fCuttrackNegNcls(0)
  ,fCuttrackPosNcls(0)
  ,fCutPostrackRap(0)
  ,fCutNegtrackRap(0)
  ,fCutRap(0)
  ,fCutPostrackEta(0)
  ,fCutNegtrackEta(0)
  ,fCutEta(0)
  ,fCutV0cosPointAngle(0)
  ,fCutChi2PosDaughter(0)
  ,fCutChi2NegDaughter(0)
  ,fKinkDaughters(0)
  ,fRequireTPCRefit(0)
  ,fCutArmenteros(0)
  ,fCutV0DecayMin(0)
  ,fCutV0DecayMax(0)
  ,fCutV0totMom(0)
  ,fCutDcaV0Daughters(0)
  ,fCutDcaPosToPrimVertex(0)
  ,fCutDcaNegToPrimVertex(0)
  ,fCutV0RadiusMin(0)
  ,fCutV0RadiusMax(0)
  ,fCutBetheBloch(0)
  ,fCutRatio(0)  
  ,fK0Type(0)  
  ,fFilterMaskK0(0)
  ,fListK0s(0)
  ,fPIDResponse(0)
  ,fV0QAK0(0)
  ,fFFHistosRecCutsK0Evt(0)      
  ,fFFHistosIMK0AllEvt(0)        
  ,fFFHistosIMK0Jet(0)           
  ,fFFHistosIMK0Cone(0)
  ,fFFHistosPhiCorrIMK0(0)
  ,fLaType(0)  
  ,fFilterMaskLa(0)
  ,fListLa(0)
  ,fFFHistosIMLaAllEvt(0)        
  ,fFFHistosIMLaJet(0)           
  ,fFFHistosIMLaCone(0)
  ,fFFHistosPhiCorrIMLa(0)
  ,fALaType(0)  
  ,fFilterMaskALa(0)
  ,fListALa(0)
  ,fListFeeddownLaCand(0)
  ,fListFeeddownALaCand(0)
  ,jetConeFDLalist(0)
  ,jetConeFDALalist(0)
  ,fListMCgenK0s(0)
  ,fListMCgenLa(0)
  ,fListMCgenALa(0)
  ,fListMCgenK0sCone(0)
  ,fListMCgenLaCone(0)
  ,fListMCgenALaCone(0)
  ,IsArmenterosSelected(0)
  ,fFFHistosIMALaAllEvt(0)        
  ,fFFHistosIMALaJet(0)           
  ,fFFHistosIMALaCone(0)
  ,fFFHistosPhiCorrIMALa(0)
  ,fFFIMNBinsJetPt(0)    
  ,fFFIMJetPtMin(0) 
  ,fFFIMJetPtMax(0)
  ,fFFIMNBinsInvM(0) 
  ,fFFIMInvMMin(0)   
  ,fFFIMInvMMax(0)   
  ,fFFIMNBinsPt(0)      
  ,fFFIMPtMin(0)        
  ,fFFIMPtMax(0)        
  ,fFFIMNBinsXi(0)      
  ,fFFIMXiMin(0)        
  ,fFFIMXiMax(0)        
  ,fFFIMNBinsZ(0)       
  ,fFFIMZMin(0)         
  ,fFFIMZMax(0) 
  ,fFFIMLaNBinsJetPt(0)    
  ,fFFIMLaJetPtMin(0) 
  ,fFFIMLaJetPtMax(0)
  ,fFFIMLaNBinsInvM(0) 
  ,fFFIMLaInvMMin(0)   
  ,fFFIMLaInvMMax(0)   
  ,fFFIMLaNBinsPt(0)      
  ,fFFIMLaPtMin(0)        
  ,fFFIMLaPtMax(0)        
  ,fFFIMLaNBinsXi(0)      
  ,fFFIMLaXiMin(0)        
  ,fFFIMLaXiMax(0)        
  ,fFFIMLaNBinsZ(0)       
  ,fFFIMLaZMin(0)         
  ,fFFIMLaZMax(0)
  ,fPhiCorrIMNBinsPt(0)   
  ,fPhiCorrIMPtMin(0)
  ,fPhiCorrIMPtMax(0)
  ,fPhiCorrIMNBinsPhi(0)
  ,fPhiCorrIMPhiMin(0)
  ,fPhiCorrIMPhiMax(0)
  ,fPhiCorrIMNBinsInvM(0)
  ,fPhiCorrIMInvMMin(0)
  ,fPhiCorrIMInvMMax(0)
  ,fPhiCorrIMLaNBinsPt(0)   
  ,fPhiCorrIMLaPtMin(0)
  ,fPhiCorrIMLaPtMax(0)
  ,fPhiCorrIMLaNBinsPhi(0)
  ,fPhiCorrIMLaPhiMin(0)
  ,fPhiCorrIMLaPhiMax(0)
  ,fPhiCorrIMLaNBinsInvM(0)
  ,fPhiCorrIMLaInvMMin(0)
  ,fPhiCorrIMLaInvMMax(0)
  ,fh1EvtAllCent(0)
  ,fh1Evt(0)
  ,fh1K0Mult(0)
  ,fh1dPhiJetK0(0) 
  ,fh1LaMult(0)
  ,fh1dPhiJetLa(0) 
  ,fh1ALaMult(0)
  ,fh1dPhiJetALa(0)  
  ,fh1JetEta(0)         
  ,fh1JetPhi(0)                 
  ,fh2JetEtaPhi(0)
  ,fh1V0JetPt(0)
  ,fh2FFJetTrackEta(0)  
  ,fh1trackPosNCls(0)           
  ,fh1trackNegNCls(0) 
  ,fh1trackPosRap(0)            
  ,fh1trackNegRap(0)          
  ,fh1V0Rap(0)          
  ,fh1trackPosEta(0)            
  ,fh1trackNegEta(0)          
  ,fh1V0Eta(0)  
  ,fh1V0totMom(0)            
  ,fh1CosPointAngle(0)        
  ,fh1Chi2Pos(0)                 
  ,fh1Chi2Neg(0) 
  ,fh1DecayLengthV0(0) 
  ,fh2ProperLifetimeK0sVsPtBeforeCut(0)  
  ,fh2ProperLifetimeK0sVsPtAfterCut(0)            
  ,fh1ProperLifetimeV0BeforeCut(0)  
  ,fh1ProperLifetimeV0AfterCut(0)  
  ,fh1V0Radius(0)       
  ,fh1DcaV0Daughters(0)        
  ,fh1DcaPosToPrimVertex(0)   
  ,fh1DcaNegToPrimVertex(0)    
  ,fh2ArmenterosBeforeCuts(0)
  ,fh2ArmenterosAfterCuts(0)
  ,fh2BB3SigProton(0)
  ,fh2BBLaPos(0)
  ,fh2BBLaNeg(0)
  ,fh1CrossedRowsOverFindableNeg(0)
  ,fh1CrossedRowsOverFindablePos(0)
  ,fh1PosDaughterCharge(0)
  ,fh1NegDaughterCharge(0)
  ,fh1PtMCK0s(0)
  ,fh1PtMCLa(0)
  ,fh1PtMCALa(0)
  ,fh1EtaK0s(0)
  ,fh1EtaLa(0)
  ,fh1EtaALa(0)
  ,fh3InvMassEtaTrackPtK0s(0)
  ,fh3InvMassEtaTrackPtLa(0)
  ,fh3InvMassEtaTrackPtALa(0)
  ,fh1noAssociatedK0s(0)
  ,fh1TrackMultCone(0)
  ,fh2TrackMultCone(0)
  ,fh2MCgenK0Cone(0)
  ,fh2MCgenLaCone(0)
  ,fh2MCgenALaCone(0)
  ,fh2MCEtagenK0Cone(0)
  ,fh2MCEtagenLaCone(0)
  ,fh2MCEtagenALaCone(0)
  ,fh1FFIMK0ConeSmear(0)
  ,fh1FFIMLaConeSmear(0)
  ,fh1FFIMALaConeSmear(0)
  ,fh3MCrecK0Cone(0)
  ,fh3MCrecLaCone(0)
  ,fh3MCrecALaCone(0) 
  ,fh3MCrecK0ConeSmear(0) 
  ,fh3MCrecLaConeSmear(0)   
  ,fh3MCrecALaConeSmear(0)
  ,fh3SecContinCone(0)
  ,fh3StrContinCone(0)
  ,fh3IMK0PerpCone(0)
  ,fh3IMLaPerpCone(0)
  ,fh3IMALaPerpCone(0)
  ,fh3IMK0MedianCone(0)
  ,fh3IMLaMedianCone(0)
  ,fh3IMALaMedianCone(0)
  ,fh1MCMultiplicityPrimary(0)
  ,fh1MCMultiplicityTracks(0)
  ,fh1MCmotherLa(0)
  ,fh1MCmotherALa(0)
  ,fh3FeedDownLa(0)
  ,fh3FeedDownALa(0)
  ,fh1MCProdRadiusK0s(0)
  ,fh1MCProdRadiusLambda(0)
  ,fh1MCProdRadiusAntiLambda(0)
  ,fh1MCPtV0s(0)
  ,fh1MCPtK0s(0)
  ,fh1MCPtLambda(0) 
  ,fh1MCPtAntiLambda(0) 
  ,fh1MCXiPt(0)
  ,fh1MCXibarPt(0)
  ,fh2MCEtaVsPtK0s(0)
  ,fh2MCEtaVsPtLa(0)
  ,fh2MCEtaVsPtALa(0)
  ,fh1MCRapK0s(0) 
  ,fh1MCRapLambda(0)
  ,fh1MCRapAntiLambda(0)
  ,fh1MCEtaAllK0s(0) 
  ,fh1MCEtaK0s(0) 
  ,fh1MCEtaLambda(0)
  ,fh1MCEtaAntiLambda(0)


{
  // constructor
  
  DefineOutput(1,TList::Class());  
}

//__________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem(const  AliAnalysisTaskJetChem &copy)
  : AliAnalysisTaskFragmentationFunction()
  
  ,fAnalysisMC(copy.fAnalysisMC)
  ,fDeltaVertexZ(copy.fDeltaVertexZ)
  ,fCuttrackNegNcls(copy.fCuttrackNegNcls)
  ,fCuttrackPosNcls(copy.fCuttrackPosNcls)
  ,fCutPostrackRap(copy.fCutPostrackRap)
  ,fCutNegtrackRap(copy.fCutNegtrackRap)
  ,fCutRap(copy.fCutRap)
  ,fCutPostrackEta(copy.fCutPostrackEta)
  ,fCutNegtrackEta(copy.fCutNegtrackEta)
  ,fCutEta(copy.fCutEta)
  ,fCutV0cosPointAngle(copy.fCutV0cosPointAngle)
  ,fCutChi2PosDaughter(copy.fCutChi2PosDaughter)
  ,fCutChi2NegDaughter(copy.fCutChi2NegDaughter)
  ,fKinkDaughters(copy.fKinkDaughters)
  ,fRequireTPCRefit(copy.fRequireTPCRefit)
  ,fCutArmenteros(copy.fCutArmenteros)
  ,fCutV0DecayMin(copy.fCutV0DecayMin)
  ,fCutV0DecayMax(copy.fCutV0DecayMax)
  ,fCutV0totMom(copy.fCutV0totMom)
  ,fCutDcaV0Daughters(copy.fCutDcaV0Daughters)
  ,fCutDcaPosToPrimVertex(copy.fCutDcaPosToPrimVertex)
  ,fCutDcaNegToPrimVertex(copy.fCutDcaNegToPrimVertex)
  ,fCutV0RadiusMin(copy.fCutV0RadiusMin)
  ,fCutV0RadiusMax(copy.fCutV0RadiusMax)
  ,fCutBetheBloch(copy.fCutBetheBloch)
  ,fCutRatio(copy.fCutRatio)
  ,fK0Type(copy.fK0Type)              
  ,fFilterMaskK0(copy.fFilterMaskK0)
  ,fListK0s(copy.fListK0s)
  ,fPIDResponse(copy.fPIDResponse)
  ,fV0QAK0(copy.fV0QAK0)
  ,fFFHistosRecCutsK0Evt(copy.fFFHistosRecCutsK0Evt)      
  ,fFFHistosIMK0AllEvt(copy.fFFHistosIMK0AllEvt)        
  ,fFFHistosIMK0Jet(copy.fFFHistosIMK0Jet)           
  ,fFFHistosIMK0Cone(copy.fFFHistosIMK0Cone)          
  ,fFFHistosPhiCorrIMK0(copy.fFFHistosPhiCorrIMK0)     
  ,fLaType(copy.fLaType)                  
  ,fFilterMaskLa(copy.fFilterMaskLa)
  ,fListLa(copy.fListLa)
  ,fFFHistosIMLaAllEvt(copy.fFFHistosIMLaAllEvt)        
  ,fFFHistosIMLaJet(copy.fFFHistosIMLaJet)           
  ,fFFHistosIMLaCone(copy.fFFHistosIMLaCone)          
  ,fFFHistosPhiCorrIMLa(copy.fFFHistosPhiCorrIMLa) 
  ,fALaType(copy.fALaType)                 
  ,fFilterMaskALa(copy.fFilterMaskALa)
  ,fListALa(copy.fListALa)
  ,fListFeeddownLaCand(copy.fListFeeddownLaCand)
  ,fListFeeddownALaCand(copy.fListFeeddownALaCand)
  ,jetConeFDLalist(copy.jetConeFDLalist)
  ,jetConeFDALalist(copy.jetConeFDALalist)
  ,fListMCgenK0s(copy.fListMCgenK0s)
  ,fListMCgenLa(copy.fListMCgenLa)
  ,fListMCgenALa(copy.fListMCgenALa)
  ,fListMCgenK0sCone(copy.fListMCgenK0sCone)
  ,fListMCgenLaCone(copy.fListMCgenLaCone)
  ,fListMCgenALaCone(copy.fListMCgenALaCone)
  ,IsArmenterosSelected(copy.IsArmenterosSelected)
  ,fFFHistosIMALaAllEvt(copy.fFFHistosIMALaAllEvt)        
  ,fFFHistosIMALaJet(copy.fFFHistosIMALaJet)           
  ,fFFHistosIMALaCone(copy.fFFHistosIMALaCone)          
  ,fFFHistosPhiCorrIMALa(copy.fFFHistosPhiCorrIMALa) 
  ,fFFIMNBinsJetPt(copy.fFFIMNBinsJetPt) 
  ,fFFIMJetPtMin(copy.fFFIMJetPtMin)     
  ,fFFIMJetPtMax(copy.fFFIMJetPtMax)     
  ,fFFIMNBinsInvM(copy.fFFIMNBinsInvM)  
  ,fFFIMInvMMin(copy.fFFIMInvMMin)    
  ,fFFIMInvMMax(copy.fFFIMInvMMax)    
  ,fFFIMNBinsPt(copy.fFFIMNBinsPt)      
  ,fFFIMPtMin(copy.fFFIMPtMin)        
  ,fFFIMPtMax(copy.fFFIMPtMax)        
  ,fFFIMNBinsXi(copy.fFFIMNBinsXi)      
  ,fFFIMXiMin(copy.fFFIMXiMin)        
  ,fFFIMXiMax(copy.fFFIMXiMax)        
  ,fFFIMNBinsZ(copy.fFFIMNBinsZ)       
  ,fFFIMZMin(copy.fFFIMZMin)         
  ,fFFIMZMax(copy.fFFIMZMax) 
  ,fFFIMLaNBinsJetPt(copy.fFFIMLaNBinsJetPt)   
  ,fFFIMLaJetPtMin(copy.fFFIMLaJetPtMin)     
  ,fFFIMLaJetPtMax(copy.fFFIMLaJetPtMax)     
  ,fFFIMLaNBinsInvM(copy.fFFIMLaNBinsInvM)  
  ,fFFIMLaInvMMin(copy.fFFIMLaInvMMin)    
  ,fFFIMLaInvMMax(copy.fFFIMLaInvMMax)    
  ,fFFIMLaNBinsPt(copy.fFFIMLaNBinsPt)      
  ,fFFIMLaPtMin(copy.fFFIMLaPtMin)        
  ,fFFIMLaPtMax(copy.fFFIMLaPtMax)        
  ,fFFIMLaNBinsXi(copy.fFFIMLaNBinsXi)      
  ,fFFIMLaXiMin(copy.fFFIMLaXiMin)        
  ,fFFIMLaXiMax(copy.fFFIMLaXiMax)        
  ,fFFIMLaNBinsZ(copy.fFFIMLaNBinsZ)       
  ,fFFIMLaZMin(copy.fFFIMLaZMin)         
  ,fFFIMLaZMax(copy.fFFIMLaZMax) 
  ,fPhiCorrIMNBinsPt(copy.fPhiCorrIMNBinsPt)  
  ,fPhiCorrIMPtMin(copy.fPhiCorrIMPtMin)
  ,fPhiCorrIMPtMax(copy.fPhiCorrIMPtMax)
  ,fPhiCorrIMNBinsPhi(copy.fPhiCorrIMNBinsPhi)
  ,fPhiCorrIMPhiMin(copy.fPhiCorrIMPhiMin)
  ,fPhiCorrIMPhiMax(copy.fPhiCorrIMPhiMax)
  ,fPhiCorrIMNBinsInvM(copy.fPhiCorrIMNBinsInvM)
  ,fPhiCorrIMInvMMin(copy.fPhiCorrIMInvMMin)
  ,fPhiCorrIMInvMMax(copy.fPhiCorrIMInvMMax)
  ,fPhiCorrIMLaNBinsPt(copy.fPhiCorrIMLaNBinsPt)   
  ,fPhiCorrIMLaPtMin(copy.fPhiCorrIMLaPtMin)
  ,fPhiCorrIMLaPtMax(copy.fPhiCorrIMLaPtMax)
  ,fPhiCorrIMLaNBinsPhi(copy.fPhiCorrIMLaNBinsPhi)
  ,fPhiCorrIMLaPhiMin(copy.fPhiCorrIMLaPhiMin)
  ,fPhiCorrIMLaPhiMax(copy.fPhiCorrIMLaPhiMax)
  ,fPhiCorrIMLaNBinsInvM(copy.fPhiCorrIMLaNBinsInvM)
  ,fPhiCorrIMLaInvMMin(copy.fPhiCorrIMLaInvMMin)
  ,fPhiCorrIMLaInvMMax(copy.fPhiCorrIMLaInvMMax)
  ,fh1EvtAllCent(copy.fh1EvtAllCent)
  ,fh1Evt(copy.fh1Evt)
  ,fh1K0Mult(copy.fh1K0Mult)
  ,fh1dPhiJetK0(copy.fh1dPhiJetK0)
  ,fh1LaMult(copy.fh1LaMult)
  ,fh1dPhiJetLa(copy.fh1dPhiJetLa)
  ,fh1ALaMult(copy.fh1ALaMult)
  ,fh1dPhiJetALa(copy.fh1dPhiJetALa)
  ,fh1JetEta(copy.fh1JetEta)         
  ,fh1JetPhi(copy.fh1JetPhi)                 
  ,fh2JetEtaPhi(copy.fh2JetEtaPhi)
  ,fh1V0JetPt(copy.fh1V0JetPt)
  ,fh2FFJetTrackEta(copy.fh2FFJetTrackEta) 
  ,fh1trackPosNCls(copy.fh1trackPosNCls)           
  ,fh1trackNegNCls(copy.fh1trackNegNCls)
  ,fh1trackPosRap(copy.fh1trackPosRap)            
  ,fh1trackNegRap(copy.fh1trackNegRap)          
  ,fh1V0Rap(copy.fh1V0Rap)         
  ,fh1trackPosEta(copy.fh1trackPosEta)            
  ,fh1trackNegEta(copy.fh1trackNegEta)          
  ,fh1V0Eta(copy.fh1V0Eta)   
  ,fh1V0totMom(copy.fh1V0totMom)           
  ,fh1CosPointAngle(copy.fh1CosPointAngle)        
  ,fh1Chi2Pos(copy.fh1Chi2Pos)                 
  ,fh1Chi2Neg(copy.fh1Chi2Neg)    
  ,fh1DecayLengthV0(copy.fh1DecayLengthV0)  
  ,fh2ProperLifetimeK0sVsPtBeforeCut(copy.fh2ProperLifetimeK0sVsPtBeforeCut)  
  ,fh2ProperLifetimeK0sVsPtAfterCut(copy.fh2ProperLifetimeK0sVsPtAfterCut)    
  ,fh1ProperLifetimeV0BeforeCut(copy.fh1ProperLifetimeV0BeforeCut) 
  ,fh1ProperLifetimeV0AfterCut(copy.fh1ProperLifetimeV0AfterCut) 
  ,fh1V0Radius(copy.fh1V0Radius)          
  ,fh1DcaV0Daughters(copy.fh1DcaV0Daughters)        
  ,fh1DcaPosToPrimVertex(copy.fh1DcaPosToPrimVertex)   
  ,fh1DcaNegToPrimVertex(copy.fh1DcaNegToPrimVertex)    
  ,fh2ArmenterosBeforeCuts(copy.fh2ArmenterosBeforeCuts)
  ,fh2ArmenterosAfterCuts(copy.fh2ArmenterosAfterCuts)
  ,fh2BB3SigProton(copy.fh2BB3SigProton)
  ,fh2BBLaPos(copy.fh2BBLaPos)
  ,fh2BBLaNeg(copy.fh2BBLaPos)
  ,fh1CrossedRowsOverFindableNeg(copy.fh1CrossedRowsOverFindableNeg)
  ,fh1CrossedRowsOverFindablePos(copy.fh1CrossedRowsOverFindablePos)
  ,fh1PosDaughterCharge(copy.fh1PosDaughterCharge)
  ,fh1NegDaughterCharge(copy.fh1NegDaughterCharge)
  ,fh1PtMCK0s(copy.fh1PtMCK0s)
  ,fh1PtMCLa(copy.fh1PtMCLa)
  ,fh1PtMCALa(copy.fh1PtMCALa)
  ,fh1EtaK0s(copy.fh1EtaK0s)
  ,fh1EtaLa(copy.fh1EtaLa)
  ,fh1EtaALa(copy.fh1EtaALa)
  ,fh3InvMassEtaTrackPtK0s(copy.fh3InvMassEtaTrackPtK0s)
  ,fh3InvMassEtaTrackPtLa(copy.fh3InvMassEtaTrackPtLa)
  ,fh3InvMassEtaTrackPtALa(copy.fh3InvMassEtaTrackPtALa)
  ,fh1noAssociatedK0s(copy.fh1noAssociatedK0s)
  ,fh1TrackMultCone(copy.fh1TrackMultCone)
  ,fh2TrackMultCone(copy.fh2TrackMultCone)
  ,fh2MCgenK0Cone(copy.fh2MCgenK0Cone)
  ,fh2MCgenLaCone(copy.fh2MCgenLaCone)
  ,fh2MCgenALaCone(copy.fh2MCgenALaCone)
  ,fh2MCEtagenK0Cone(copy.fh2MCEtagenK0Cone)
  ,fh2MCEtagenLaCone(copy.fh2MCEtagenLaCone)
  ,fh2MCEtagenALaCone(copy.fh2MCEtagenALaCone)
  ,fh1FFIMK0ConeSmear(copy.fh1FFIMK0ConeSmear)
  ,fh1FFIMLaConeSmear(copy.fh1FFIMLaConeSmear)
  ,fh1FFIMALaConeSmear(copy.fh1FFIMALaConeSmear)  
  ,fh3MCrecK0Cone(copy.fh3MCrecK0Cone)
  ,fh3MCrecLaCone(copy.fh3MCrecLaCone)
  ,fh3MCrecALaCone(copy.fh3MCrecALaCone) 
  ,fh3MCrecK0ConeSmear(copy.fh3MCrecK0ConeSmear)
  ,fh3MCrecLaConeSmear(copy.fh3MCrecLaConeSmear)
  ,fh3MCrecALaConeSmear(copy.fh3MCrecALaConeSmear)
  ,fh3SecContinCone(copy.fh3SecContinCone)
  ,fh3StrContinCone(copy.fh3StrContinCone)
  ,fh3IMK0PerpCone(copy.fh3IMK0PerpCone)
  ,fh3IMLaPerpCone(copy.fh3IMLaPerpCone)
  ,fh3IMALaPerpCone(copy.fh3IMALaPerpCone)  
  ,fh3IMK0MedianCone(copy.fh3IMK0MedianCone)
  ,fh3IMLaMedianCone(copy.fh3IMLaMedianCone)
  ,fh3IMALaMedianCone(copy.fh3IMALaMedianCone)  
  ,fh1MCMultiplicityPrimary(copy.fh1MCMultiplicityPrimary)
  ,fh1MCMultiplicityTracks(copy.fh1MCMultiplicityTracks)
  ,fh1MCmotherLa(copy.fh1MCmotherLa)
  ,fh1MCmotherALa(copy.fh1MCmotherALa)
  ,fh3FeedDownLa(copy.fh3FeedDownLa)
  ,fh3FeedDownALa(copy.fh3FeedDownALa)
  ,fh1MCProdRadiusK0s(copy.fh1MCProdRadiusK0s)
  ,fh1MCProdRadiusLambda(copy.fh1MCProdRadiusLambda)
  ,fh1MCProdRadiusAntiLambda(copy.fh1MCProdRadiusAntiLambda)
  ,fh1MCPtV0s(copy.fh1MCPtV0s)
  ,fh1MCPtK0s(copy.fh1MCPtK0s) 
  ,fh1MCPtLambda(copy.fh1MCPtLambda) 
  ,fh1MCPtAntiLambda(copy.fh1MCPtAntiLambda) 
  ,fh1MCXiPt(copy.fh1MCXiPt)
  ,fh1MCXibarPt(copy.fh1MCXibarPt)
  ,fh2MCEtaVsPtK0s(copy.fh2MCEtaVsPtK0s)
  ,fh2MCEtaVsPtLa(copy.fh2MCEtaVsPtLa)
  ,fh2MCEtaVsPtALa(copy.fh2MCEtaVsPtALa)
  ,fh1MCRapK0s(copy.fh1MCRapK0s) 
  ,fh1MCRapLambda(copy.fh1MCRapLambda)
  ,fh1MCRapAntiLambda(copy.fh1MCRapAntiLambda)
  ,fh1MCEtaAllK0s(copy.fh1MCEtaAllK0s) 
  ,fh1MCEtaK0s(copy.fh1MCEtaK0s) 
  ,fh1MCEtaLambda(copy.fh1MCEtaLambda)
  ,fh1MCEtaAntiLambda(copy.fh1MCEtaAntiLambda)

{
  // copy constructor
  
}

// _________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem& AliAnalysisTaskJetChem::operator=(const AliAnalysisTaskJetChem& o)
{
  // assignment
  
  if(this!=&o){
    AliAnalysisTaskFragmentationFunction::operator=(o);

    fAnalysisMC                     = o.fAnalysisMC;
    fDeltaVertexZ                   = o.fDeltaVertexZ;
    fCuttrackNegNcls                = o.fCuttrackNegNcls;
    fCuttrackPosNcls                = o.fCuttrackPosNcls;
    fCutPostrackRap                 = o.fCutPostrackRap;
    fCutNegtrackRap                 = o.fCutNegtrackRap;  
    fCutRap                         = o.fCutRap;
    fCutPostrackEta                 = o.fCutPostrackEta;
    fCutNegtrackEta                 = o.fCutNegtrackEta;  
    fCutEta                         = o.fCutEta;
    fCutV0cosPointAngle             = o.fCutV0cosPointAngle;
    fCutChi2PosDaughter             = o.fCutChi2PosDaughter;
    fCutChi2NegDaughter             = o.fCutChi2NegDaughter;
    fKinkDaughters                  = o.fKinkDaughters;
    fRequireTPCRefit                = o.fRequireTPCRefit;
    fCutArmenteros                  = o.fCutArmenteros;
    fCutV0DecayMin                  = o.fCutV0DecayMin;
    fCutV0DecayMax                  = o.fCutV0DecayMax;
    fCutV0totMom                    = o.fCutV0totMom;
    fCutDcaV0Daughters              = o.fCutDcaV0Daughters;
    fCutDcaPosToPrimVertex          = o.fCutDcaPosToPrimVertex;
    fCutDcaNegToPrimVertex          = o.fCutDcaNegToPrimVertex;
    fCutV0RadiusMin                 = o.fCutV0RadiusMin;
    fCutV0RadiusMax                 = o.fCutV0RadiusMax;
    fCutBetheBloch                  = o.fCutBetheBloch; 
    fCutRatio                       = o.fCutRatio;
    fK0Type                         = o.fK0Type;
    fFilterMaskK0                   = o.fFilterMaskK0;
    fListK0s                        = o.fListK0s;
    fPIDResponse                    = o.fPIDResponse;
    fV0QAK0                         = o.fV0QAK0;
    fFFHistosRecCutsK0Evt           = o.fFFHistosRecCutsK0Evt;      
    fFFHistosIMK0AllEvt             = o.fFFHistosIMK0AllEvt;        
    fFFHistosIMK0Jet                = o.fFFHistosIMK0Jet;           
    fFFHistosIMK0Cone               = o.fFFHistosIMK0Cone;          
    fFFHistosPhiCorrIMK0            = o.fFFHistosPhiCorrIMK0;
    fLaType                         = o.fLaType;
    fFilterMaskLa                   = o.fFilterMaskLa;
    fListLa                         = o.fListLa;
    fFFHistosIMLaAllEvt             = o.fFFHistosIMLaAllEvt;        
    fFFHistosIMLaJet                = o.fFFHistosIMLaJet;           
    fFFHistosIMLaCone               = o.fFFHistosIMLaCone;          
    fFFHistosPhiCorrIMLa            = o.fFFHistosPhiCorrIMLa;
    fALaType                        = o.fALaType;
    fFilterMaskALa                  = o.fFilterMaskALa;
    fListFeeddownLaCand             = o.fListFeeddownLaCand;
    fListFeeddownALaCand            = o.fListFeeddownALaCand;
    jetConeFDLalist                 = o.jetConeFDLalist;
    jetConeFDALalist                = o.jetConeFDALalist;
    fListMCgenK0s                   = o.fListMCgenK0s;
    fListMCgenLa                    = o.fListMCgenLa;
    fListMCgenALa                   = o.fListMCgenALa;
    fListMCgenK0sCone               = o.fListMCgenK0sCone;
    fListMCgenLaCone                = o.fListMCgenLaCone;
    fListMCgenALaCone               = o.fListMCgenALaCone;
    IsArmenterosSelected            = o.IsArmenterosSelected;
    fFFHistosIMALaAllEvt            = o.fFFHistosIMALaAllEvt;        
    fFFHistosIMALaJet               = o.fFFHistosIMALaJet;           
    fFFHistosIMALaCone              = o.fFFHistosIMALaCone;          
    fFFHistosPhiCorrIMALa           = o.fFFHistosPhiCorrIMALa;
    fFFIMNBinsJetPt                 = o.fFFIMNBinsJetPt;   
    fFFIMJetPtMin                   = o.fFFIMJetPtMin; 
    fFFIMJetPtMax                   = o.fFFIMJetPtMax;
    fFFIMNBinsPt                    = o.fFFIMNBinsPt;      
    fFFIMPtMin                      = o.fFFIMPtMin;        
    fFFIMPtMax                      = o.fFFIMPtMax;        
    fFFIMNBinsXi                    = o.fFFIMNBinsXi;      
    fFFIMXiMin                      = o.fFFIMXiMin;        
    fFFIMXiMax                      = o.fFFIMXiMax;        
    fFFIMNBinsZ                     = o.fFFIMNBinsZ;       
    fFFIMZMin                       = o.fFFIMZMin;         
    fFFIMZMax                       = o.fFFIMZMax;  
    fFFIMLaNBinsJetPt               = o.fFFIMLaNBinsJetPt;    
    fFFIMLaJetPtMin                 = o.fFFIMLaJetPtMin; 
    fFFIMLaJetPtMax                 = o.fFFIMLaJetPtMax;
    fFFIMLaNBinsPt                  = o.fFFIMLaNBinsPt;      
    fFFIMLaPtMin                    = o.fFFIMLaPtMin;        
    fFFIMLaPtMax                    = o.fFFIMLaPtMax;        
    fFFIMLaNBinsXi                  = o.fFFIMLaNBinsXi;      
    fFFIMLaXiMin                    = o.fFFIMLaXiMin;        
    fFFIMLaXiMax                    = o.fFFIMLaXiMax;        
    fFFIMLaNBinsZ                   = o.fFFIMLaNBinsZ;       
    fFFIMLaZMin                     = o.fFFIMLaZMin;         
    fFFIMLaZMax                     = o.fFFIMLaZMax;
    fPhiCorrIMNBinsPt               = o.fPhiCorrIMNBinsPt;
    fPhiCorrIMPtMin                 = o.fPhiCorrIMPtMin;
    fPhiCorrIMPtMax                 = o.fPhiCorrIMPtMax;
    fPhiCorrIMNBinsPhi              = o.fPhiCorrIMNBinsPhi;
    fPhiCorrIMPhiMin                = o.fPhiCorrIMPhiMin;
    fPhiCorrIMPhiMax                = o.fPhiCorrIMPhiMax;
    fPhiCorrIMNBinsInvM             = o.fPhiCorrIMNBinsInvM;
    fPhiCorrIMInvMMin               = o.fPhiCorrIMInvMMin;
    fPhiCorrIMInvMMax               = o.fPhiCorrIMInvMMax;
    fPhiCorrIMLaNBinsPt             = o.fPhiCorrIMLaNBinsPt;   
    fPhiCorrIMLaPtMin               = o.fPhiCorrIMLaPtMin;
    fPhiCorrIMLaPtMax               = o.fPhiCorrIMLaPtMax;
    fPhiCorrIMLaNBinsPhi            = o.fPhiCorrIMLaNBinsPhi;
    fPhiCorrIMLaPhiMin              = o.fPhiCorrIMLaPhiMin;
    fPhiCorrIMLaPhiMax              = o.fPhiCorrIMLaPhiMax;
    fPhiCorrIMLaNBinsInvM           = o.fPhiCorrIMLaNBinsInvM;
    fPhiCorrIMLaInvMMin             = o.fPhiCorrIMLaInvMMin;
    fPhiCorrIMLaInvMMax             = o.fPhiCorrIMLaInvMMax;
    fh1EvtAllCent                   = o.fh1EvtAllCent;
    fh1Evt                          = o.fh1Evt;
    fh1K0Mult                       = o.fh1K0Mult;
    fh1dPhiJetK0                    = o.fh1dPhiJetK0;
    fh1LaMult                       = o.fh1LaMult;
    fh1dPhiJetLa                    = o.fh1dPhiJetLa;
    fh1ALaMult                      = o.fh1ALaMult;
    fh1dPhiJetALa                   = o.fh1dPhiJetALa;
    fh1JetEta                       = o.fh1JetEta;         
    fh1JetPhi                       = o.fh1JetPhi;                 
    fh2JetEtaPhi                    = o.fh2JetEtaPhi;
    fh1V0JetPt                     = o.fh1V0JetPt;
    fh2FFJetTrackEta                = o.fh2FFJetTrackEta; 
    fh1trackPosNCls                 = o.fh1trackPosNCls;           
    fh1trackNegNCls                 = o.fh1trackNegNCls;    
    fh1trackPosRap                  = o.fh1trackPosRap;            
    fh1trackNegRap                  = o.fh1trackNegRap;        
    fh1V0Rap                        = o.fh1V0Rap;        
    fh1trackPosEta                  = o.fh1trackPosEta;            
    fh1trackNegEta                  = o.fh1trackNegEta;        
    fh1V0Eta                        = o.fh1V0Eta;  
    fh1V0totMom                     = o.fh1V0totMom;            
    fh1CosPointAngle                = o.fh1CosPointAngle;        
    fh1Chi2Pos                      = o.fh1Chi2Pos;                 
    fh1Chi2Neg                      = o.fh1Chi2Neg;              
    fh1DecayLengthV0                = o.fh1DecayLengthV0;  
    fh2ProperLifetimeK0sVsPtBeforeCut = o.fh2ProperLifetimeK0sVsPtBeforeCut;
    fh2ProperLifetimeK0sVsPtAfterCut= o.fh2ProperLifetimeK0sVsPtAfterCut; 
    fh1ProperLifetimeV0BeforeCut    = o.fh1ProperLifetimeV0BeforeCut; 
    fh1ProperLifetimeV0AfterCut     = o.fh1ProperLifetimeV0AfterCut; 
    fh1V0Radius                     = o.fh1V0Radius;         
    fh1DcaV0Daughters               = o.fh1DcaV0Daughters;        
    fh1DcaPosToPrimVertex           = o.fh1DcaPosToPrimVertex;   
    fh1DcaNegToPrimVertex           = o.fh1DcaNegToPrimVertex;    
    fh2ArmenterosBeforeCuts         = o.fh2ArmenterosBeforeCuts;
    fh2ArmenterosAfterCuts          = o.fh2ArmenterosAfterCuts;
    fh2BB3SigProton                 = o.fh2BB3SigProton;
    fh2BBLaPos                      = o.fh2BBLaPos;
    fh2BBLaNeg                      = o.fh2BBLaPos;
    fh1CrossedRowsOverFindableNeg   = o.fh1CrossedRowsOverFindableNeg;
    fh1CrossedRowsOverFindablePos   = o.fh1CrossedRowsOverFindablePos;
    fh1PosDaughterCharge            = o.fh1PosDaughterCharge;
    fh1NegDaughterCharge            = o.fh1NegDaughterCharge;
    fh1PtMCK0s                      = o.fh1PtMCK0s;
    fh1PtMCLa                       = o.fh1PtMCLa;
    fh1PtMCALa                      = o.fh1PtMCALa;
    fh1EtaK0s                       = o.fh1EtaK0s;
    fh1EtaLa                        = o.fh1EtaLa;
    fh1EtaALa                       = o.fh1EtaALa;
    fh3InvMassEtaTrackPtK0s         = o.fh3InvMassEtaTrackPtK0s;
    fh3InvMassEtaTrackPtLa          = o.fh3InvMassEtaTrackPtLa;
    fh3InvMassEtaTrackPtALa         = o.fh3InvMassEtaTrackPtALa;
    fh1noAssociatedK0s              = o.fh1noAssociatedK0s;
    fh1TrackMultCone                = o.fh1TrackMultCone;
    fh2TrackMultCone                = o.fh2TrackMultCone;
    fh2MCgenK0Cone                  = o.fh2MCgenK0Cone;
    fh2MCgenLaCone                  = o.fh2MCgenLaCone;
    fh2MCgenALaCone                 = o.fh2MCgenALaCone; 
    fh2MCEtagenK0Cone               = o.fh2MCEtagenK0Cone;
    fh2MCEtagenLaCone               = o.fh2MCEtagenLaCone;
    fh2MCEtagenALaCone              = o.fh2MCEtagenALaCone;
    fh1FFIMK0ConeSmear              = o.fh1FFIMK0ConeSmear;
    fh1FFIMLaConeSmear              = o.fh1FFIMLaConeSmear;
    fh1FFIMALaConeSmear             = o.fh1FFIMALaConeSmear;
    fh3MCrecK0Cone                  = o.fh3MCrecK0Cone;
    fh3MCrecLaCone                  = o.fh3MCrecLaCone;
    fh3MCrecALaCone                 = o.fh3MCrecALaCone;
    fh3MCrecK0ConeSmear             = o.fh3MCrecK0ConeSmear;
    fh3MCrecLaConeSmear             = o.fh3MCrecLaConeSmear;
    fh3MCrecALaConeSmear            = o.fh3MCrecALaConeSmear;
    fh3SecContinCone                = o.fh3SecContinCone;
    fh3StrContinCone                = o.fh3StrContinCone;
    fh3IMK0PerpCone                 = o.fh3IMK0PerpCone;
    fh3IMLaPerpCone                 = o.fh3IMLaPerpCone;
    fh3IMALaPerpCone                = o.fh3IMALaPerpCone;
    fh3IMK0MedianCone               = o.fh3IMK0MedianCone;
    fh3IMLaMedianCone               = o.fh3IMLaMedianCone;
    fh3IMALaMedianCone              = o.fh3IMALaMedianCone; 
    fh1MCMultiplicityPrimary        = o.fh1MCMultiplicityPrimary;
    fh1MCMultiplicityTracks         = o.fh1MCMultiplicityTracks;
    fh1MCmotherLa                   = o.fh1MCmotherLa;
    fh1MCmotherALa                  = o.fh1MCmotherALa;
    fh3FeedDownLa                   = o.fh3FeedDownLa;
    fh3FeedDownALa                  = o.fh3FeedDownALa;
    fh1MCProdRadiusK0s              = o.fh1MCProdRadiusK0s;
    fh1MCProdRadiusLambda           = o.fh1MCProdRadiusLambda;
    fh1MCProdRadiusAntiLambda       = o.fh1MCProdRadiusAntiLambda;
    fh1MCPtV0s                      = o.fh1MCPtV0s;
    fh1MCPtK0s                      = o.fh1MCPtK0s; 
    fh1MCPtLambda                   = o.fh1MCPtLambda;
    fh1MCPtAntiLambda               = o.fh1MCPtAntiLambda; 
    fh1MCXiPt                       = o.fh1MCXiPt;
    fh1MCXibarPt                    = o.fh1MCXibarPt;
    fh2MCEtaVsPtK0s                 = o.fh2MCEtaVsPtK0s;
    fh2MCEtaVsPtLa                  = o.fh2MCEtaVsPtLa;
    fh2MCEtaVsPtALa                 = o.fh2MCEtaVsPtALa;
    fh1MCRapK0s                     = o.fh1MCRapK0s; 
    fh1MCRapLambda                  = o.fh1MCRapLambda;
    fh1MCRapAntiLambda              = o.fh1MCRapAntiLambda;
    fh1MCEtaAllK0s                  = o.fh1MCEtaAllK0s; 
    fh1MCEtaK0s                     = o.fh1MCEtaK0s; 
    fh1MCEtaLambda                  = o.fh1MCEtaLambda;
    fh1MCEtaAntiLambda              = o.fh1MCEtaAntiLambda;
}
    
  return *this;
}

//_______________________________________________
AliAnalysisTaskJetChem::~AliAnalysisTaskJetChem()
{
  // destructor  


  if(fListK0s) delete fListK0s;
  if(fListLa) delete fListLa;
  if(fListALa) delete fListALa;
  if(fListFeeddownLaCand) delete fListFeeddownLaCand;
  if(fListFeeddownALaCand) delete fListFeeddownALaCand;
  if(jetConeFDLalist) delete jetConeFDLalist;
  if(jetConeFDALalist) delete jetConeFDALalist;   
  if(fListMCgenK0s) delete fListMCgenK0s;
  if(fListMCgenLa) delete fListMCgenLa;
  if(fListMCgenALa) delete fListMCgenALa;



}

//________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AliFragFuncHistosInvMass(const char* name, 
									   Int_t nJetPt, Float_t jetPtMin, Float_t jetPtMax,  
									   Int_t nInvMass, Float_t invMassMin, Float_t invMassMax,
									   Int_t nPt, Float_t ptMin, Float_t ptMax,
									   Int_t nXi, Float_t xiMin, Float_t xiMax,
									   Int_t nZ , Float_t zMin , Float_t zMax )
  : TObject()
  ,fNBinsJetPt(nJetPt)
  ,fJetPtMin(jetPtMin)
  ,fJetPtMax(jetPtMax)
  ,fNBinsInvMass(nInvMass)
  ,fInvMassMin(invMassMin)  
  ,fInvMassMax(invMassMax)
  ,fNBinsPt(nPt) 
  ,fPtMin(ptMin)   
  ,fPtMax(ptMax)   
  ,fNBinsXi(nXi) 
  ,fXiMin(xiMin)   
  ,fXiMax(xiMax)   
  ,fNBinsZ(nZ)  
  ,fZMin(zMin)    
  ,fZMax(zMax)    
  ,fh3TrackPt(0)
  ,fh3Xi(0)
  ,fh3Z(0)
  ,fh1JetPt(0)
  ,fNameFF(name)
{
  // default constructor

}

//______________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AliFragFuncHistosInvMass(const AliFragFuncHistosInvMass& copy)
  : TObject()
  ,fNBinsJetPt(copy.fNBinsJetPt)
  ,fJetPtMin(copy.fJetPtMin)
  ,fJetPtMax(copy.fJetPtMax)
  ,fNBinsInvMass(copy.fNBinsInvMass)
  ,fInvMassMin(copy.fInvMassMin)  
  ,fInvMassMax(copy.fInvMassMax)
  ,fNBinsPt(copy.fNBinsPt) 
  ,fPtMin(copy.fPtMin)   

  ,fPtMax(copy.fPtMax)   
  ,fNBinsXi(copy.fNBinsXi) 
  ,fXiMin(copy.fXiMin)   
  ,fXiMax(copy.fXiMax)   
  ,fNBinsZ(copy.fNBinsZ)  
  ,fZMin(copy.fZMin)    
  ,fZMax(copy.fZMax)    
  ,fh3TrackPt(copy.fh3TrackPt)
  ,fh3Xi(copy.fh3Xi)
  ,fh3Z(copy.fh3Z)
  ,fh1JetPt(copy.fh1JetPt)
  ,fNameFF(copy.fNameFF)
{
  // copy constructor
}

//______________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass& AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::operator=(const AliAnalysisTaskJetChem::AliFragFuncHistosInvMass& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsJetPt   = o.fNBinsJetPt;
    fJetPtMin     = o.fJetPtMin;
    fJetPtMax     = o.fJetPtMax;
    fNBinsInvMass = o.fNBinsInvMass;
    fInvMassMin   = o.fInvMassMin;  
    fInvMassMax   = o.fInvMassMax;
    fNBinsPt      = o.fNBinsPt; 
    fPtMin        = o.fPtMin;   
    fPtMax        = o.fPtMax;   
    fNBinsXi      = o.fNBinsXi; 
    fXiMin        = o.fXiMin;   
    fXiMax        = o.fXiMax;   
    fNBinsZ       = o.fNBinsZ;  
    fZMin         = o.fZMin;    
    fZMax         = o.fZMax;    
    fh3TrackPt    = o.fh3TrackPt;
    fh3Xi         = o.fh3Xi;
    fh3Z          = o.fh3Z;
    fh1JetPt      = o.fh1JetPt;
    fNameFF       = o.fNameFF;
  }
    
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::~AliFragFuncHistosInvMass()
{ 
  // destructor 

  if(fh1JetPt)   delete fh1JetPt;
  if(fh3TrackPt) delete fh3TrackPt;
  if(fh3Xi)      delete fh3Xi;
  if(fh3Z)       delete fh3Z;
}

//_________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::DefineHistos()
{
  // book FF histos

  fh1JetPt   = new TH1F(Form("fh1FFJetPtIM%s", fNameFF.Data()),"",fNBinsJetPt,fJetPtMin,fJetPtMax);
  fh3TrackPt = new TH3F(Form("fh3FFTrackPtIM%s",fNameFF.Data()),"",fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsPt, fPtMin, fPtMax);
  fh3Xi      = new TH3F(Form("fh3FFXiIM%s", fNameFF.Data()),"", fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsXi, fXiMin, fXiMax);
  fh3Z       = new TH3F(Form("fh3FFZIM%s", fNameFF.Data()),"", fNBinsJetPt, fJetPtMin, fJetPtMax, fNBinsInvMass, fInvMassMin, fInvMassMax, fNBinsZ, fZMin, fZMax);

  AliAnalysisTaskFragmentationFunction::SetProperties(fh1JetPt, "p_{t} (GeV/c)", "entries"); 
  AliAnalysisTaskJetChem::SetProperties(fh3TrackPt,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","p_{t} (GeV/c)");
  AliAnalysisTaskJetChem::SetProperties(fh3Xi,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","#xi");
  AliAnalysisTaskJetChem::SetProperties(fh3Z,"jet p_{t} (GeV/c)","inv Mass (GeV/c^2)","z");
}

//________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::FillFF(Float_t trackPt, Float_t invM, Float_t jetPt, Bool_t incrementJetPt)
{
  // fill FF
 
  if(incrementJetPt) fh1JetPt->Fill(jetPt);    
  fh3TrackPt->Fill(jetPt,invM,trackPt);//Fill(x,y,z)
  
  Double_t z = 0.;
  if(jetPt>0) z = trackPt / jetPt;
  Double_t xi = 0;
  if(z>0) xi = TMath::Log(1/z);
  
  fh3Xi->Fill(jetPt,invM,xi);
  fh3Z->Fill(jetPt,invM,z);
}

//___________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  list->Add(fh3TrackPt);
  list->Add(fh3Xi);
  list->Add(fh3Z);
}

// ---


//_______________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AliFragFuncHistosPhiCorrInvMass(const char* name,
											 Int_t nPt,   Float_t ptMin,   Float_t ptMax,
											 Int_t nPhi,  Float_t phiMin,  Float_t phiMax,
											 Int_t nInvMass,  Float_t invMassMin,  Float_t invMassMax)
  : TObject()
  ,fNBinsPt(nPt)
  ,fPtMin(ptMin)
  ,fPtMax(ptMax)
  ,fNBinsPhi(nPhi)
  ,fPhiMin(phiMin)
  ,fPhiMax(phiMax)
  ,fNBinsInvMass(nInvMass)
  ,fInvMassMin(invMassMin)  
  ,fInvMassMax(invMassMax)
  ,fh3PhiCorr(0) 
  ,fNamePhiCorr(name) 
{
  // default constructor
}

//____________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AliFragFuncHistosPhiCorrInvMass(const AliFragFuncHistosPhiCorrInvMass& copy)
  : TObject()
  ,fNBinsPt(copy.fNBinsPt)
  ,fPtMin(copy.fPtMin)
  ,fPtMax(copy.fPtMax)
  ,fNBinsPhi(copy.fNBinsPhi)
  ,fPhiMin(copy.fPhiMin)
  ,fPhiMax(copy.fPhiMax)
  ,fNBinsInvMass(copy.fNBinsInvMass)
  ,fInvMassMin(copy.fInvMassMin)  
  ,fInvMassMax(copy.fInvMassMax)
  ,fh3PhiCorr(copy.fh3PhiCorr)
  ,fNamePhiCorr(copy.fNamePhiCorr)
{
  // copy constructor
}

//________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass& AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::operator=(const AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fNBinsPt      = o.fNBinsPt;
    fPtMin        = o.fPtMin;
    fPtMax        = o.fPtMax;
    fNBinsPhi     = o.fNBinsPhi;
    fPhiMin       = o.fPhiMin;
    fPhiMax       = o.fPhiMax;
    fNBinsInvMass = o.fNBinsInvMass;
    fInvMassMin   = o.fInvMassMin;  
    fInvMassMax   = o.fInvMassMax;
    fh3PhiCorr    = o.fh3PhiCorr;
    fNamePhiCorr  = o.fNamePhiCorr;
  }
  
  return *this;
}

//_________________________________________________________________________________________
AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::~AliFragFuncHistosPhiCorrInvMass()
{
  // destructor 
  
  if(fh3PhiCorr) delete fh3PhiCorr;
}

//__________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::DefineHistos()
{
  // book jet QA histos

  fh3PhiCorr  = new TH3F(Form("fh3PhiCorrIM%s", fNamePhiCorr.Data()), 
			 Form("%s: p_{t} - #phi - m_{inv} distribution",fNamePhiCorr.Data()), 
			 fNBinsPt, fPtMin, fPtMax, 
			 fNBinsPhi, fPhiMin, fPhiMax,
			 fNBinsInvMass, fInvMassMin, fInvMassMax);
  
  AliAnalysisTaskJetChem::SetProperties(fh3PhiCorr, "p_{t} (GeV/c)", "#phi", "m_{inv} (GeV/c^2)"); 
}

//___________________________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::FillPhiCorr(Float_t pt, Float_t phi, Float_t invM)
{
  // fill jet QA histos 

  fh3PhiCorr->Fill(pt, phi, invM);
}

//______________________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosPhiCorrInvMass::AddToOutput(TList* list) const 
{
  // add histos to list

  list->Add(fh3PhiCorr);
}

//____________________________________________________
void AliAnalysisTaskJetChem::UserCreateOutputObjects()
{
  // create output objects
   
  if(fDebug > 1) Printf("AliAnalysisTaskJetChem::UserCreateOutputObjects()");
 
  // create list of tracks and jets 
  
  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE); //objects in TList wont be deleted when TList is deleted 
  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);
  fBckgJetsRec = new TList();
  fBckgJetsRec->SetOwner(kFALSE);
  fListK0s = new TList(); 
  fListK0s->SetOwner(kFALSE);
  fListLa = new TList(); 
  fListLa->SetOwner(kFALSE);
  fListALa = new TList(); 
  fListALa->SetOwner(kFALSE);
  fListFeeddownLaCand = new TList();    //feeddown Lambda candidates
  fListFeeddownLaCand->SetOwner(kFALSE);
  fListFeeddownALaCand = new TList();   //feeddown Antilambda candidates
  fListFeeddownALaCand->SetOwner(kFALSE);
  jetConeFDLalist = new TList();     
  jetConeFDLalist->SetOwner(kFALSE);  //feeddown Lambda candidates in jet cone
  jetConeFDALalist = new TList();     
  jetConeFDALalist->SetOwner(kFALSE); //feeddown Antilambda candidates in jet cone
  fListMCgenK0s = new TList();          //MC generated K0s 
  fListMCgenK0s->SetOwner(kFALSE);
  fListMCgenLa = new TList();           //MC generated Lambdas
  fListMCgenLa->SetOwner(kFALSE);
  fListMCgenALa = new TList();          //MC generated Antilambdas
  fListMCgenALa->SetOwner(kFALSE);

  
  // Create histograms / output container
 
  //for AliPIDResponse:
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  OpenFile(1);
  fCommonHistList = new TList();
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);//By default (fAddDirectory = kTRUE), histograms are automatically added to the list of objects in memory
  	
  // histograms inherited from AliAnalysisTaskFragmentationFunction

  fh1EvtSelection            = new TH1F("fh1EvtSelection", "Event Selection", 6, -0.5, 5.5);
  fh1EvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
  fh1EvtSelection->GetXaxis()->SetBinLabel(2,"event trigger selection: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(3,"event class: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(4,"vertex Ncontr: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(5,"vertex z: rejected");
  fh1EvtSelection->GetXaxis()->SetBinLabel(6,"vertex type: rejected");


  fh1EvtCent 	             = new TH1F("fh1EvtCent","centrality",100,0.,100.);
  fh1VertexNContributors     = new TH1F("fh1VertexNContributors", "Vertex N contributors", 11,-.5, 10.5);
  fh1VertexZ                 = new TH1F("fh1VertexZ", "Vertex z distribution", 30, -15., 15.);
  fh1Xsec                    = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Trials                  = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1PtHard                  = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fh1PtHardTrials            = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);
  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",10,-0.5,9.5);
 
  // histograms JetChem task
 
  fh1EvtAllCent	                = new TH1F("fh1EvtAllCent","before centrality selection",100,0.,100.);
  fh1Evt                        = new TH1F("fh1Evt", "All events runned over", 3, 0.,1.);
  fh1EvtMult 	                = new TH1F("fh1EvtMult","multiplicity",1200,0.,12000.);
  fh1K0Mult 	                = new TH1F("fh1K0Mult","K0 multiplicity",1000,0.,1000.);//500. all
  fh1dPhiJetK0                  = new TH1F("fh1dPhiJetK0","",640,-1,5.4);
  fh1LaMult 	                = new TH1F("fh1LaMult","La multiplicity",1000,0.,1000.);
  fh1dPhiJetLa                  = new TH1F("fh1dPhiJetLa","",640,-1,5.4);
  fh1ALaMult 	                = new TH1F("fh1ALaMult","ALa multiplicity",1000,0.,1000.);
  fh1dPhiJetALa                 = new TH1F("fh1dPhiJetALa","",640,-1,5.4);
  fh1JetEta                     = new TH1F("fh1JetEta","#eta distribution of all jets",400,-2.,2.);
  fh1JetPhi                     = new TH1F("fh1JetPhi","#phi distribution of all jets",630,0.,6.3);
  fh2JetEtaPhi                  = new TH2F("fh2JetEtaPhi","#eta and #phi distribution of all jets",400,-2.,2.,630,0.,6.3);
  fh1V0JetPt                    = new TH1F("fh1V0JetPt","#p_{T} distribution of all jets containing v0s",200,0.,200.);
  fh2FFJetTrackEta              = new TH2F("fh2FFJetTrackEta","charged track eta distr. in jet cone",200,-1.,1.,40,0.,200.);  
  fh1trackPosNCls               = new TH1F("fh1trackPosNCls","NTPC clusters positive daughters",250,0.,250.);
  fh1trackNegNCls               = new TH1F("fh1trackNegNCls","NTPC clusters negative daughters",250,0.,250.);
  fh1trackPosEta                = new TH1F("fh1trackPosEta","eta positive daughters",100,-2.,2.);
  fh1trackNegEta                = new TH1F("fh1trackNegEta","eta negative daughters",100,-2.,2.);
  fh1V0Eta                      = new TH1F("fh1V0Eta","V0 eta",60,-1.5,1.5);
  fh1V0totMom                   = new TH1F("fh1V0totMom","V0 tot mom",240,0.,20.); 
  fh1CosPointAngle              = new TH1F("fh1CosPointAngle", "Cosine of V0's pointing angle",1000,0.99,1.0);
  fh1Chi2Pos                    = new TH1F("fh1Chi2Pos", "V0s chi2",100,0.,5.);
  fh1Chi2Neg                    = new TH1F("fh1Chi2Neg", "V0s chi2",100,0.,5.);
  fh1DecayLengthV0              = new TH1F("fh1DecayLengthV0", "V0s decay Length;decay length(cm)",1200,0.,120.);
  fh2ProperLifetimeK0sVsPtBeforeCut = new TH2F("fh2ProperLifetimeK0sVsPtBeforeCut"," K0s ProperLifetime vs Pt; p_{T} (GeV/#it{c})",1500,0.,15.,500,0.,250.);
  fh2ProperLifetimeK0sVsPtAfterCut = new TH2F("fh2ProperLifetimeK0sVsPtAfterCut"," K0s ProperLifetime vs Pt; p_{T} (GeV/#it{c})",1500,0.,15.,500,0.,250.);
  fh1ProperLifetimeV0BeforeCut  = new TH1F("fh1ProperLifetimeV0BeforeCut", "V0s 2D distance over transerse mom.;(cm)",1200,0.,120.);
  fh1ProperLifetimeV0AfterCut   = new TH1F("fh1ProperLifetimeV0AfterCut", "V0s 2D distance over transverse mom.;(cm)",1200,0.,120.);
  fh1V0Radius                   = new TH1F("fh1V0Radius", "V0s Radius;Radius(cm)",400,0.,40.);
  fh1DcaV0Daughters             = new TH1F("fh1DcaV0Daughters", "DCA between daughters;dca(cm)",200,0.,2.);
  fh1DcaPosToPrimVertex         = new TH1F("fh1DcaPosToPrimVertex", "Positive V0 daughter;dca(cm)",1000,0.,10.);
  fh1DcaNegToPrimVertex         = new TH1F("fh1DcaNegToPrimVertex", "Negative V0 daughter;dca(cm)",1000,0.,10.);
  fh2ArmenterosBeforeCuts       = new TH2F("fh2ArmenterosBeforeCuts","Armenteros Podolanski Plot for K0s Candidates;#alpha;(p^{arm})_{T}/(GeV/#it{c})",200,-1.2,1.2,600,0.,0.35);
  fh2ArmenterosAfterCuts        = new TH2F("fh2ArmenterosAfterCuts","Armenteros Podolanski Plot for K0s Candidates;#alpha;(p^{arm})_{T}/(GeV/#it{c});",200,-1.2,1.2,600,0.,0.35);
  fh2BB3SigProton               = new TH2F("fh2BB3SigProton","-dE/dX against Momentum for Protons @3sigma from TPC; P (GeV); -dE/dx (keV/cm ?)",1000,0.,10.,1000,0.,200.);
  fh2BBLaPos                    = new TH2F("fh2BBLaPos","PID of the positive daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
  fh2BBLaNeg                    = new TH2F("fh2BBLaNeg","PID of the negative daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
  fh1CrossedRowsOverFindableNeg = new TH1F("fh1CrossedRowsOverFindableNeg","pos daughter crossed rows over findable in TPC;counts",200,0.,2.);
  fh1CrossedRowsOverFindablePos = new TH1F("fh1CrossedRowsOverFindablePos","neg daughter crossed rows over findable in TPC;counts",200,0.,2.);
  fh1PosDaughterCharge          = new TH1F("fh1PosDaughterCharge","charge of V0 positive daughters; V0 daughters",3,-2.,2.);
  fh1NegDaughterCharge          = new TH1F("fh1NegDaughterCharge","charge of V0 negative daughters; V0 daughters",3,-2.,2.);
  fh1PtMCK0s                    = new TH1F("fh1PtMCK0s","Pt of MC rec K0s; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1PtMCLa                     = new TH1F("fh1PtMCLa","Pt of MC rec La; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1PtMCALa                    = new TH1F("fh1PtMCALa","Pt of MC rec ALa; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1EtaK0s                     = new TH1F("fh1EtaK0s","K^{0}_{s} entries ;#eta",200,-1.,1.);
  fh1EtaLa                      = new TH1F("fh1EtaLa","#Lambda entries ;#eta",200,-1.,1.);
  fh1EtaALa                     = new TH1F("fh1EtaALa","#bar{#Lambda} entries ;#eta",200,-1.,1.);
  fh3InvMassEtaTrackPtK0s       = new TH3F("fh3InvMassEtaTrackPtK0s","#eta; invMass (GeV/{#it{c}}^{2}); #it{p}_{T} (GeV/#it{c})", 200, -1., 1., 240, 0.4, 0.6, 140, 0., 14.);
  fh3InvMassEtaTrackPtLa        = new TH3F("fh3InvMassEtaTrackPtLa", "#eta; invMass (GeV/{#it{c}}^{2}; #it{p}_{T} (GeV/#it{c}))",  200, -1., 1., 140, 1.06, 1.2, 140, 0., 14.);
  fh3InvMassEtaTrackPtALa       = new TH3F("fh3InvMassEtaTrackPtALa","#eta; invMass (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})",  200, -1., 1., 140, 1.06, 1.2, 140, 0., 14.);
  fh3IMK0PerpCone               = new TH3F("fh3IMK0PerpCone","{K_{0}}^{s} content in perpendicular cone",39,5.,200., 400,0.3,0.7, 200,0.,20.);
  fh3IMLaPerpCone               = new TH3F("fh3IMLaPerpCone","#Lambda content in perpendicular cone",39,5.,200., 140,1.06,1.2, 200,0.,20.);
  fh3IMALaPerpCone              = new TH3F("fh3IMALaPerpCone","#Antilambda content in perpendicular cone",39,5.,200., 140,1.06,1.2, 200,0.,20.);
  fh3IMK0MedianCone             = new TH3F("fh3IMK0MedianCone","{K_{0}}^{s} content in median cluster cone",39,5.,200., 400,0.3,0.7, 200,0.,20.);
  fh3IMLaMedianCone             = new TH3F("fh3IMLaMedianCone","#Lambda content in median cluster cone",39,5.,200., 140,1.06,1.2, 200,0.,20.);
  fh3IMALaMedianCone            = new TH3F("fh3IMALaMedianCone","#Antilambda content in median cluster cone",39,5.,200., 140,1.06,1.2, 200,0.,20.);

  fh1noAssociatedK0s            = new TH1F("fh1noAssociatedK0s","not selected as associated particle",12,0.,12.);
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(1,"K0s: accepted as associated particle");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(2,"K0s: v0 not K0s pdg code (310)");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(3,"K0s: v0 is not primary particle");

  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(4,"K0s: pos daughter is pion");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(5,"K0s: neg daughter is pion");

  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(6,"K0s: pos daughter particle is proton");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(7,"K0s: pos daughter particle is electron"); 
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(8,"K0s: pos daughter particle is myon");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(9,"K0s: neg daughter particle is proton");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(10,"K0s: neg daughter particle is electron");  
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(11,"K0s: neg daughter particle is myon");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(12,"K0s: pos daughter particle is something else");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(13,"K0s: neg daughter particle is something else");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(14,"K0s: pos daughter not pion pdg code (211)");
  fh1noAssociatedK0s->GetXaxis()->SetBinLabel(15,"K0s: neg daughter not pion pdg code (211)");

  fh1TrackMultCone          = new TH1F("fh1TrackMultCone","track multiplicity in jet cone; number of tracks",200,0.,1000.);

  fh2TrackMultCone          = new TH2F("fh2TrackMultCone","track multiplicity in jet cone vs. jet momentum; number of tracks; jet it{p}_{T} (GeV/it{c})",200,0.,1000.,39,5.,200.);

  fFFHistosRecCuts   	    = new AliFragFuncHistos("RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  fV0QAK0                   = new AliFragFuncQATrackHistos("V0QAK0",fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							    fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							    fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							    fQATrackHighPtThreshold);
  
  fFFHistosRecCutsK0Evt      = new AliFragFuncHistos("RecCutsK0Evt", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  
  fFFHistosIMK0AllEvt        = new AliFragFuncHistosInvMass("K0AllEvt", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
  
  fFFHistosIMK0Jet           = new AliFragFuncHistosInvMass("K0Jet", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
    
  fFFHistosIMK0Cone          = new AliFragFuncHistosInvMass("K0Cone", fFFIMNBinsJetPt, fFFIMJetPtMin, fFFIMJetPtMax, 
							    fFFIMNBinsInvM,fFFIMInvMMin,fFFIMInvMMax,
							    fFFIMNBinsPt, fFFIMPtMin, fFFIMPtMax, 
							    fFFIMNBinsXi, fFFIMXiMin, fFFIMXiMax,  
							    fFFIMNBinsZ , fFFIMZMin , fFFIMZMax);
  
  fFFHistosPhiCorrIMK0       = new AliFragFuncHistosPhiCorrInvMass("K0",fPhiCorrIMNBinsPt, fPhiCorrIMPtMin, fPhiCorrIMPtMax, 
								   fPhiCorrIMNBinsPhi, fPhiCorrIMPhiMin, fPhiCorrIMPhiMax,  
								   fPhiCorrIMNBinsInvM , fPhiCorrIMInvMMin , fPhiCorrIMInvMMax);
  
  fFFHistosIMLaAllEvt        = new AliFragFuncHistosInvMass("LaAllEvt", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  fFFHistosIMLaJet           = new AliFragFuncHistosInvMass("LaJet", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  
  fFFHistosIMLaCone          = new AliFragFuncHistosInvMass("LaCone", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  fFFHistosPhiCorrIMLa       = new AliFragFuncHistosPhiCorrInvMass("La",fPhiCorrIMLaNBinsPt, fPhiCorrIMLaPtMin, fPhiCorrIMLaPtMax, 
								   fPhiCorrIMLaNBinsPhi, fPhiCorrIMLaPhiMin, fPhiCorrIMLaPhiMax,  
								   fPhiCorrIMLaNBinsInvM , fPhiCorrIMLaInvMMin , fPhiCorrIMLaInvMMax);

 
  fFFHistosIMALaAllEvt        = new AliFragFuncHistosInvMass("ALaAllEvt", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  fFFHistosIMALaJet           = new AliFragFuncHistosInvMass("ALaJet", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  fFFHistosIMALaCone          = new AliFragFuncHistosInvMass("ALaCone", fFFIMLaNBinsJetPt, fFFIMLaJetPtMin, fFFIMLaJetPtMax, 
							    fFFIMLaNBinsInvM,fFFIMLaInvMMin,fFFIMLaInvMMax,
							    fFFIMLaNBinsPt, fFFIMLaPtMin, fFFIMLaPtMax, 
							    fFFIMLaNBinsXi, fFFIMLaXiMin, fFFIMLaXiMax,  
							    fFFIMLaNBinsZ , fFFIMLaZMin , fFFIMLaZMax);
  
  fFFHistosPhiCorrIMALa       = new AliFragFuncHistosPhiCorrInvMass("ALa",fPhiCorrIMLaNBinsPt, fPhiCorrIMLaPtMin, fPhiCorrIMLaPtMax, 
								   fPhiCorrIMLaNBinsPhi, fPhiCorrIMLaPhiMin, fPhiCorrIMLaPhiMax,  
								   fPhiCorrIMLaNBinsInvM , fPhiCorrIMLaInvMMin , fPhiCorrIMLaInvMMax);

  //***************
  // MC histograms
  //***************

  fh2MCgenK0Cone                = new TH2F("fh2MCgenK0Cone", "MC gen {K^{0}}^{s} #it{p}_{T}  in cone around rec jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",39,5.,200.,200,0.,20.);
  fh2MCgenLaCone                = new TH2F("fh2MCgenLaCone", "MC gen #Lambda #it{p}_{T} in cone around rec jet axis versus jet #it{p}_{T} ; jet #it{p}_{T}",39,5.,200.,200,0.,20.);
  fh2MCgenALaCone               = new TH2F("fh2MCgenALaCone", "MC gen #Antilambda #it{p}_{T} in cone around rec jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",39,5.,200.,200,0.,20.);

  fh2MCgenK0Cone->GetYaxis()->SetTitle("MC gen K^{0}}^{s} #it{p}_{T}");
  fh2MCgenLaCone->GetYaxis()->SetTitle("MC gen #Lambda #it{p}_{T}");
  fh2MCgenALaCone->GetYaxis()->SetTitle("MC gen #Antilambda #it{p}_{T}");

  fh2MCEtagenK0Cone             = new TH2F("fh2MCEtagenK0Cone","MC gen {K^{0}}^{s} #it{p}_{T} #eta distribution in jet cone;#eta",39,5.,200.,200,-1.,1.);
  fh2MCEtagenLaCone             = new TH2F("fh2MCEtagenLaCone","MC gen #Lambda #it{p}_{T} #eta distribution in jet cone;#eta",39,5.,200.,200,-1.,1.);
  fh2MCEtagenALaCone            = new TH2F("fh2MCEtagenALaCone","MC gen #Antilambda #it{p}_{T} #eta distribution in jet cone;#eta",39,5.,200.,200,-1.,1.);
  fh1FFIMK0ConeSmear            = new TH1F("fh1FFIMK0ConeSmear","Smeared jet pt study for K0s-in-cone-jets; smeared jet #it{p}_{T}", 39,5.,200.);
  fh1FFIMLaConeSmear            = new TH1F("fh1FFIMLaConeSmear","Smeared jet pt study for La-in-cone-jets; smeared jet #it{p}_{T}", 39,5.,200.);
  fh1FFIMALaConeSmear           = new TH1F("fh1FFIMALaConeSmear","Smeared jet pt study for ALa-in-cone-jets; smeared jet #it{p}_{T}", 39,5.,200.);
  
  fh3MCrecK0Cone                = new TH3F("fh3MCrecK0Cone", "MC rec {K^{0}}^{s} #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",39,5.,200., 400,0.3,0.7, 200,0.,20.);  
  fh3MCrecLaCone                = new TH3F("fh3MCrecLaCone", "MC rec {#Lambda #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2});#it{p}_{T}",39,5.,200., 140,1.06,1.2, 200,0.,20.);            
  fh3MCrecALaCone               = new TH3F("fh3MCrecALaCone", "MC rec {#Antilambda #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2});#it{p}_{T}",39,5.,200.,140,1.06,1.2, 200,0.,20.);
  fh3MCrecK0ConeSmear           = new TH3F("fh3MCrecK0ConeSmear", "MC rec {K^{0}}^{s} #it{p}_{T}  in cone around jet axis matching MC gen particle, with jet p_{T} smeared; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",39,5.,200., 400,0.3,0.7, 200,0.,20.);  
  fh3MCrecLaConeSmear           = new TH3F("fh3MCrecLaConeSmear", "MC rec {#Lambda #it{p}_{T}  in cone around jet axis matching MC gen particle, with jet p_{T} smeared; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2});#it{p}_{T}",39,5.,200., 140,1.06,1.2, 200,0.,20.);            
  fh3MCrecALaConeSmear          = new TH3F("fh3MCrecALaConeSmear", "MC rec {#Antilambda #it{p}_{T}  in cone around jet axis matching MC gen particle, with jet p_{T} smeared; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2});#it{p}_{T}",39,5.,200.,140,1.06,1.2, 200,0.,20.);
  fh3SecContinCone              = new TH3F("fh3SecContinCone","secondary contamination of jet cones; jet #it{p}_{T}; track #it{p}_{T}, #eta",39,5.,200.,200,0.,20.,200,-1.,1.);
  fh3StrContinCone              = new TH3F("fh3StrContinCone","strange particle contamination of jet cones; jet #it{p}_{T}; track #it{p}_{T}, #eta",39,5.,200.,200,0.,20.,200,-1.,1.);

  fh1MCMultiplicityPrimary      = new TH1F("fh1MCMultiplicityPrimary", "MC Primary Particles;NPrimary;Count", 201, -0.5, 200.5);
  fh1MCMultiplicityTracks       = new TH1F("h1MCMultiplicityTracks", "MC Tracks;Ntracks;Count", 201, -0.5, 200.5);
  // fh1MCmotherK0s             = new TH1F("fh1MCmotherK0s","K0s mother pdg codes",10,0.,10.);
  fh1MCmotherLa                 = new TH1F("fh1MCmotherLa","Lambdas mother pdg codes",10,0.,10.);
  fh1MCmotherLa->GetXaxis()->SetBinLabel(1.,"#Sigma^{-}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(2.,"#Sigma^{0}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(3.,"#Sigma^{+}");  
  fh1MCmotherLa->GetXaxis()->SetBinLabel(4.,"#Omega^{-}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(5.,"#Xi^{0}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(6.,"#Xi^{-}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(7.,"#Xi^{+}");
  fh1MCmotherLa->GetXaxis()->SetBinLabel(8.,"primary particle");
  fh1MCmotherALa                = new TH1F("fh1MCmotherALa","Antilambdas mother pdg codes",10,0.,10.);
  fh1MCmotherALa->GetXaxis()->SetBinLabel(1.,"#bar{#Sigma^{-}}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(2.,"#bar{#Sigma^{0}}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(3.,"#bar{#Sigma^{+}}");  
  fh1MCmotherALa->GetXaxis()->SetBinLabel(4.,"#bar{#Omega^{-}}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(5.,"#bar{#Xi^{0}}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(6.,"#Xi^{-}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(7.,"#Xi^{+}");
  fh1MCmotherALa->GetXaxis()->SetBinLabel(8.,"primary particle");
  fh3FeedDownLa                 = new TH3F("fh3FeedDownLa","#Lambda stemming from feeddown from Xi(0/-)", 39, 5., 200., 200, 1.05, 1.25, 200,0.,20.);
  fh3FeedDownALa                = new TH3F("fh3FeedDownALa","#bar#Lambda stemming from feeddown from Xibar(0/+)", 39, 5., 200., 200, 1.05, 1.25, 200, 0., 20.);
  fh1MCProdRadiusK0s            = new TH1F("fh1MCProdRadiusK0s","MC gen. MC K0s prod radius",600,0.,200.);
  fh1MCProdRadiusLambda         = new TH1F("fh1MCProdRadiusLambda","MC gen. MC La prod radius",600,0.,200.);
  fh1MCProdRadiusAntiLambda     = new TH1F("fh1MCProdRadiusAntiLambda","MC gen. MC ALa prod radius",600,0.,200.);

  // Pt and inv mass distributions

  fh1MCPtV0s                    = new TH1F("fh1MCPtV0s", "MC gen. V^{0} in rap range;#it{p}_{T} (GeV/#it{c})",200,0,20.);// 0.1 GeV/c steps
  fh1MCPtK0s                    = new TH1F("fh1MCPtK0s", "MC gen. K^{0}_{s} in eta range;#it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1MCPtLambda                 = new TH1F("fh1MCPtLambda", "MC gen. #Lambda in rap range;#it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1MCPtAntiLambda             = new TH1F("fh1MCPtAntiLambda", "MC gen. #AntiLambda in rap range;#it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1MCXiPt                     = new TH1F("fh1MCXiPt", "MC gen. #Xi^{-/o};#it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1MCXibarPt                  = new TH1F("fh1MCXibarPt", "MC gen. #bar{#Xi}^{+/o};#it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh2MCEtaVsPtK0s               = new TH2F("fh2MCEtaVsPtK0s","MC gen. K^{0}_{s} #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);
  fh2MCEtaVsPtLa                = new TH2F("fh2MCEtaVsPtLa","MC gen. #Lambda #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);
  fh2MCEtaVsPtALa               = new TH2F("fh2MCEtaVsPtALa","MC gen. #bar{#Lambda}  #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);

  // Rapidity
  fh1MCRapK0s                   = new TH1F("fh1MCRapK0s", "MC gen. K0s;rap with cut",200,-10,10); 
  fh1MCRapLambda                = new TH1F("fh1MCRapLambda", "MC gen. #Lambda;rap",200,-10,10);
  fh1MCRapAntiLambda            = new TH1F("fh1MCRapAntiLambda", "MC gen. #bar{#Lambda};rap",200,-10,10);
  fh1MCEtaAllK0s                = new TH1F("fh1MCEtaAllK0s", "MC gen. K0s;#eta",200,-1.,1.); 
  fh1MCEtaK0s                   = new TH1F("fh1MCEtaK0s", "MC gen. K0s;#eta with cut",200,-1.,1.); 
  fh1MCEtaLambda                = new TH1F("fh1MCEtaLambda", "MC gen. #Lambda;#eta",200,-1.,1.);
  fh1MCEtaAntiLambda            = new TH1F("fh1MCEtaAntiLambda", "MC gen. #bar{#Lambda};#eta",200,-1.,1.);

  fV0QAK0->DefineHistos();
  fFFHistosRecCuts->DefineHistos();
  fFFHistosRecCutsK0Evt->DefineHistos();
  fFFHistosIMK0AllEvt->DefineHistos();
  fFFHistosIMK0Jet->DefineHistos();
  fFFHistosIMK0Cone->DefineHistos();
  fFFHistosPhiCorrIMK0->DefineHistos();
  fFFHistosIMLaAllEvt->DefineHistos();
  fFFHistosIMLaJet->DefineHistos();
  fFFHistosIMLaCone->DefineHistos();
  fFFHistosPhiCorrIMLa->DefineHistos();
  fFFHistosIMALaAllEvt->DefineHistos();
  fFFHistosIMALaJet->DefineHistos();
  fFFHistosIMALaCone->DefineHistos();
  fFFHistosPhiCorrIMALa->DefineHistos();

  const Int_t saveLevel = 5;
  if(saveLevel>0){

    fCommonHistList->Add(fh1EvtAllCent);
    fCommonHistList->Add(fh1Evt);
    fCommonHistList->Add(fh1EvtSelection);
    fCommonHistList->Add(fh1EvtCent);
    fCommonHistList->Add(fh1VertexNContributors);
    fCommonHistList->Add(fh1VertexZ);
    fCommonHistList->Add(fh1Xsec);
    fCommonHistList->Add(fh1Trials);
    fCommonHistList->Add(fh1PtHard);
    fCommonHistList->Add(fh1PtHardTrials);
    fCommonHistList->Add(fh1nRecJetsCuts);
    fCommonHistList->Add(fh1EvtMult);
    fCommonHistList->Add(fh1K0Mult);
    fCommonHistList->Add(fh1dPhiJetK0);
    fCommonHistList->Add(fh1LaMult);
    fCommonHistList->Add(fh1dPhiJetLa);
    fCommonHistList->Add(fh1ALaMult);
    fCommonHistList->Add(fh1dPhiJetALa);
    fCommonHistList->Add(fh1JetEta);        
    fCommonHistList->Add(fh1JetPhi);               
    fCommonHistList->Add(fh2JetEtaPhi);
    fCommonHistList->Add(fh1V0JetPt); 
    fCommonHistList->Add(fh2FFJetTrackEta);   
    fCommonHistList->Add(fh1trackPosNCls);           
    fCommonHistList->Add(fh1trackNegNCls);          
    fCommonHistList->Add(fh1trackPosEta);            
    fCommonHistList->Add(fh1trackNegEta);          
    fCommonHistList->Add(fh1V0Eta); 
    fCommonHistList->Add(fh1V0totMom);        
    fCommonHistList->Add(fh1CosPointAngle);        
    fCommonHistList->Add(fh1Chi2Pos);                 
    fCommonHistList->Add(fh1Chi2Neg);              
    fCommonHistList->Add(fh1DecayLengthV0); 
    fCommonHistList->Add(fh2ProperLifetimeK0sVsPtBeforeCut);
    fCommonHistList->Add(fh2ProperLifetimeK0sVsPtAfterCut);
    fCommonHistList->Add(fh1ProperLifetimeV0BeforeCut);
    fCommonHistList->Add(fh1ProperLifetimeV0AfterCut);
    fCommonHistList->Add(fh1V0Radius);     
    fCommonHistList->Add(fh1DcaV0Daughters);        
    fCommonHistList->Add(fh1DcaPosToPrimVertex);   
    fCommonHistList->Add(fh1DcaNegToPrimVertex);    
    fCommonHistList->Add(fh2ArmenterosBeforeCuts);
    fCommonHistList->Add(fh2ArmenterosAfterCuts);
    fCommonHistList->Add(fh2BB3SigProton);
    fCommonHistList->Add(fh2BBLaPos);
    fCommonHistList->Add(fh2BBLaNeg);
    fCommonHistList->Add(fh1CrossedRowsOverFindableNeg);
    fCommonHistList->Add(fh1CrossedRowsOverFindablePos);
    fCommonHistList->Add(fh1PosDaughterCharge);
    fCommonHistList->Add(fh1NegDaughterCharge);
    fCommonHistList->Add(fh1PtMCK0s);
    fCommonHistList->Add(fh1PtMCLa);
    fCommonHistList->Add(fh1PtMCALa);
    fCommonHistList->Add(fh1EtaK0s);
    fCommonHistList->Add(fh1EtaLa);
    fCommonHistList->Add(fh1EtaALa);  
    fCommonHistList->Add(fh3InvMassEtaTrackPtK0s);
    fCommonHistList->Add(fh3InvMassEtaTrackPtLa);
    fCommonHistList->Add(fh3InvMassEtaTrackPtALa);
    fCommonHistList->Add(fh1noAssociatedK0s); 
    fCommonHistList->Add(fh1TrackMultCone);
    fCommonHistList->Add(fh2TrackMultCone);
    fCommonHistList->Add(fh2MCgenK0Cone);
    fCommonHistList->Add(fh2MCgenLaCone);
    fCommonHistList->Add(fh2MCgenALaCone);
    fCommonHistList->Add(fh2MCEtagenK0Cone);
    fCommonHistList->Add(fh2MCEtagenLaCone);
    fCommonHistList->Add(fh2MCEtagenALaCone);
    fCommonHistList->Add(fh1FFIMK0ConeSmear);
    fCommonHistList->Add(fh1FFIMLaConeSmear);
    fCommonHistList->Add(fh1FFIMALaConeSmear);
    fCommonHistList->Add(fh3MCrecK0Cone);
    fCommonHistList->Add(fh3MCrecLaCone);
    fCommonHistList->Add(fh3MCrecALaCone); 
    fCommonHistList->Add(fh3MCrecK0ConeSmear);
    fCommonHistList->Add(fh3MCrecLaConeSmear);
    fCommonHistList->Add(fh3MCrecALaConeSmear); 
    fCommonHistList->Add(fh3SecContinCone);
    fCommonHistList->Add(fh3StrContinCone);
    fCommonHistList->Add(fh3IMK0PerpCone);
    fCommonHistList->Add(fh3IMLaPerpCone);
    fCommonHistList->Add(fh3IMALaPerpCone);
    fCommonHistList->Add(fh3IMK0MedianCone);
    fCommonHistList->Add(fh3IMLaMedianCone);
    fCommonHistList->Add(fh3IMALaMedianCone);
    fCommonHistList->Add(fh1MCMultiplicityPrimary);       
    fCommonHistList->Add(fh1MCMultiplicityTracks);       
    fCommonHistList->Add(fh1MCmotherLa);
    fCommonHistList->Add(fh1MCmotherALa);
    fCommonHistList->Add(fh3FeedDownLa);
    fCommonHistList->Add(fh3FeedDownALa);
    fCommonHistList->Add(fh1MCProdRadiusK0s);
    fCommonHistList->Add(fh1MCProdRadiusLambda);
    fCommonHistList->Add(fh1MCProdRadiusAntiLambda);
    fCommonHistList->Add(fh1MCPtV0s);                    
    fCommonHistList->Add(fh1MCPtK0s);
    fCommonHistList->Add(fh1MCPtLambda);    
    fCommonHistList->Add(fh1MCPtAntiLambda);
    fCommonHistList->Add(fh1MCXiPt);
    fCommonHistList->Add(fh1MCXibarPt);
    fCommonHistList->Add(fh2MCEtaVsPtK0s); 
    fCommonHistList->Add(fh2MCEtaVsPtLa);
    fCommonHistList->Add(fh2MCEtaVsPtALa);     
    fCommonHistList->Add(fh1MCRapK0s);
    fCommonHistList->Add(fh1MCRapLambda);
    fCommonHistList->Add(fh1MCRapAntiLambda);   
    fCommonHistList->Add(fh1MCEtaAllK0s);
    fCommonHistList->Add(fh1MCEtaK0s);
    fCommonHistList->Add(fh1MCEtaLambda);
    fCommonHistList->Add(fh1MCEtaAntiLambda);         

    fV0QAK0->AddToOutput(fCommonHistList);
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecCutsK0Evt->AddToOutput(fCommonHistList);
    fFFHistosIMK0AllEvt->AddToOutput(fCommonHistList);
    fFFHistosIMK0Jet->AddToOutput(fCommonHistList);
    fFFHistosIMK0Cone->AddToOutput(fCommonHistList);
    fFFHistosPhiCorrIMK0->AddToOutput(fCommonHistList);
    fFFHistosIMLaAllEvt->AddToOutput(fCommonHistList);
    fFFHistosIMLaJet->AddToOutput(fCommonHistList);
    fFFHistosIMLaCone->AddToOutput(fCommonHistList);
    fFFHistosPhiCorrIMLa->AddToOutput(fCommonHistList);
    fFFHistosIMALaAllEvt->AddToOutput(fCommonHistList);
    fFFHistosIMALaJet->AddToOutput(fCommonHistList);
    fFFHistosIMALaCone->AddToOutput(fCommonHistList);
    fFFHistosPhiCorrIMALa->AddToOutput(fCommonHistList);
    
 
  }

   // =========== Switch on Sumw2 for all histos ===========
   for (Int_t i=0; i<fCommonHistList->GetEntries(); ++i){

   TH1 *h1 = dynamic_cast<TH1*>(fCommonHistList->At(i));
 
   if (h1) h1->Sumw2();//The error per bin will be computed as sqrt(sum of squares of weight) for each bin
    else{   
      THnSparse *hnSparse = dynamic_cast<THnSparse*>(fCommonHistList->At(i));
      if(hnSparse) hnSparse->Sumw2();
    }

  }
  TH1::AddDirectory(oldStatus);
 
}

//_______________________________________________
void AliAnalysisTaskJetChem::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if(fDebug > 1) Printf("AliAnalysisTaskJetChem::UserExec()");
   
  if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);

   // Trigger selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  
  //std::cout<<"inputHandler->IsEventSelected(): "<<inputHandler->IsEventSelected()<<std::endl;
  //std::cout<<"fEvtSelectionMask: "<<fEvtSelectionMask<<std::endl;
  
  if(!(inputHandler->IsEventSelected() & fEvtSelectionMask)){
    //std::cout<<"########event rejected!!############"<<std::endl;
    fh1EvtSelection->Fill(1.);
    if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
    PostData(1, fCommonHistList);
    return;
  } 
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());//casting of pointers for inherited class, only for ESDs
  if(!fESD){
    if(fDebug>3) Printf("%s:%d ESDEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  fMCEvent = MCEvent();
  if(!fMCEvent){
    if(fDebug>3) Printf("%s:%d MCEvent not found in the input", (char*)__FILE__,__LINE__);
  }
  
  // get AOD event from input/output         
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
    if( outHandler && outHandler->InheritsFrom("AliAODHandler") ){
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
  
  //primary vertex position:
  AliAODVertex *myPrimaryVertex = NULL;
  myPrimaryVertex = (AliAODVertex*)fAOD->GetPrimaryVertex();
  if (!myPrimaryVertex) return;
  fh1Evt->Fill(1.);//fill in every event that was accessed with InputHandler

  // event selection  *****************************************
  
  // *** vertex cut ***
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  Int_t nTracksPrim = primVtx->GetNContributors();
  fh1VertexNContributors->Fill(nTracksPrim);
  
  if (fDebug > 1) Printf("%s:%d primary vertex selection: %d", (char*)__FILE__,__LINE__,nTracksPrim);
  //if(!nTracksPrim){
  if(nTracksPrim <= 2){
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
  
  // accepts only events that have same "primary" and SPD vertex, special issue of LHC11h PbPb data

  //fAOD: pointer to global primary vertex
  
  const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
  
  if (TMath::Abs(spdVtx->GetZ() - primVtx->GetZ())>fDeltaVertexZ) { if (fDebug > 1) Printf("deltaZVertex: event REJECTED..."); return;}


  //check for vertex radius to be smaller than 1 cm, (that was first applied by Vit Kucera in his analysis)

  Double_t vtxX = primVtx->GetX();
  Double_t vtxY = primVtx->GetY();
 
  if(TMath::Sqrt(vtxX*vtxX + vtxY*vtxY)>=1){
    if (fDebug > 1) Printf("%s:%d primary vertex r = %f: event REJECTED...",(char*)__FILE__,__LINE__,TMath::Sqrt(vtxX*vtxX + vtxY*vtxY)); 
    return; 
  }
  

  TString primVtxName(primVtx->GetName());
  
  if(primVtxName.CompareTo("TPCVertex",TString::kIgnoreCase) == 1){
    if (fDebug > 1) Printf("%s:%d primary vertex selection: TPC vertex, event REJECTED...",(char*)__FILE__,__LINE__);
    fh1EvtSelection->Fill(5.);
    PostData(1, fCommonHistList);
    return;
  }
  
  Bool_t selectedHelper = AliAnalysisHelperJetTasks::Selected();
  if(!selectedHelper){
    fh1EvtSelection->Fill(6.);
    PostData(1, fCommonHistList);
    return;
  }

  Double_t centPercent = -1;
  if(fEventClass>0){
    Int_t cl = 0;
    if(handler && handler->InheritsFrom("AliAODInputHandler")){ 
    
      centPercent = fAOD->GetHeader()->GetCentrality();

      fh1EvtAllCent->Fill(centPercent);
     
      if(centPercent < 0) cl = -1;
      if(centPercent >= 0)  cl = 1;
      if(centPercent > 10) cl = 2; //standard PWG-JE binning
      if(centPercent > 30) cl = 3;
      if(centPercent > 50) cl = 4;
      if(centPercent > 80) cl = 5; //takes centralities higher than my upper edge of 80%, not to be used
  
    }
    else {

      cl = AliAnalysisHelperJetTasks::EventClass();

      if(fESD) centPercent = fESD->GetCentrality()->GetCentralityPercentile("V0M"); //ESD JetServices Task has the centrality binning 0-10,10-30,30-50,50-80
      fh1EvtAllCent->Fill(centPercent);
    }
    
    if(cl!=fEventClass){ // event not in selected event class, reject event#########################################
     
      if (fDebug > 1) Printf("%s:%d event not in selected event class: event REJECTED ...",(char*)__FILE__,__LINE__);
      fh1EvtSelection->Fill(2.);
      PostData(1, fCommonHistList);
      return;
    }
  }//end if fEventClass > 0
  
  
  if (fDebug > 1) Printf("%s:%d event ACCEPTED ...",(char*)__FILE__,__LINE__); 
  
  //test test
  //Printf("Analysis event #%5d", (Int_t) fEntry);

  fh1EvtSelection->Fill(0.);
  fh1EvtCent->Fill(centPercent);
    
  //___ get MC information __________________________________________________________________

 
  Double_t ptHard = 0.; //parton energy bins -> energy of particle
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  
  if(fMCEvent){
     AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
     AliGenPythiaEventHeader*  pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);//check usage of Pythia (pp) or Hijing (PbPb)
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
	 if(fDebug>3) Printf("%s:%d no pythiaGenHeader or hjingGenHeader found", (char*)__FILE__,__LINE__);
       } else {
	 if(fDebug>3) Printf("%s:%d hijingGenHeader found", (char*)__FILE__,__LINE__);
       }
     }
     
     fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
  }
  
    //____ fetch jets _______________________________________________________________

  Int_t nJCuts = GetListOfJets(fJetsRecCuts, kJetsRecAcceptance);//fetch list with jets

  Int_t nRecJetsCuts = 0;                                        //number of reconstructed jets after jet cuts
  if(nJCuts>=0) nRecJetsCuts = fJetsRecCuts->GetEntries(); 
  if(fDebug>2)Printf("%s:%d Selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  if(nRecJetsCuts != nJCuts) Printf("%s:%d Mismatch selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  fh1nRecJetsCuts->Fill(nRecJetsCuts);


  //____ fetch background clusters ___________________________________________________
  if(fBranchRecBckgClusters.Length() != 0){

    Int_t nBJ = GetListOfBckgJets(fBckgJetsRec, kJetsRec);
    Int_t nRecBckgJets = 0;
    if(nBJ>=0) nRecBckgJets = fBckgJetsRec->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);
    if(nBJ != nRecBckgJets) Printf("%s:%d Mismatch Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);
  }

  
  //____ fetch reconstructed particles __________________________________________________________
 
  Int_t nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);//all tracks of event
  if(fDebug>2)Printf("%s:%d selected reconstructed tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());
  if(fTracksRecCuts->GetEntries() != nTCuts) 
    Printf("%s:%d Mismatch selected reconstructed tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());
  fh1EvtMult->Fill(fTracksRecCuts->GetEntries());

  Int_t nK0s = GetListOfV0s(fListK0s,fK0Type,kK0,myPrimaryVertex,fAOD);//all V0s in event with K0s assumption
  
  //std::cout<< "nK0s: "<<nK0s<<std::endl;

  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  if(nK0s != fListK0s->GetEntries()) Printf("%s:%d Mismatch selected K0s: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  fh1K0Mult->Fill(fListK0s->GetEntries());

  
  Int_t nLa = GetListOfV0s(fListLa,fLaType,kLambda,myPrimaryVertex,fAOD);//all V0s in event with Lambda particle assumption 
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nLa,fListLa->GetEntries());
  if(nLa != fListLa->GetEntries()) Printf("%s:%d Mismatch selected La: %d %d",(char*)__FILE__,__LINE__,nLa,fListLa->GetEntries());
  fh1LaMult->Fill(fListLa->GetEntries());
 
  Int_t nALa = GetListOfV0s(fListALa,fALaType,kAntiLambda,myPrimaryVertex,fAOD);//all V0s in event with Antilambda particle assumption
  if(fDebug>2)Printf("%s:%d Selected Rec tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nALa,fListALa->GetEntries());
  if(nALa != fListALa->GetEntries()) Printf("%s:%d Mismatch selected ALa: %d %d",(char*)__FILE__,__LINE__,nALa,fListALa->GetEntries());
  fh1ALaMult->Fill(fListALa->GetEntries());


    
  //fetch MC gen particles_______________________________________________________

  if(fAnalysisMC){ // here 

    //fill feeddown histo for associated particles

    // Access MC generated particles, fill TLists and histograms :
    
    Int_t nMCgenK0s = GetListOfMCParticles(fListMCgenK0s,kK0,fAOD); //fill TList with MC generated primary true K0s (list to fill, particletype, mc aod event)
    if(nMCgenK0s != fListMCgenK0s->GetEntries()) Printf("%s:%d Mismatch selected MCgenK0s: %d %d",(char*)__FILE__,__LINE__,nMCgenK0s,fListMCgenK0s->GetEntries());
    
    
    for(Int_t it=0; it<fListMCgenK0s->GetSize(); ++it){ // loop MC generated K0s, filling histograms
      
      AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0s->At(it));
      if(!mcp0) continue;
      
      //MC gen K0s                  
      
      Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      
      fh1MCEtaK0s->Fill(fEtaCurrentPart); 
      fh1MCRapK0s->Fill(fRapCurrentPart);
      fh1MCPtK0s->Fill(fPtCurrentPart);	  
      
      fh2MCEtaVsPtK0s->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered
      
    }//end of the loop
    
    
    Int_t nMCgenLa = GetListOfMCParticles(fListMCgenLa,kLambda,fAOD); //fill TList with MC generated primary true Lambdas (list to fill, particletype, mc aod event)
    if(nMCgenLa != fListMCgenLa->GetEntries()) Printf("%s:%d Mismatch selected MCgenLa: %d %d",(char*)__FILE__,__LINE__,nMCgenLa,fListMCgenLa->GetEntries());

	
    for(Int_t it=0; it<fListMCgenLa->GetSize(); ++it){ // loop MC generated La, filling histograms
      
      AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLa->At(it));
      if(!mcp0) continue;
	  
      //MC gen Lambdas  
      
      Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      
      fh1MCEtaLambda->Fill(fEtaCurrentPart); 
      fh1MCRapLambda->Fill(fRapCurrentPart);
      fh1MCPtLambda->Fill(fPtCurrentPart);	  
      fh2MCEtaVsPtLa->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered
      
    }//end of the loop


    Int_t nMCgenALa = GetListOfMCParticles(fListMCgenALa,kAntiLambda,fAOD); //fill TList with MC generated primary true Antilambdas (list to fill, particletype, mc aod event)
    if(nMCgenALa != fListMCgenALa->GetEntries()) Printf("%s:%d Mismatch selected MCgenALa: %d %d",(char*)__FILE__,__LINE__,nMCgenALa,fListMCgenALa->GetEntries());
  
   	
    for(Int_t it=0; it<fListMCgenALa->GetSize(); ++it){ // loop MC generated ALa, filling histograms
      
      AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenALa->At(it));
      if(!mcp0) continue;
      
      //MC gen Antilambdas                  
      
      Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      
      fh1MCEtaAntiLambda->Fill(fEtaCurrentPart); 
      fh1MCRapAntiLambda->Fill(fRapCurrentPart);
      fh1MCPtAntiLambda->Fill(fPtCurrentPart);	  
      fh2MCEtaVsPtALa->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered
	
    }//end of the loop

	

  //loop over MC feeddown candidates in TList

    //.... 

	
  } //end MCAnalysis part for gen particles
      
      
  // ___ V0 QA + K0s + La + ALa pt spectra all events _______________________________________________
  
  Double_t lPrimaryVtxPosition[3];
  Double_t lV0Position[3];
  lPrimaryVtxPosition[0] = primVtx->GetX();
  lPrimaryVtxPosition[1] = primVtx->GetY();
  lPrimaryVtxPosition[2] = primVtx->GetZ();
  
  //------------------------------------------
 
  for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 
         
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
    if(!v0) continue;
    
    // VO's main characteristics to check the reconstruction cuts
    
    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
    Double_t invMK0s=0;
    Double_t trackPt=0;   
    Double_t fV0Radius      = -999;
    Double_t fDcaV0Daughters = v0->DcaV0Daughters();
    Double_t fDcaPosToPrimVertex = v0->DcaPosToPrimVertex();//IP of positive charged daughter
    Double_t fDcaNegToPrimVertex = v0->DcaNegToPrimVertex();//IP of negative charged daughter
    Int_t negDaughterpdg = 0;
    Int_t posDaughterpdg = 0;
    Int_t motherType = 0;
    Int_t v0Label = -1;
    Double_t MCPt = 0;
    Bool_t fPhysicalPrimary = kFALSE;//don't use IsPhysicalPrimary() anymore for MC analysis, use instead 2D distance from primary to secondary vertex
    Int_t MCv0PdgCode = 0;

    AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));  
    AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));   
    
    Double_t PosEta = trackPos->AliAODTrack::Eta();//daughter track charge is sometimes wrong here, account for that!!!
    Double_t NegEta = trackNeg->AliAODTrack::Eta();
    
    //Double_t trackPosNcls = trackPos->GetTPCNcls();//Get number of clusters for positive charged tracks
    //Double_t trackNegNcls = trackNeg->GetTPCNcls();//Get number of clusters for negative charged tracks
    
    CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
    
    Double_t jetPt = fFFIMJetPtMin; // assign pro forma jet energy
    //Double_t fRap = v0->RapK0Short();
    Double_t fEta = v0->PseudoRapV0();
    
    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);

    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();
    
    Double_t fV0mom[3];
    
    fV0mom[0]=v0->MomV0X();
    fV0mom[1]=v0->MomV0Y();
    fV0mom[2]=v0->MomV0Z();
    Double_t fV0TotalMomentum = TMath::Sqrt(fV0mom[0]*fV0mom[0]+fV0mom[1]*fV0mom[1]+fV0mom[2]*fV0mom[2]);
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
    
    fV0QAK0->FillTrackQA(v0->Eta(), TVector2::Phi_0_2pi(v0->Phi()), v0->Pt()); 
    fFFHistosIMK0AllEvt->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaK0s->Fill(fEta);
    if(fAnalysisMC){
      TList *listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kK0, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);

      //std::cout<<"mclabelcheck: "<<mclabelcheck<<std::endl;
      //std::cout<<"IsPhysicalPrimary: "<<fPhysicalPrimary<<std::endl;

      if(mclabelcheck == kFALSE)continue;
      fh3InvMassEtaTrackPtK0s->Fill(fEta,invMK0s,trackPt);
      fh1PtMCK0s->Fill(MCPt);
    }
 

    fh1V0Eta->Fill(fEta);
    fh1V0totMom->Fill(fV0TotalMomentum);
    fh1CosPointAngle->Fill(fV0cosPointAngle);
    fh1DecayLengthV0->Fill(fV0DecayLength);
    fh1V0Radius->Fill(fV0Radius);
    fh1DcaV0Daughters->Fill(fDcaV0Daughters);
    fh1DcaPosToPrimVertex->Fill(fDcaPosToPrimVertex);
    fh1DcaNegToPrimVertex->Fill(fDcaNegToPrimVertex);
    fh1trackPosEta->Fill(PosEta);
    fh1trackNegEta->Fill(NegEta);  
  }
  

  // __La pt spectra all events _______________________________________________

    
  for(Int_t it=0; it<fListLa->GetSize(); ++it){ // loop all La 
      
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLa->At(it));
    if(!v0) continue;
    
    // VO's main characteristics to check the reconstruction cuts
    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
    Double_t invMLa =0;
    Double_t trackPt=0;
    Double_t fV0Radius      = -999;
    Double_t fDcaV0Daughters = v0->DcaV0Daughters();
    Double_t fDcaPosToPrimVertex = v0->DcaPosToPrimVertex();//IP of positive charged daughter
    Double_t fDcaNegToPrimVertex = v0->DcaNegToPrimVertex();//IP of negative charged daughter
    Int_t negDaughterpdg = 0;
    Int_t posDaughterpdg = 0;
    Int_t motherType = 0;
    Int_t v0Label = -1;
    Double_t MCPt = 0;
    Bool_t fPhysicalPrimary = kFALSE;
    Int_t MCv0PdgCode = 0;
    AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));  
    AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));   
    
    //Double_t trackPosNcls = trackPos->GetTPCNcls();//Get number of clusters for positive charged tracks
    //Double_t trackNegNcls = trackNeg->GetTPCNcls();//Get number of clusters for negative charged tracks
    
    Double_t PosEta = trackPos->AliAODTrack::Eta();//daughter track charge is sometimes wrong here, account for that!!!
    Double_t NegEta = trackNeg->AliAODTrack::Eta();
    
    CalculateInvMass(v0, kLambda, invMLa, trackPt);//function to calculate invMass with TLorentzVector class
    
    
    Double_t jetPt = fFFIMJetPtMin; // assign pro forma jet energy
    // Double_t fRap = v0->Y(3122);
    Double_t fEta = v0->PseudoRapV0();
    
    Double_t fV0mom[3];
    
    fV0mom[0]=v0->MomV0X();
    fV0mom[1]=v0->MomV0Y();
    fV0mom[2]=v0->MomV0Z();
    Double_t fV0TotalMomentum = TMath::Sqrt(fV0mom[0]*fV0mom[0]+fV0mom[1]*fV0mom[1]+fV0mom[2]*fV0mom[2]);
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();  
    
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
    
    fFFHistosIMLaAllEvt->FillFF(trackPt, invMLa, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaLa->Fill(fEta);
    if(fAnalysisMC){     
      TList* listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
      if(mclabelcheck == kFALSE)continue;
       fh3InvMassEtaTrackPtLa->Fill(fEta,invMLa,trackPt);
      fh1PtMCLa->Fill(MCPt);
    }
    fh1V0Eta->Fill(fEta);
    fh1V0totMom->Fill(fV0TotalMomentum);
    fh1CosPointAngle->Fill(fV0cosPointAngle);
    fh1DecayLengthV0->Fill(fV0DecayLength);
    fh1V0Radius->Fill(fV0Radius);
    fh1DcaV0Daughters->Fill(fDcaV0Daughters);
    fh1DcaPosToPrimVertex->Fill(fDcaPosToPrimVertex);
    fh1DcaNegToPrimVertex->Fill(fDcaNegToPrimVertex);
    fh1trackPosEta->Fill(PosEta);
    fh1trackNegEta->Fill(NegEta);
  }
  
  // __ALa pt spectra all events _______________________________________________
    
  for(Int_t it=0; it<fListALa->GetSize(); ++it){ // loop all ALa 
    
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALa->At(it));
    if(!v0) continue;
      

    //VO's main characteristics to check the reconstruction cuts
    Double_t invMALa =0;
    Double_t trackPt=0;
    Double_t fV0Radius      = -999;
    Double_t fDcaV0Daughters = v0->DcaV0Daughters();
    Double_t fDcaPosToPrimVertex = v0->DcaPosToPrimVertex();//IP of positive charged daughter
    Double_t fDcaNegToPrimVertex = v0->DcaNegToPrimVertex();//IP of negative charged daughter
    Int_t negDaughterpdg = 0;
    Int_t posDaughterpdg = 0;
    Int_t motherType = 0;
    Int_t v0Label = -1;
    Double_t MCPt = 0;
    Bool_t fPhysicalPrimary = kFALSE;
    Int_t MCv0PdgCode = 0;
    
    AliAODTrack *trackPos = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));  
    AliAODTrack *trackNeg = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));   
      
    Double_t PosEta = trackPos->AliAODTrack::Eta();//daughter track charge is sometimes wrong here, account for that!!!
    Double_t NegEta = trackNeg->AliAODTrack::Eta();
       
    //Double_t trackPosNcls = trackPos->GetTPCNcls();//Get number of clusters for positive charged tracks //not used anymore by Strangeness PAG group
    //Double_t trackNegNcls = trackNeg->GetTPCNcls();//Get number of clusters for negative charged tracks
    
    CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
      
    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
    Double_t jetPt = fFFIMJetPtMin; // assign pro forma jet energy
    //      Double_t fRap = v0->Y(-3122);
    Double_t fEta = v0->PseudoRapV0();
    
    Double_t fV0mom[3];
    
    fV0mom[0]=v0->MomV0X();
    fV0mom[1]=v0->MomV0Y();
    fV0mom[2]=v0->MomV0Z();
    Double_t fV0TotalMomentum = TMath::Sqrt(fV0mom[0]*fV0mom[0]+fV0mom[1]*fV0mom[1]+fV0mom[2]*fV0mom[2]);

    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();  
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
       
    fFFHistosIMALaAllEvt->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaALa->Fill(fEta);
    if(fAnalysisMC){
      TList* listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
      if(mclabelcheck == kFALSE)continue;
      fh3InvMassEtaTrackPtALa->Fill(fEta,invMALa,trackPt);
      fh1PtMCALa->Fill(MCPt);
    }
    fh1V0Eta->Fill(fEta);
    fh1V0totMom->Fill(fV0TotalMomentum);
    fh1CosPointAngle->Fill(fV0cosPointAngle);
    fh1DecayLengthV0->Fill(fV0DecayLength);
    fh1V0Radius->Fill(fV0Radius);
    fh1DcaV0Daughters->Fill(fDcaV0Daughters);
    fh1DcaPosToPrimVertex->Fill(fDcaPosToPrimVertex);
    fh1DcaNegToPrimVertex->Fill(fDcaNegToPrimVertex);
    fh1trackPosEta->Fill(PosEta);
    fh1trackNegEta->Fill(NegEta);
  }
  

  //____ fill all jet related histos  ________________________________________________________________________________________________________________________
 
  
  //fill jet histos in general
  for(Int_t ij=0; ij<nRecJetsCuts; ++ij){                               // ij is an index running over the list of the reconstructed jets after cuts, all jets in event
    
    AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));

    Double_t jetPt  = jet->Pt();
    Double_t jetEta = jet->Eta();
    Double_t jetPhi = jet->Phi();

    //if(ij==0){ // loop over leading jets for ij = 0, for ij>= 0 look into all jets

    if(ij>=0){//all jets in event

      TList* jettracklist = new TList();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
      Int_t njetTracks    = 0;
 
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);// list of jet tracks from trackrefs
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);  // fill list of tracks in cone around jet axis with cone Radius (= 0.4 standard)
      }
      
      if(GetFFMinNTracks()>0 && jettracklist->GetSize() <= GetFFMinNTracks()) isBadJet = kTRUE; // reject jets with less tracks than fFFMinNTracks
      if(isBadJet) continue; // rejects jets in which no track has a track pt higher than 5 GeV/c (see AddTask macro)
     
      Float_t fJetAreaMin = 0.6*TMath::Pi()*GetFFRadius()*GetFFRadius(); // minimum jet area cut
      if (jet->EffectiveAreaCharged() < fJetAreaMin)continue;
      //std::cout<<"GetFFRadius(): "<<GetFFRadius()<<std::endl;
      //std::cout<<"fJetAreaMin: "<<fJetAreaMin<<std::endl;
      //std::cout<<"jet->EffectiveAreaCharged()"<<jet->EffectiveAreaCharged()<<std::endl;


      fh1JetEta->Fill(jetEta);        
      fh1JetPhi->Fill(jetPhi);                
      fh2JetEtaPhi->Fill(jetEta,jetPhi);  
  
      // printf("pT = %f, eta = %f, phi = %f, leadtr pt = %f\n, ",jetPt,jetEta,jetphi,leadtrack);

      for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all particles in jet
	
	AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//all tracks in jet cone	
	if(!trackVP)continue;
	
	Float_t trackPt = trackVP->Pt();//transversal momentum of jet particle
	Float_t trackEta = trackVP->Eta();

	Float_t leadtrPt = 0;

	if(trackPt > 5.){leadtrPt = trackPt;}

	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	fFFHistosRecCuts->FillFF(trackPt, jetPt, incrementJetPt);
	if(nK0s>0) fFFHistosRecCutsK0Evt->FillFF(trackPt, jetPt, incrementJetPt);
	fh2FFJetTrackEta->Fill(trackEta,jetPt);


      }
     
      njetTracks = jettracklist->GetSize();

      //____________________________________________________________________________________________________________________      
      //alternative method to estimate secondary constribution in jet cone (second method you can see below in rec. K0s loop & rec. Lambdas loop & rec. Antilambdas loop)

      if(fAnalysisMC){

	TList *list = fAOD->GetList();  
	AliAODMCHeader *mcHeadr=(AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());	  
	if(!mcHeadr)continue;

	Double_t mcXv=0., mcYv=0., mcZv=0.;//MC primary vertex position

	mcXv=mcHeadr->GetVtxX(); mcYv=mcHeadr->GetVtxY(); mcZv=mcHeadr->GetVtxZ(); // position of the MC primary vertex

	for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all tracks in the jet
 	  
	  AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//track in jet cone	
	  if(!trackVP)continue;
	  AliAODTrack *tr = dynamic_cast<AliAODTrack*> (trackVP);                   //fetch one jet track from the TList
	  if(!tr)continue;
	  
	  //get MC label information
	  TList *mclist = fAOD->GetList();                                           
	 
	  //fetch the MC stack
	  TClonesArray* stackMC = (TClonesArray*)mclist->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
	  if (!stackMC) {Printf("ERROR: stack not available");}

	  else {
	    
	    Int_t trackLabel = TMath::Abs(tr->GetLabel());                       //fetch jet track label in MC stack
	    
	    AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(stackMC->At(trackLabel));  //fetch MC gen. particle for rec. jet track

	    if(!part)continue;  //skip non-existing objects     
	    

	    //Bool_t IsPhysicalPrimary = part->IsPhysicalPrimary();//not recommended to check, better use distance between primary vertex and secondary vertex
	    
	    Float_t fDistPrimaryMax = 0.01;
	    // Get the distance between production point of the MC mother particle and the primary vertex
	    
	    Double_t dx = mcXv-part->Xv();//mc primary vertex - mc gen. v0 vertex 
	    Double_t dy = mcYv-part->Yv();
	    Double_t dz = mcZv-part->Zv();
	    
	    Float_t fDistPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	    Bool_t fPhysicalPrimary = (fDistPrimary < fDistPrimaryMax);
 
	    // std::cout<<"fDistPrimary"<<fDistPrimary<<std::endl;
	    // std::cout<<"fPhysicalPrimary"<<fPhysicalPrimary<<std::endl;

	    if(!fPhysicalPrimary)continue;//rejects Kstar and other strong decaying particles from Secondary Contamination
	    
	    Bool_t isFromStrange = kFALSE;// flag to check whether particle has strange mother
	    
	    Int_t iMother = part->GetMother(); //get mother MC gen. particle label
	    
	    if(iMother >= 0){
	      AliAODMCParticle *partM = dynamic_cast<AliAODMCParticle*>(stackMC->At(iMother)); //fetch mother of MC gen. particle
	      if(!partM) continue;

	      Int_t codeM =  TMath::Abs(partM->GetPdgCode());                                 //mothers pdg code
	      
	      Int_t mfl = Int_t (codeM/ TMath::Power(10, Int_t(TMath::Log10(codeM))));        //asks for first number of mothers pdg code (strange particles always start with 3..)
	      
	      if  (mfl == 3 && codeM != 3) isFromStrange = kTRUE;
	    }
    
	    //cut on primary particles:




	    if(isFromStrange == kTRUE){

	      Double_t trackPt = part->Pt();
	      Double_t trackEta = part->Eta();
	      fh3StrContinCone->Fill(jetPt, trackPt, trackEta);//MC gen. particle parameters, but rec. jet pt
	     	      
	      }//isFromStrange is kTRUE  
	  } //end else
	}//end loop over jet tracks
	
      }// end fAnalysisMC
      

      fh1TrackMultCone->Fill(njetTracks);
      fh2TrackMultCone->Fill(njetTracks,jetPt);      
     	  
      // ---- K0s ---- 
      
      // fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
      
      for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
	if(!v0) continue;//rejection of events with no V0 vertex

	Double_t v0Mom[3];
	v0->PxPyPz(v0Mom);
	TVector3 v0MomVect(v0Mom);
	
	Double_t dPhiJetK0 = (jet->MomentumVector()->Vect()).DeltaPhi(v0MomVect);
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	Double_t invMK0s =0;
	Double_t trackPt=0;	
	CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fFFHistosIMK0Jet->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
 	fFFHistosPhiCorrIMK0->FillPhiCorr(trackPt,TVector2::Phi_0_2pi(dPhiJetK0),invMK0s);
	
	if(dPhiJetK0<fh1dPhiJetK0->GetXaxis()->GetXmin()) dPhiJetK0 += 2*TMath::Pi();
	fh1dPhiJetK0->Fill(dPhiJetK0);
	
      }

      if(fListK0s->GetSize() == 0){ // no K0: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMK0Jet->FillFF(-1, -1, jetPt, incrementJetPt);
      }
      
      //____fetch reconstructed K0s in cone around jet axis:_______________________________________________________________________________
      
      TList* jetConeK0list = new TList();

      Double_t sumPtK0     = 0.;
      
      Bool_t isBadJetK0    = kFALSE; // dummy, do not use


      GetTracksInCone(fListK0s, jetConeK0list, jet, GetFFRadius(), sumPtK0, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0); //reconstructed K0s in cone around jet axis
    
      if(fDebug>2)Printf("%s:%d nK0s total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetConeK0list->GetEntries(),GetFFRadius());
      
      
      for(Int_t it=0; it<jetConeK0list->GetSize(); ++it){ // loop for K0s in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeK0list->At(it));
	if(!v0) continue;
	
	Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
	Double_t invMK0s =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class

	if(fAnalysisMC){
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE){fh1FFIMK0ConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	}

	fFFHistosIMK0Cone->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
      }
      
      
      if(jetConeK0list->GetSize() == 0){ // no K0: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMK0Cone->FillFF(-1, -1, jetPt, incrementJetPt);
	if(fAnalysisMC){
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE){fh1FFIMK0ConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	}
      }    

      //fetch particles in perpendicular cone to estimate UE event contribution to particle spectrum
      //these perpendicular cone particle spectra serve to subtract the particles in jet cones, that are stemming from the Underlying event, on a statistical basis
      //for normalization the common jet pT spectrum is used: fh1FFIMK0Cone, fh1FFIMLaCone and fh1FFIMALaCone
      
      //____fetch reconstructed K0s in cone perpendicular to jet axis:_______________________________________________________________________________
      
      TList* jetPerpConeK0list = new TList();
      
      Double_t sumPerpPtK0     = 0.;
      
      GetTracksInPerpCone(fListK0s, jetPerpConeK0list, jet, GetFFRadius(), sumPerpPtK0); //reconstructed K0s in cone around jet axis
      
      if(fDebug>2)Printf("%s:%d nK0s total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetPerpConeK0list->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetPerpConeK0list->GetSize(); ++it){ // loop for K0s in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeK0list->At(it));
	if(!v0) continue;
	
	Double_t invMPerpK0s =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kK0, invMPerpK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMK0PerpCone->Fill(jetPt, invMPerpK0s, trackPt); //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
      }
      
      
      if(jetPerpConeK0list->GetSize() == 0){ // no K0s in jet cone 
	
	//Bool_t incrementPerpJetPt = kTRUE;
	fh3IMK0PerpCone->Fill(jetPt, -1, -1);
      }
   
      // ____ rec K0s in median cluster___________________________________________________________________________________________________________ 
      
      TList* jetMedianConeK0list = new TList();
      
      AliAODJet* medianCluster = GetMedianCluster();

      Double_t sumMedianPtK0     = 0.;

      Bool_t isBadJetK0Median    = kFALSE; // dummy, do not use
     
      GetTracksInCone(fListK0s, jetMedianConeK0list, medianCluster, GetFFRadius(), sumMedianPtK0, 0., 0., isBadJetK0Median); //reconstructed K0s in median cone around jet axis
      //GetTracksInCone(fListK0s, jetConeK0list, jet, GetFFRadius(), sumPtK0, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0); //original use of function
      
      //cut parameters from Fragmentation Function task:
      //Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value, use GetFFMinLTrackPt()
      //Float_t fFFMaxTrackPt;    // reject jetscontaining any track with pt larger than this value, use GetFFMaxTrackPt()
      
      for(Int_t it=0; it<jetMedianConeK0list->GetSize(); ++it){ // loop for K0s in median cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetMedianConeK0list->At(it));
	if(!v0) continue;
	
	Double_t invMMedianK0s =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kK0, invMMedianK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMK0MedianCone->Fill(jetPt, invMMedianK0s, trackPt); //(x,y,z)
      }
      
      if(jetMedianConeK0list->GetSize() == 0){ // no K0s in median cluster cone 
	
	fh3IMK0MedianCone->Fill(jetPt, -1, -1);
      }
      
      //_________________________________________________________________________________________________________________________________________
      
      //____fetch reconstructed Lambdas in cone perpendicular to jet axis:__________________________________________________________________________
      
      TList* jetPerpConeLalist = new TList();
      
      Double_t sumPerpPtLa     = 0.;
      
      GetTracksInPerpCone(fListLa, jetPerpConeLalist, jet, GetFFRadius(), sumPerpPtLa); //reconstructed Lambdas in cone around jet axis //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
      
      if(fDebug>2)Printf("%s:%d nLa total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetPerpConeLalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetPerpConeLalist->GetSize(); ++it){ // loop for Lambdas in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeLalist->At(it));
	if(!v0) continue;
	
	Double_t invMPerpLa =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kLambda, invMPerpLa, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMLaPerpCone->Fill(jetPt, invMPerpLa, trackPt);
      }
      
      
      if(jetPerpConeLalist->GetSize() == 0){ // no Lambdas in jet
	
	fh3IMLaPerpCone->Fill(jetPt, -1, -1);
	
      }
      
      //__________________________________________________________________________________________________________________________________________
          // ____ rec Lambdas in median cluster___________________________________________________________________________________________________________ 
      
      TList* jetMedianConeLalist = new TList();
      
      //AliAODJet* medianCluster = GetMedianCluster(); //already loaded at part for K0s ??

      Double_t sumMedianPtLa     = 0.;
      Bool_t isBadJetLaMedian    = kFALSE; // dummy, do not use
     
      GetTracksInCone(fListLa, jetMedianConeLalist, medianCluster, GetFFRadius(), sumMedianPtLa, 0, 0, isBadJetLaMedian); //reconstructed Lambdas in median cone around jet axis
     
         //cut parameters from Fragmentation Function task:
      //Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value, use GetFFMinLTrackPt()
      //Float_t fFFMaxTrackPt;    // reject jets containing any track with pt larger than this value, use GetFFMaxTrackPt()
      
      for(Int_t it=0; it<jetMedianConeLalist->GetSize(); ++it){ // loop for Lambdas in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetMedianConeLalist->At(it));
	if(!v0) continue;
	
	Double_t invMMedianLa =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kLambda, invMMedianLa, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMLaMedianCone->Fill(jetPt, invMMedianLa, trackPt); //(x,y,z)
      }
      
      if(jetMedianConeLalist->GetSize() == 0){ // no Lambdas in median cluster cone 
	
	fh3IMLaMedianCone->Fill(jetPt, -1, -1);
      }
      
      //_________________________________________________________________________________________________________________________________________
      
      
      //____fetch reconstructed Antilambdas in cone perpendicular to jet axis:___________________________________________________________________
      
      TList* jetPerpConeALalist = new TList();
      
      Double_t sumPerpPtALa     = 0.;
      
      GetTracksInPerpCone(fListALa, jetPerpConeALalist, jet, GetFFRadius(), sumPerpPtALa); //reconstructed Antilambdas in cone around jet axis //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
      
      if(fDebug>2)Printf("%s:%d nALa total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetPerpConeALalist->GetEntries(),GetFFRadius());
            
      for(Int_t it=0; it<jetPerpConeALalist->GetSize(); ++it){ // loop for ALa in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeALalist->At(it));
	if(!v0) continue;
	
	Double_t invMPerpALa =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kAntiLambda, invMPerpALa, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMALaPerpCone->Fill(jetPt, invMPerpALa, trackPt);
	
      }
      
      
      if(jetPerpConeALalist->GetSize() == 0){ // no Antilambda 

	fh3IMALaPerpCone->Fill(jetPt, -1, -1);
	
      }
      

          // ____ rec Antilambdas in median cluster___________________________________________________________________________________________________________ 
      
      TList* jetMedianConeALalist = new TList();
      
      //AliAODJet* medianCluster = GetMedianCluster(); //already loaded at part for K0s

      Double_t sumMedianPtALa     = 0.;
      
      Bool_t isBadJetALaMedian    = kFALSE; // dummy, do not use
     
      GetTracksInCone(fListALa, jetMedianConeALalist, medianCluster, GetFFRadius(), sumMedianPtALa, 0, 0, isBadJetALaMedian); //reconstructed Antilambdas in median cone around jet axis
     
         
      //cut parameters from Fragmentation Function task:
      //Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value, use GetFFMinLTrackPt()
      //Float_t fFFMaxTrackPt;    // reject jets containing any track with pt larger than this value, use GetFFMaxTrackPt()
      
      for(Int_t it=0; it<jetMedianConeALalist->GetSize(); ++it){ // loop for Antilambdas in median cluster cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetMedianConeALalist->At(it));
	if(!v0) continue;
	
	Double_t invMMedianALa =0;
	Double_t trackPt=0;
	
	CalculateInvMass(v0, kAntiLambda, invMMedianALa, trackPt);  //function to calculate invMass with TLorentzVector class
	
	fh3IMALaMedianCone->Fill(jetPt, invMMedianALa, trackPt); //(x,y,z)
      }
      
      if(jetMedianConeALalist->GetSize() == 0){ // no Antilambdas in median cluster cone 
	
	fh3IMALaMedianCone->Fill(jetPt, -1, -1);
      }
      
      //_________________________________________________________________________________________________________________________________________
      


      //MC Analysis, fetch MC truth in jet cones, denominator of rec. efficiency in jet cones
      //__________________________________________________________________________________________________________________________________________
      
      if(fAnalysisMC){    
      
	//fill feeddown candidates from TList	
	//std::cout<<"fListFeeddownLaCand entries: "<<fListFeeddownLaCand->GetSize()<<std::endl;


	Double_t sumPtFDLa     = 0.;
	Bool_t isBadJetFDLa    = kFALSE; // dummy, do not use
	
	GetTracksInCone(fListFeeddownLaCand, jetConeFDLalist, jet, GetFFRadius(), sumPtFDLa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetFDLa);
	
	Double_t sumPtFDALa     = 0.;
	Bool_t isBadJetFDALa    = kFALSE; // dummy, do not use
	
	GetTracksInCone(fListFeeddownALaCand, jetConeFDALalist, jet, GetFFRadius(), sumPtFDALa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetFDALa);
	
	for(Int_t it=0; it<jetConeFDLalist->GetSize(); ++it){ 
	  
	  AliAODv0* mcfd = dynamic_cast<AliAODv0*>(jetConeFDLalist->At(it));
	  if(!mcfd) continue;

	  Double_t invMLaFDcand = 0;
	  Double_t trackPt = 0;//pt of ass. particle, not used for the histos
	  
	  CalculateInvMass(mcfd, kLambda, invMLaFDcand, trackPt);
	  
	  //Get MC gen. Lambda transverse momentum
	  TClonesArray *st = 0x0;
	 	 	  
	  TList *lt = fAOD->GetList();
	  if(!lt)continue;
	 
	  st = (TClonesArray*)lt->FindObject(AliAODMCParticle::StdBranchName());
	  
	  Int_t v0lab= TMath::Abs(mcfd->GetLabel());

	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	    
	  Double_t genLaPt = mcp->Pt();
	  	  
	  fh3FeedDownLa->Fill(jetPt, invMLaFDcand, genLaPt);
	  
	}//end loop over feeddown candidates for Lambda particles in jet cone
	


	//feeddown for Antilambdas from Xi(bar)+ and Xi(bar)0 :

	for(Int_t it=0; it<jetConeFDALalist->GetSize(); ++it){ 
	  
	  AliAODv0* mcfd = dynamic_cast<AliAODv0*>(jetConeFDALalist->At(it));
	  if(!mcfd) continue;

	  Double_t invMALaFDcand = 0;
	  Double_t trackPt = 0;//pt of ass. particle, not used for the histos
	  
	  CalculateInvMass(mcfd, kAntiLambda, invMALaFDcand, trackPt);
	  
	  //Get MC gen. Antilambda transverse momentum
	  TClonesArray *st = 0x0;
	 	 	  
	  TList *lt = fAOD->GetList();
	  if(!lt)continue;
	 
	  st = (TClonesArray*)lt->FindObject(AliAODMCParticle::StdBranchName());
	  
	  Int_t v0lab= TMath::Abs(mcfd->GetLabel());

	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	    
	  Double_t genALaPt = mcp->Pt();
		  
	  fh3FeedDownALa->Fill(jetPt, invMALaFDcand, genALaPt);
	  
	}//end loop over feeddown candidates for Antilambda particles in jet cone
	
      	

	//____fetch MC generated K0s in cone around jet axis__(note: particles can stem from fragmentation but also from underlying event)________
	
	Double_t sumPtMCgenK0s   = 0.;
	Bool_t isBadJetMCgenK0s  = kFALSE; // dummy, do not use
	
	
	fListMCgenK0sCone = new TList();      //MC generated K0s in (only geometrical) jet cone (these are MC gen K0s falling geometrically into jet cone (R = 0.4) around jet axis, that was found by anti-kt jet finder, particles can stem from fragmentation but also from underlying event!!)
	
	//first: sampling MC gen K0s       
	
	GetTracksInCone(fListMCgenK0s, fListMCgenK0sCone, jet, GetFFRadius(), sumPtMCgenK0s, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMCgenK0s); //MC generated K0s in cone around jet axis 
	
	if(fDebug>2)Printf("%s:%d nMCgenK0s in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,fListMCgenK0sCone->GetEntries(),GetFFRadius());
	
	
	for(Int_t it=0; it<fListMCgenK0sCone->GetSize(); ++it){ // loop MC generated K0s in cone around jet axis
	  
	  AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0sCone->At(it));
	  if(!mcp0) continue;
	  
	  //Double_t fRapMCgenK0s   = MyRapidity(mcp0->E(),mcp0->Pz());//get rec. particle in cone information
	  Double_t fEtaMCgenK0s   = mcp0->Eta();
	  Double_t fPtMCgenK0s    = mcp0->Pt();
	  
	  fh2MCgenK0Cone->Fill(jetPt,fPtMCgenK0s); 
	  fh2MCEtagenK0Cone->Fill(jetPt,fEtaMCgenK0s);
	  
	}
	
	//check whether the reconstructed K0s in jet cone are stemming from MC gen K0s (on MCgenK0s list):__________________________________________________
	
	for(Int_t ic=0; ic<jetConeK0list->GetSize(); ++ic){    //loop over all reconstructed K0s in jet cone
 
	//for(Int_t ic=0; ic<fListK0s->GetSize(); ++ic){     //loop over all reconstructed K0s -> previous definition of reconstruction efficiency, not sure what is the better one to choose
	   
	  Int_t negDaughterpdg;
	  Int_t posDaughterpdg;
	  Int_t motherType;
	  Int_t v0Label;
	  Double_t fPtMCrecK0Match;
	  Double_t invMK0Match;
	  Double_t MCPt;
	  Int_t nnum =-1;
	  Int_t pnum =-1;
	  Bool_t fPhysicalPrimary = -1;
          Int_t MCv0PDGCode =0;
	  Double_t jetPtSmear = -1;

	  AliAODv0* v0c = dynamic_cast<AliAODv0*>(jetConeK0list->At(ic));//pointer to reconstructed K0s inside jet cone (cone is placed around reconstructed jet axis)
	  
	  //AliAODv0* v0c = dynamic_cast<AliAODv0*>(fListK0s->At(ic));//pointer to reconstructed K0s
	  if(!v0c) continue;
	  
	  Bool_t daughtercheck = DaughterTrackCheck(v0c, nnum, pnum);//check daughter tracks have proper sign
	  if(daughtercheck == kFALSE)continue;
	  
	  const AliAODTrack *trackMCNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));
	  const AliAODTrack *trackMCPos=(AliAODTrack *)(v0c->GetDaughter(pnum));
	  
	  TList *listmc = fAOD->GetList();
	 	  
	  Bool_t mclabelcheck = MCLabelCheck(v0c, kK0, trackMCNeg, trackMCPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode);
	  		  
	  if(mclabelcheck == kFALSE)continue;  //requirements for rec. V0 associated to MC true primary particle

	  for(Int_t it=0; it<fListMCgenK0s->GetSize(); ++it){                                    // loop over MC generated K0s in event, check whether associated MC particle is part of it
	  
	  //for(Int_t it=0; it<fListMCgenK0sCone->GetSize(); ++it){//belongs to previous definition of rec. eff. of V0s within jet cone  

	    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    //AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0sCone->At(it));
	    AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0s->At(it));
	    if(!mcp0) continue;
	    
	    Bool_t particleMatching = IsParticleMatching(mcp0, v0Label);
	    
	    if(particleMatching == kFALSE)continue;                                              //if reconstructed V0 particle doesn't match to the associated MC particle go to next stack entry	    
	    CalculateInvMass(v0c, kK0, invMK0Match, fPtMCrecK0Match);
	    
	    Double_t fPtMCgenK0s    = mcp0->Pt();
	    
	    fh3MCrecK0Cone->Fill(jetPt,invMK0Match,fPtMCgenK0s);                                 //fill matching rec. K0s in 3D histogram

	    SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);           //jetPt, cent, jetRadius, ptmintrack, &jetPtSmear   	 
	    
	    fh3MCrecK0ConeSmear->Fill(jetPtSmear,invMK0Match,fPtMCgenK0s);    //fill matching rec. K0s in 3D histogram, jet pT smeared according to deltaptjet distribution width  
  

	  } // end MCgenK0s / MCgenK0sCone loop
  
	  //___________
	  //check the K0s daughters contamination of the jet tracks:
	  	
	  TClonesArray *stackMC = 0x0;
	  
	  for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all particles in the jet
	      
	    AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//all tracks in jet cone	
	    if(!trackVP)continue;
	    AliAODTrack *tr = dynamic_cast<AliAODTrack*> (trackVP);                   //fetch one jet track from the TList
	    if(!tr)continue;
	      	 
	    //get MC label information
	    TList *mclist = fAOD->GetList();                                           //fetch the MC stack
	    if(!mclist)continue;

	    stackMC = (TClonesArray*)mclist->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
	    if (!stackMC) {Printf("ERROR: stack not available");}
	    else {
	      	      
	      Int_t particleLabel = TMath::Abs(tr->GetLabel());                       //fetch jet track label in MC stack
	      if(!tr)continue;
	      //v0c is pointer to K0s candidate, is fetched already above, here it is just checked again whether daughters are properly ordered by their charge   

	      Bool_t daughterchecks = DaughterTrackCheck(v0c, nnum, pnum);
		
	      if(daughterchecks == kFALSE)continue;                                   //make sure that daughters are properly ordered

	      const AliAODTrack *trackNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));    //fetch v0 daughters of reconstructed K0s
	      const AliAODTrack *trackPos=(AliAODTrack *)(v0c->GetDaughter(pnum));
	      
	      if(!trackNeg)continue;
	      if(!trackPos)continue;

	      Int_t negAssLabel = TMath::Abs(trackNeg->GetLabel());                   //negative (reconstructed) charged track label in MC stack
	      Int_t posAssLabel = TMath::Abs(trackPos->GetLabel());                   //positive (reconstructed) charged track label in MC stack

	    
	      if(particleLabel == posAssLabel){                                       //check whether jet track and each of the rec. K0s daughters have same MC stack label -> are identical
                AliAODMCParticle* mctrackPos = dynamic_cast<AliAODMCParticle*>(stackMC->At(posAssLabel));
		if(!mctrackPos) continue;
		Double_t trackPosPt = mctrackPos->Pt();
		Double_t trackPosEta = mctrackPos->Eta();
		fh3SecContinCone->Fill(jetPt, trackPosPt, trackPosEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo 
	      
	      if(particleLabel == negAssLabel){
		AliAODMCParticle* mctrackNeg = dynamic_cast<AliAODMCParticle*>(stackMC->At(negAssLabel));
		if(!mctrackNeg) continue;
		Double_t trackNegPt = mctrackNeg->Pt();
		Double_t trackNegEta = mctrackNeg->Eta();
		fh3SecContinCone->Fill(jetPt, trackNegPt, trackNegEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo
	    }
	  }
	  	  
	    
	  //_______________
	  
	  
	} //end rec-K0-in-cone loop
	
	//________________________________________________________________________________________________________________________________________________________
	  


	delete fListMCgenK0sCone;
	
	
      }//end fAnalysisMC
      
      delete jetConeK0list;      
      delete jetPerpConeK0list;
      delete jetPerpConeLalist;
      delete jetPerpConeALalist;
      delete jetMedianConeK0list;
      delete jetMedianConeLalist;
      delete jetMedianConeALalist;

      //---------------La--------------------------------------------------------------------------------------------------------------------------------------------
      
      // fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
      
      for(Int_t it=0; it<fListLa->GetSize(); ++it){ // loop all La 
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLa->At(it));
	if(!v0) continue;
	
	Double_t v0Mom[3];
	v0->PxPyPz(v0Mom);
	TVector3 v0MomVect(v0Mom);

	Double_t dPhiJetLa = (jet->MomentumVector()->Vect()).DeltaPhi(v0MomVect);
	
        Double_t invMLa =0;
        Double_t trackPt=0;

        CalculateInvMass(v0, kLambda, invMLa, trackPt); //function to calculate invMass with TLorentzVector class
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	//if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	fFFHistosIMLaJet->FillFF(trackPt, invMLa, jetPt, incrementJetPt);
 	fFFHistosPhiCorrIMLa->FillPhiCorr(trackPt,TVector2::Phi_0_2pi(dPhiJetLa),invMLa);
	
	if(dPhiJetLa<fh1dPhiJetLa->GetXaxis()->GetXmin()) dPhiJetLa += 2*TMath::Pi();
	fh1dPhiJetLa->Fill(dPhiJetLa);
      }

      if(fListLa->GetSize() == 0){ // no La: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMLaJet->FillFF(-1, -1, jetPt, incrementJetPt);
      }
	
  
      // ____fetch rec. Lambdas in cone around jet axis_______________________________________________________________________________________
      
      TList* jetConeLalist = new TList();
      Double_t sumPtLa     = 0.;
      Bool_t isBadJetLa    = kFALSE; // dummy, do not use

      GetTracksInCone(fListLa, jetConeLalist, jet, GetFFRadius(), sumPtLa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetLa);//method inherited from FF

      if(fDebug>2)Printf("%s:%d nLa total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetConeLalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetConeLalist->GetSize(); ++it){ // loop La in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeLalist->At(it));
	if(!v0) continue;                                                                                                                                                                                                         
	Double_t invMLa =0;
	Double_t trackPt=0;
	
        CalculateInvMass(v0, kLambda, invMLa, trackPt); //function to calculate invMass with TLorentzVector class
	
	Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	if(fAnalysisMC){
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE){fh1FFIMLaConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	}

	fFFHistosIMLaCone->FillFF(trackPt, invMLa, jetPt, incrementJetPt);
      }

      if(jetConeLalist->GetSize() == 0){ // no La: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMLaCone->FillFF(-1, -1, jetPt, incrementJetPt);

	if(fAnalysisMC){ 
	  Double_t jetPtSmear;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE)fh1FFIMLaConeSmear->Fill(jetPtSmear);}

      }
      
      if(fAnalysisMC){
	
	//____fetch MC generated Lambdas in cone around jet axis__(particles can stem from fragmentation but also from underlying event)_____________
	
	Double_t sumPtMCgenLa      = 0.;
	Bool_t isBadJetMCgenLa  = kFALSE; // dummy, do not use 
	
	//sampling MC gen. Lambdas in cone around reconstructed jet axis      
	fListMCgenLaCone = new TList(); 
	
	GetTracksInCone(fListMCgenLa, fListMCgenLaCone, jet, GetFFRadius(), sumPtMCgenLa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMCgenLa);//fetch MC generated Lambdas in cone of resolution parameter R around jet axis 
	
	if(fDebug>2)Printf("%s:%d nMCgenLa in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,fListMCgenLaCone->GetEntries(),GetFFRadius());
	
	for(Int_t it=0; it<fListMCgenLaCone->GetSize(); ++it){ // loop MC generated La in cone around jet axis
	  
	  AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLaCone->At(it));
	  if(!mcp0) continue;
	  
	  //Double_t fRapMCgenLa   = MyRapidity(mcp0->E(),mcp0->Pz());
	  Double_t fEtaMCgenLa   = mcp0->Eta();
	  Double_t fPtMCgenLa    = mcp0->Pt();
    
	  fh2MCgenLaCone->Fill(jetPt,fPtMCgenLa);
	  fh2MCEtagenLaCone->Fill(jetPt,fEtaMCgenLa);
	}
	
	
	//check whether the reconstructed La are stemming from MC gen La on fListMCgenLa List:__________________________________________________

	for(Int_t ic=0; ic<jetConeLalist->GetSize(); ++ic){//loop over all reconstructed La within jet cone, new definition

	  //for(Int_t ic=0; ic<fListLa->GetSize(); ++ic){//old definition
	  
	  Int_t negDaughterpdg;
	  Int_t posDaughterpdg;
	  Int_t motherType;
	  Int_t v0Label;
	  Double_t fPtMCrecLaMatch;
	  Double_t invMLaMatch;
	  Double_t MCPt;
	  Int_t nnum;
	  Int_t pnum;
	  Bool_t fPhysicalPrimary = -1;
	  Int_t MCv0PDGCode =0;
	  Double_t jetPtSmear = -1;

	  AliAODv0* v0c = dynamic_cast<AliAODv0*>(jetConeLalist->At(ic));//new definition

	  //AliAODv0* v0c = dynamic_cast<AliAODv0*>(fListLa->At(ic));//old definition
	  if(!v0c) continue;
	  
	  Bool_t daughtercheck = DaughterTrackCheck(v0c, nnum, pnum);
	  if(daughtercheck == kFALSE)continue;
	  
	  const AliAODTrack *trackMCNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));
	  const AliAODTrack *trackMCPos=(AliAODTrack *)(v0c->GetDaughter(pnum));	

	  TList *listmc = fAOD->GetList();
	  
	  Bool_t mclabelcheck = MCLabelCheck(v0c, kLambda, trackMCNeg, trackMCPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode);

	  if(mclabelcheck == kFALSE)continue;
	  	  
	  for(Int_t it=0; it<fListMCgenLa->GetSize(); ++it){//new definition                                  // loop over MC generated K0s in cone around jet axis

	    // for(Int_t it=0; it<fListMCgenLaCone->GetSize(); ++it){//old definition                                  // loop over MC generated La in cone around jet axis

	    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    
	    AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLa->At(it));//new definition
	    //AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLaCone->At(it));//old definition
	    
	    if(!mcp0) continue;
	    
	    Bool_t particleMatching = IsParticleMatching(mcp0, v0Label);
	    	    
	    
	    if(particleMatching == kFALSE)continue; //particle doesn't match on any associated MC gen particle in cone around rec jet axis
	    
	    CalculateInvMass(v0c, kLambda, invMLaMatch, fPtMCrecLaMatch);
	  
	    Double_t fPtMCgenLa    = mcp0->Pt();
	    
	    fh3MCrecLaCone->Fill(jetPt,invMLaMatch,fPtMCgenLa);                        //fill matching rec. K0s 3D histogram

	    SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);

	    fh3MCrecLaConeSmear->Fill(jetPtSmear,invMLaMatch,fPtMCgenLa);              //fill matching rec. Lambdas in 3D histogram, jet pT smeared according to deltaptjet distribution width     
	        

	  } // end MCgenLa loop
	  
	    //check the Lambda daughters contamination of the jet tracks://///////////////////////////////////////////////////////////////////////////////////////////
	  	
	  TClonesArray *stackMC = 0x0;
	  
	  for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all particles in the jet
	      
	    AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//all tracks in jet cone	
	    if(!trackVP)continue;
	    AliAODTrack *tr = dynamic_cast<AliAODTrack*> (trackVP);                   //fetch one jet track from the TList
	    if(!tr)continue;
	      	 
	    //get MC label information
	    TList *mclist = fAOD->GetList();                                           //fetch the MC stack
	    
	    stackMC = (TClonesArray*)mclist->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
	    if (!stackMC) {Printf("ERROR: stack not available");}
	    else {
	      	      
	      Int_t particleLabel = TMath::Abs(tr->GetLabel());                       //fetch jet track label in MC stack
	      
	      Bool_t daughterchecks = DaughterTrackCheck(v0c, nnum, pnum);
		
	      if(daughterchecks == kFALSE)continue;                                   //make sure that daughters are properly ordered

	      const AliAODTrack *trackNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));    //fetch v0 daughters of reconstructed K0s
	      const AliAODTrack *trackPos=(AliAODTrack *)(v0c->GetDaughter(pnum));
	      
	      Int_t negAssLabel = TMath::Abs(trackNeg->GetLabel());                   //negative (reconstructed) charged track label in MC stack
	      Int_t posAssLabel = TMath::Abs(trackPos->GetLabel());                   //positive (reconstructed) charged track label in MC stack

	    
	      if(particleLabel == posAssLabel){                                       //check whether jet track and each of the rec. K0s daughters have same MC stack label -> are identical

		AliAODMCParticle* mctrackPos = dynamic_cast<AliAODMCParticle*>(stackMC->At(posAssLabel));
		if(!mctrackPos) continue;

		Double_t trackPosPt = trackPos->Pt();
		Double_t trackPosEta = trackPos->Eta();
		fh3SecContinCone->Fill(jetPt, trackPosPt, trackPosEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo 
	      
	      if(particleLabel == negAssLabel){

		AliAODMCParticle* mctrackNeg = dynamic_cast<AliAODMCParticle*>(stackMC->At(negAssLabel));
		if(!mctrackNeg) continue;

		Double_t trackNegPt = trackNeg->Pt();
		Double_t trackNegEta = trackNeg->Eta();
		fh3SecContinCone->Fill(jetPt, trackNegPt, trackNegEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo
	    }
	  }

	  	  
	    
	} //end rec-La-in-cone loop
	//________________________________________________________________________________________________________________________________________________________
	
	delete fListMCgenLaCone;
	
      }//end fAnalysisMC
      
      delete jetConeLalist;
         
      
 
      //---------------ALa-----------
    
      
      // fQAJetHistosRecCutsLeading->FillJetQA( jet->Eta(), TVector2::Phi_0_2pi(jet->Phi()), jet->Pt() );
      
      for(Int_t it=0; it<fListALa->GetSize(); ++it){ // loop all ALa 
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALa->At(it));
	if(!v0) continue;
	
	Double_t v0Mom[3];
	v0->PxPyPz(v0Mom);
	TVector3 v0MomVect(v0Mom);

	Double_t dPhiJetALa = (jet->MomentumVector()->Vect()).DeltaPhi(v0MomVect);
	
        Double_t invMALa =0;
        Double_t trackPt=0;

        CalculateInvMass(v0, kAntiLambda, invMALa, trackPt); //function to calculate invMass with TLorentzVector class
	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	//if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	fFFHistosIMALaJet->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
 	fFFHistosPhiCorrIMALa->FillPhiCorr(trackPt,TVector2::Phi_0_2pi(dPhiJetALa),invMALa);
	
	if(dPhiJetALa<fh1dPhiJetALa->GetXaxis()->GetXmin()) dPhiJetALa += 2*TMath::Pi();
	fh1dPhiJetALa->Fill(dPhiJetALa);
      }

      if(fListALa->GetSize() == 0){ // no ALa: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMALaJet->FillFF(-1, -1, jetPt, incrementJetPt);
      }
	
  
      // ____fetch rec. Antilambdas in cone around jet axis_______________________________________________________________________________________
      
      TList* jetConeALalist = new TList();
      Double_t sumPtALa     = 0.;
      Bool_t isBadJetALa    = kFALSE; // dummy, do not use

      GetTracksInCone(fListALa, jetConeALalist, jet, GetFFRadius(), sumPtALa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetALa);//method inherited from FF
      
      if(fDebug>2)Printf("%s:%d nALa total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetConeALalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetConeALalist->GetSize(); ++it){ // loop ALa in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeALalist->At(it));
	if(!v0) continue;                                                                                                                                                                                                         
	Double_t invMALa =0;
	Double_t trackPt=0;
	
        CalculateInvMass(v0, kAntiLambda, invMALa, trackPt); //function to calculate invMass with TLorentzVector class
	
	Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;

	if(fAnalysisMC){    //jet pt smearing study for Antilambdas
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE){fh1FFIMALaConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	}
	
	fFFHistosIMALaCone->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
      }

      if(jetConeALalist->GetSize() == 0){ // no ALa: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	fFFHistosIMALaCone->FillFF(-1, -1, jetPt, incrementJetPt);

	if(fAnalysisMC){ 
	  Double_t jetPtSmear;  
	  SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);  
	  if(incrementJetPt == kTRUE)fh1FFIMALaConeSmear->Fill(jetPtSmear);}

      }
      
      if(fAnalysisMC){
	
	//____fetch MC generated Antilambdas in cone around jet axis__(particles can stem from fragmentation but also from underlying event)_____________
	
	Double_t sumPtMCgenALa      = 0.;
	Bool_t isBadJetMCgenALa  = kFALSE; // dummy, do not use 
	
	//sampling MC gen Antilambdas in cone around reconstructed jet axis      
	fListMCgenALaCone = new TList(); 
	
	GetTracksInCone(fListMCgenALa, fListMCgenALaCone, jet, GetFFRadius(), sumPtMCgenALa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMCgenALa);//MC generated K0s in cone around jet axis 
	
	if(fDebug>2)Printf("%s:%d nMCgenALa in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,fListMCgenALaCone->GetEntries(),GetFFRadius());
	
	for(Int_t it=0; it<fListMCgenALaCone->GetSize(); ++it){ // loop MC generated La in cone around jet axis
	  
	  AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenALaCone->At(it));
	  if(!mcp0) continue;
	  
	  //Double_t fRapMCgenALa   = MyRapidity(mcp0->E(),mcp0->Pz());
	  Double_t fEtaMCgenALa   = mcp0->Eta();
	  Double_t fPtMCgenALa    = mcp0->Pt();
    
	  fh2MCgenALaCone->Fill(jetPt,fPtMCgenALa);
	  fh2MCEtagenALaCone->Fill(jetPt,fEtaMCgenALa);
	}
	
	
	//check whether the reconstructed ALa are stemming from MC gen ALa on MCgenALa List:__________________________________________________

	for(Int_t ic=0; ic<jetConeALalist->GetSize(); ++ic){//loop over all reconstructed ALa
   
	  Int_t negDaughterpdg;
	  Int_t posDaughterpdg;
	  Int_t motherType;
	  Int_t v0Label;
	  Double_t fPtMCrecALaMatch;
	  Double_t invMALaMatch;
	  Double_t MCPt;
	  Int_t nnum;
	  Int_t pnum;
	  Bool_t fPhysicalPrimary = -1;
	  Int_t MCv0PDGCode =0;
	  Double_t jetPtSmear = -1;
	  
	  AliAODv0* v0c = dynamic_cast<AliAODv0*>(jetConeALalist->At(ic));
	  if(!v0c) continue;
	  
	  Bool_t daughtercheck = DaughterTrackCheck(v0c, nnum, pnum);
	  if(daughtercheck == kFALSE)continue;
	  
	  const AliAODTrack *trackMCNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));
	  const AliAODTrack *trackMCPos=(AliAODTrack *)(v0c->GetDaughter(pnum));	

	  TList *listmc = fAOD->GetList();
  	  if(!listmc)continue;

	  Bool_t mclabelcheck = MCLabelCheck(v0c, kAntiLambda, trackMCNeg, trackMCPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode);

	  if(mclabelcheck == kFALSE)continue;
	 
	  for(Int_t it=0; it<fListMCgenALa->GetSize(); ++it){                                  // loop over MC generated Antilambdas in cone around jet axis

	    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	    AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenALa->At(it));
	    if(!mcp0) continue;
	    
	    Bool_t particleMatching = IsParticleMatching(mcp0, v0Label);
    	    
	    if(particleMatching == kFALSE)continue; //particle doesn't match on any associated MC gen particle in cone around rec jet axis
	    
	    CalculateInvMass(v0c, kAntiLambda, invMALaMatch, fPtMCrecALaMatch);
	  
	    Double_t fPtMCgenALa  = mcp0->Pt();
	    
	    fh3MCrecALaCone->Fill(jetPt,invMALaMatch,fPtMCgenALa);                          //fill matching rec. K0s 3D histogram

	    SmearJetPt(jetPt,centPercent,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);
	    
	    fh3MCrecALaConeSmear->Fill(jetPtSmear,invMALaMatch,fPtMCgenALa); 
 	    

	  } // end MCgenALa loop


	   //___________
	  //check the Antilambda daughters contamination of the jet tracks:
	  	
	  TClonesArray *stackMC = 0x0;
	  
	  for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all particles in the jet
	      
	    AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//all tracks in jet cone	
	    if(!trackVP)continue;
	    AliAODTrack *tr = dynamic_cast<AliAODTrack*> (trackVP);                   //fetch one jet track from the TList
	    if(!tr)continue;
	      	 
	    //get MC label information
	    TList *mclist = fAOD->GetList();                                           //fetch the MC stack
	    if(!mclist)continue;
	    
	    stackMC = (TClonesArray*)mclist->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
	    if (!stackMC) {Printf("ERROR: stack not available");}
	    else {
	      	      
	      Int_t particleLabel = TMath::Abs(tr->GetLabel());                       //fetch jet track label in MC stack
	      
	      Bool_t daughterchecks = DaughterTrackCheck(v0c, nnum, pnum);
		
	      if(daughterchecks == kFALSE)continue;                                   //make sure that daughters are properly ordered

	      const AliAODTrack *trackNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));    //fetch v0 daughters of reconstructed K0s
	      const AliAODTrack *trackPos=(AliAODTrack *)(v0c->GetDaughter(pnum));
	      if(!trackPos)continue;
	      if(!trackNeg)continue; 

	      Int_t negAssLabel = TMath::Abs(trackNeg->GetLabel());                   //negative (reconstructed) charged track label in MC stack
	      Int_t posAssLabel = TMath::Abs(trackPos->GetLabel());                   //positive (reconstructed) charged track label in MC stack

	      if(!negAssLabel)continue;
	      if(!posAssLabel)continue;
	    
	      if(particleLabel == posAssLabel){                                       //check whether jet track and each of the rec. K0s daughters have same MC stack label -> are identical
		AliAODMCParticle* mctrackPos = dynamic_cast<AliAODMCParticle*>(stackMC->At(posAssLabel));
		if(!mctrackPos) continue;

		Double_t trackPosPt = trackPos->Pt();
		Double_t trackPosEta = trackPos->Eta();
		if(!trackPosPt)continue;
		if(!trackPosEta)continue;

		fh3SecContinCone->Fill(jetPt, trackPosPt, trackPosEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo 
	      
	      if(particleLabel == negAssLabel){

		AliAODMCParticle* mctrackNeg = dynamic_cast<AliAODMCParticle*>(stackMC->At(negAssLabel));
		if(!mctrackNeg) continue;

		Double_t trackNegPt = trackNeg->Pt();
		Double_t trackNegEta = trackNeg->Eta();
		
		if(!trackNegPt)continue;
		if(!trackNegEta)continue;

		fh3SecContinCone->Fill(jetPt, trackNegPt, trackNegEta);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo
	    }
	  }
	  
	} //end rec-ALa-in-cone loop
	//________________________________________________________________________________________________________________________________________________________
	
	delete fListMCgenALaCone;
	
      }//end fAnalysisMC
      
      delete jetConeALalist;
      delete jettracklist; //had been initialised at jet loop beginning


      }//end of if 'leading' or 'all jet' requirement
  }//end of jet loop

  


  fTracksRecCuts->Clear();
  fJetsRecCuts->Clear();
  fBckgJetsRec->Clear();
  fListK0s->Clear();
  fListLa->Clear();
  fListALa->Clear();
  fListFeeddownLaCand->Clear();
  fListFeeddownALaCand->Clear();
  jetConeFDLalist->Clear();
  jetConeFDALalist->Clear();
  fListMCgenK0s->Clear();
  fListMCgenLa->Clear();
  fListMCgenALa->Clear();
  
  //Post output data.
  PostData(1, fCommonHistList);    
}

// ____________________________________________________________________________________________
void AliAnalysisTaskJetChem::SetProperties(TH3F* h,const char* x, const char* y, const char* z)
{
  //Set properties of histos (x,y and z title)

  h->SetXTitle(x);
  h->SetYTitle(y);
  h->SetZTitle(z);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->GetZaxis()->SetTitleColor(1);
}


//________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetChem::AcceptBetheBloch(AliAODv0 *v0, AliPIDResponse *PIDResponse, const Int_t particletype) //dont use for MC Analysis
{ 
  
	Int_t nnum = 1; 
	Int_t pnum = 0;
	//---
	const AliAODTrack *ntracktest=(AliAODTrack *)v0->GetDaughter(nnum); 
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}
	
	const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
	
	//Check if both tracks are available
	if (!trackPos || !trackNeg) {
	  Printf("strange analysis::UserExec:: Error:Could not retrieve one of the daughter tracks\n");
	  return kFALSE;
	}
	
	//remove like sign V0s
	if ( trackPos->Charge() == trackNeg->Charge() ){
	  //if(fDebug>1) Printf("%s:%d found like-sign V0", (char*)__FILE__,__LINE__);
	  return kFALSE;
	  }  
	//--

        Double_t nsig_p = 0; //number of sigmas that positive daughter track has got in TPC pid information
        Double_t nsig_n = 0;

        const AliAODPid *pid_p=trackPos->GetDetPid();  // returns fDetPID, more detailed or detector specific pid information
        const AliAODPid *pid_n=trackNeg->GetDetPid();

	if(!pid_p)return kFALSE;
	if(!pid_n)return kFALSE;
	
        if (pid_p)
	  {
	    if(particletype == 1) //PID cut on positive charged Lambda daughters (only those with pt < 1 GeV/c)
	      {	
	
		nsig_p=PIDResponse->NumberOfSigmasTPC(trackPos,AliPID::kProton);
		Double_t protonPt = trackPos->Pt();
		if ((TMath::Abs(nsig_p) >= fCutBetheBloch) && (fCutBetheBloch >0) && (protonPt < 1)) return kFALSE;
                
	      }
	    
	    
	  }
	
        if (pid_n)
	  {
	    if(particletype == 2)
	      {	
		nsig_n=PIDResponse->NumberOfSigmasTPC(trackNeg,AliPID::kProton);
		Double_t antiprotonPt = trackNeg->Pt();
	        if ((TMath::Abs(nsig_n) >= fCutBetheBloch) && (fCutBetheBloch >0) && (antiprotonPt < 1)) return kFALSE;
	      }
	    	      
	  }

        return kTRUE;
}

//___________________________________________________________________
Bool_t AliAnalysisTaskJetChem::IsK0InvMass(const Double_t mass) const
{
  // K0 mass ? Use FF histo limits
  
  if(fFFIMInvMMin <= mass && mass < fFFIMInvMMax) return kTRUE;

  return kFALSE;
}
//___________________________________________________________________
Bool_t AliAnalysisTaskJetChem::IsLaInvMass(const Double_t mass) const
{
  // La mass ? Use FF histo limits

  
  if(fFFIMLaInvMMin <= mass && mass < fFFIMLaInvMMax) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Int_t AliAnalysisTaskJetChem::GetListOfV0s(TList *list, const Int_t type, const Int_t particletype, AliAODVertex* primVertex, AliAODEvent* aod)
{
  // fill list of V0s selected according to type
  
  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }
  
  if(type==kTrackUndef) return 0;

  if(!primVertex) return 0;

  Double_t lPrimaryVtxPosition[3];
  Double_t lV0Position[3];
  lPrimaryVtxPosition[0] = primVertex->GetX();
  lPrimaryVtxPosition[1] = primVertex->GetY();
  lPrimaryVtxPosition[2] = primVertex->GetZ();

  for(int i=0; i<aod->GetNumberOfV0s(); i++){ // loop over V0s
    

    AliAODv0* v0 = aod->GetV0(i);
  
    if(!v0)
      {
	std::cout << std::endl
		  << "Warning in AliAnalysisTaskJetChem::GetListOfV0s:" << std::endl
		  << "v0 = " << v0 << std::endl;
	continue;
      }

    Bool_t isOnFly = v0->GetOnFlyStatus();
     
    if(!isOnFly &&  (type == kOnFly || type == kOnFlyPID || type == kOnFlydEdx || type == kOnFlyPrim)) continue; 
    if( isOnFly &&  (type == kOffl  || type == kOfflPID  || type == kOffldEdx  || type == kOfflPrim))  continue; 
  
    Int_t motherType = -1;
     //Double_t v0CalcMass = 0;   //mass of MC v0
    Double_t MCPt = 0;         //pt of MC v0
 
    Double_t pp[3]={0,0,0}; //3-momentum positive charged track
    Double_t pm[3]={0,0,0}; //3-momentum negative charged track
    Double_t v0mom[3]={0,0,0};
     
    Double_t invM = 0;
    Double_t invMK0s=0;
    Double_t invMLa=0;
    Double_t invMALa=0;
    Double_t trackPt=0;
    Int_t nnum = -1;
    Int_t pnum = -1;

 
    Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);

    if(daughtercheck == kFALSE)continue;

    const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
    const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
 
   
     ///////////////////////////////////////////////////////////////////////////////////

    //calculate InvMass for every V0 particle assumption (Kaon=1,Lambda=2,Antilambda=3)
    switch(particletype){
    case kK0: 
      CalculateInvMass(v0, kK0, invM, trackPt); //function to calculate invMass with TLorentzVector class
      invMK0s=invM;
      break; 
    case kLambda: 
      CalculateInvMass(v0, kLambda, invM, trackPt); 
      invMLa=invM;
      break;   
    case kAntiLambda: 
      CalculateInvMass(v0, kAntiLambda, invM, trackPt); 
      invMALa=invM; 
      break;
    default: 
      std::cout<<"***NO VALID PARTICLETYPE***"<<std::endl; 
      return 0;   
    }


    /////////////////////////////////////////////////////////////
    //V0 and track Cuts:
    
    if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa)))continue; 
    
    //  Double_t PosEta = trackPos->AliAODTrack::Eta();//daughter track charge is sometimes wrong here, account for that!!!
    // Double_t NegEta = trackNeg->AliAODTrack::Eta();
      
   Double_t PosEta = trackPos->Eta();//daughter track charge is sometimes wrong here, account for that!!!
   Double_t NegEta = trackNeg->Eta();

    Double_t PosCharge = trackPos->Charge();
    Double_t NegCharge = trackNeg->Charge();

    if((trackPos->Charge() == 1) && (trackNeg->Charge() == -1)) //Fill daughters charge into histo to check if they are symmetric distributed
      { fh1PosDaughterCharge->Fill(PosCharge);
	fh1NegDaughterCharge->Fill(NegCharge);
      }

    //DistOverTotMom_in_2D___________
 
    Float_t fMassK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    Float_t fMassLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
   
   
    AliAODVertex* primVtx = fAOD->GetPrimaryVertex(); // get the primary vertex
    Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
    primVtx->GetXYZ(dPrimVtxPos);
    
    Float_t fPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0
    Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
    v0->GetSecondaryVtx(dSecVtxPos);
    Double_t dDecayPath[3];
    for (Int_t iPos = 0; iPos < 3; iPos++)
      dDecayPath[iPos] = dSecVtxPos[iPos]-dPrimVtxPos[iPos]; // vector of the V0 path
    Float_t fDecLen2D = TMath::Sqrt(dDecayPath[0]*dDecayPath[0]+dDecayPath[1]*dDecayPath[1]); //transverse path length R
    Float_t fROverPt = fDecLen2D/fPtV0; // R/pT

    Float_t fMROverPtK0s = fMassK0s*fROverPt; // m*R/pT
    Float_t fMROverPtLambda = fMassLambda*fROverPt; // m*R/pT

    //___________________
    Double_t fRap = -999;//init values
    Double_t fEta = -999;
    Double_t fV0cosPointAngle = -999;
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);

    Double_t fV0mom[3];
 
    fV0mom[0]=v0->MomV0X();
    fV0mom[1]=v0->MomV0Y();
    fV0mom[2]=v0->MomV0Z();

    Double_t fV0TotalMomentum = TMath::Sqrt(fV0mom[0]*fV0mom[0]+fV0mom[1]*fV0mom[1]+fV0mom[2]*fV0mom[2]);
    //  const Double_t K0sPDGmass = 0.497614; 
    // const Double_t LambdaPDGmass = 1.115683; 

    const Double_t K0sPDGmass = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); 
    const Double_t LambdaPDGmass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

    Double_t fDistOverTotMomK0s = 0;
    Double_t fDistOverTotMomLa = 0;
   
    //calculate proper lifetime of particles in 3D (not recommended anymore)
    
    if(particletype == kK0){
 
      fDistOverTotMomK0s = fV0DecayLength * K0sPDGmass;
      fDistOverTotMomK0s /= (fV0TotalMomentum+1e-10);
    }
  
    if((particletype == kLambda)||(particletype == kAntiLambda)){
      
      fDistOverTotMomLa = fV0DecayLength * LambdaPDGmass;
      fDistOverTotMomLa /= (fV0TotalMomentum+1e-10);
    } 
    
     //TPC cluster (not used anymore) and TPCRefit cuts
    
    //Double_t trackPosNcls = trackPos->GetTPCNcls();//Get number of clusters for positive charged tracks
    //Double_t trackNegNcls = trackNeg->GetTPCNcls();//Get number of clusters for negative charged tracks

    if(fRequireTPCRefit==kTRUE){//if kTRUE: accept only if daughter track is refitted in TPC!!
      Bool_t isPosTPCRefit = (trackPos->AliAODTrack::IsOn(AliESDtrack::kTPCrefit));
      Bool_t isNegTPCRefit = (trackNeg->AliAODTrack::IsOn(AliESDtrack::kTPCrefit));
      if (!isPosTPCRefit)continue;
      if (!isNegTPCRefit)continue;
    }
   
    if(fKinkDaughters==kFALSE){//if kFALSE: no acceptance of kink daughters
      AliAODVertex* ProdVtxDaughtersPos = (AliAODVertex*) (trackPos->AliAODTrack::GetProdVertex());
      Char_t isAcceptKinkDaughtersPos  = ProdVtxDaughtersPos->GetType();
      if(isAcceptKinkDaughtersPos==AliAODVertex::kKink)continue;
      
      AliAODVertex* ProdVtxDaughtersNeg = (AliAODVertex*) (trackNeg->AliAODTrack::GetProdVertex());
      Char_t isAcceptKinkDaughtersNeg  = ProdVtxDaughtersNeg->GetType();
      if(isAcceptKinkDaughtersNeg==AliAODVertex::kKink)continue;

    }

    Double_t fV0Radius      = -999;
    Double_t fDcaV0Daughters = v0->DcaV0Daughters();
    Double_t fDcaPosToPrimVertex = v0->DcaPosToPrimVertex();//IP of positive charged daughter
    Double_t fDcaNegToPrimVertex = v0->DcaNegToPrimVertex();//IP of negative charged daughter
    Double_t avDecayLengthK0s = 2.6844;
    Double_t avDecayLengthLa = 7.89;

    //Float_t fCTauK0s = 2.6844; // [cm] c tau of K0S
    //Float_t fCTauLambda = 7.89; // [cm] c tau of Lambda and Antilambda
    
    fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();  

    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
    
    if(particletype == kK0)         {fRap = v0->RapK0Short();
                                     fEta = v0->PseudoRapV0();}
    if(particletype == kLambda)     {fRap = v0->RapLambda();
                                     fEta = v0->PseudoRapV0();}
    if(particletype == kAntiLambda) {fRap = v0->Y(-3122);
                                     fEta = v0->PseudoRapV0();}


    //cut on 3D DistOverTotMom: (not used anymore)
    //if((particletype == kLambda)||(particletype == kAntiLambda)){if(fDistOverTotMomLa >= (fCutV0DecayMax * avDecayLengthLa))  continue;}

    //cut on K0s applied below after all other cuts for histo fill purposes..

    //cut on 2D DistOverTransMom: (recommended from Iouri)
    if((particletype == kLambda)||(particletype == kAntiLambda)){if(fMROverPtLambda > (fCutV0DecayMax * avDecayLengthLa))continue;}//fCutV0DecayMax set to 5 in AddTask macro

 //Armenteros Podolanski Plot for K0s:////////////////////////////
    
    Double_t ArmenterosAlpha=-999;  
    Double_t ArmenterosPt=-999;
    Double_t PosPl;
    Double_t NegPl;
    Double_t PosPt;
    Double_t NegPt;   
    
    if(particletype == kK0){
      
      pp[0]=v0->MomPosX();
      pp[1]=v0->MomPosY();
      pp[2]=v0->MomPosZ();
      
      pm[0]=v0->MomNegX();
      pm[1]=v0->MomNegY();
      pm[2]=v0->MomNegZ();
      
      
      v0mom[0]=v0->MomV0X();
      v0mom[1]=v0->MomV0Y();
      v0mom[2]=v0->MomV0Z();
      
      TVector3 v0Pos(pp[0],pp[1],pp[2]);
      TVector3 v0Neg(pm[0],pm[1],pm[2]);
      TVector3 v0totMom(v0mom[0], v0mom[1], v0mom[2]); //vector for tot v0 momentum
      
      PosPt = v0Pos.Perp(v0totMom);             //longitudinal momentum of positive charged daughter track
      PosPl = v0Pos.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of positive charged daughter track
	  
      NegPt = v0Neg.Perp(v0totMom);             //longitudinal momentum of negative charged daughter track
      NegPl = v0Neg.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of nergative charged daughter track
      
      ArmenterosAlpha = 1.-2./(1+(PosPl/NegPl));  
      ArmenterosPt= v0->PtArmV0();
      
    }      

    if(particletype == kK0){//only cut on K0s histos
      if(IsArmenterosSelected == 1){// Armenteros Cut to reject Lambdas contamination in K0s inv. massspectrum
	fh2ArmenterosBeforeCuts->Fill(ArmenterosAlpha,ArmenterosPt);
      }
    }
    
    //some more cuts on v0s and daughter tracks:
    
 
    if((TMath::Abs(PosEta)>fCutPostrackEta) || (TMath::Abs(NegEta)>fCutNegtrackEta))continue;   //Daughters pseudorapidity cut
    if (fV0cosPointAngle < fCutV0cosPointAngle)	continue;                                       //cospointangle cut

    //if(TMath::Abs(fRap) > fCutRap)continue;                                                     //V0 Rapidity Cut
    if(TMath::Abs(fEta) > fCutEta) continue;                                                  //V0 Eta Cut
    if (fDcaV0Daughters > fCutDcaV0Daughters)continue;
    if ((fDcaPosToPrimVertex < fCutDcaPosToPrimVertex) || (fDcaNegToPrimVertex < fCutDcaNegToPrimVertex))continue;
    if ((fV0Radius < fCutV0RadiusMin) || (fV0Radius > fCutV0RadiusMax))continue;
    
    const AliAODPid *pid_p1=trackPos->GetDetPid();
    const AliAODPid *pid_n1=trackNeg->GetDetPid();

   
      if(particletype == kLambda){
	//	if(AcceptBetheBloch(v0, fPIDResponse, 1) == kFALSE){std::cout<<"******PID cut rejects Lambda!!!************"<<std::endl;}
        if(AcceptBetheBloch(v0, fPIDResponse, 1) == kFALSE)continue;
	fh2BBLaPos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());//positive lambda daughter
	fh2BBLaNeg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());//negative lambda daughter
    
	//Double_t phi = v0->Phi();
        //Double_t massLa = v0->MassLambda();

        //printf("La: i = %d, m = %f, pT = %f, eta = %f, phi = %f\n, ",i,massLa,trackPt,fEta,phi);

      } 
      
      if(particletype == kAntiLambda){
	
	if(AcceptBetheBloch(v0, fPIDResponse, 2) == kFALSE)continue;
	fh2BBLaPos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());//positive antilambda daughter
	fh2BBLaNeg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());//negative antilambda daughter
	
      }
    
    
    //Armenteros cut on K0s:
    if(particletype == kK0){
      if(IsArmenterosSelected == 1){// Armenteros Cut to reject Lambdas contamination in K0s inv. massspectrum
	
	if((ArmenterosPt<=(TMath::Abs(fCutArmenteros*ArmenterosAlpha))) && (fCutArmenteros!=-999))continue; //Cuts out Lambda contamination in K0s histos
	fh2ArmenterosAfterCuts->Fill(ArmenterosAlpha,ArmenterosPt);  
      }
    }

    //not used anymore in 3D, z component of total momentum has bad resolution, cut instead in 2D and use pT
    //Proper Lifetime Cut: DecayLength3D * PDGmass / |p_tot| < 3*2.68cm (ctau(betagamma=1))  ;  |p|/mass = beta*gamma
    //////////////////////////////////////////////
    
    //cut on 3D DistOverTotMom
    /*  if(particletype == kK0){
      
      fh2ProperLifetimeK0sVsPtBeforeCut->Fill(trackPt,fDistOverTotMomK0s); //fill these histos after all other cuts
      fh1ProperLifetimeV0BeforeCut->Fill(fDistOverTotMomK0s);
      if(fDistOverTotMomK0s >= (fCutV0DecayMax * avDecayLengthK0s))continue;
      fh1ProperLifetimeV0AfterCut->Fill(fDistOverTotMomK0s);
      fh2ProperLifetimeK0sVsPtAfterCut->Fill(trackPt,fDistOverTotMomK0s); 
    } 
    */
    
//cut on 2D DistOverTransMom
    if(particletype == kK0){//the cut on Lambdas you can find above
      
      fh2ProperLifetimeK0sVsPtBeforeCut->Fill(trackPt,fMROverPtK0s); //fill these histos after all other cuts
      fh1ProperLifetimeV0BeforeCut->Fill(fMROverPtK0s);
      if(fMROverPtK0s > (fCutV0DecayMax * avDecayLengthK0s))continue;
      fh1ProperLifetimeV0AfterCut->Fill(fMROverPtK0s);
      fh2ProperLifetimeK0sVsPtAfterCut->Fill(trackPt,fMROverPtK0s); 
      
      //Double_t phi = v0->Phi();
      // Double_t massK0s = v0->MassK0Short();
      //printf("K0S: i = %d, m = %f, pT = %f, eta = %f, phi = %f\n",i,invMK0s,trackPt,fEta,phi);
      
      //test std::cout<<" Index accepted K0s candidate in list of V0s in event: "<<i<<" m: "<<invMK0s<<" pT: "<<trackPt<<" eta: "<<fEta<<" phi: "<<v0->Phi()<<std::endl;    
      
    } 
    //MC Associated V0 particles: (reconstructed particles associated with MC truth (MC truth: true primary MC generated particle))
    
    
    if(fAnalysisMC){// begin MC part
      
      Int_t negDaughterpdg = 0;
      Int_t posDaughterpdg = 0;
      Int_t v0Label = -1;
      Bool_t fPhysicalPrimary = -1;   //v0 physical primary check
      Int_t MCv0PdgCode = 0;
      Bool_t mclabelcheck = kFALSE;

      TList *listmc = aod->GetList(); //AliAODEvent* is inherited from AliVEvent*, listmc is pointer to reconstructed event in MC list, member of AliAODEvent
      
      if(!listmc)continue;

      if((particletype == kLambda) || (particletype == kAntiLambda)){// at this point the v0 candidates already survived all V0 cuts, for the MC analysis they still have to survive the association checks in the following block
	
	//feeddown-correction for Lambda/Antilambda particles
	//feedddown comes mainly from charged and neutral Xi particles
	//feeddown from Sigma decays so quickly that it's not possible to distinguish from primary Lambdas with detector
	//feeddown for K0s from phi decays is neglectible
	//TH2F* fh2FeedDownMatrix = 0x0; //histo for feeddown already decleared above
	
	
	//first for all Lambda and Antilambda candidates____________________________________________________________________
	
	if(particletype == kLambda){
	  
	  mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
	  
	  
	  if((motherType == 3312)||(motherType == 3322)){//mother of v0 is neutral or negative Xi
	    	    
	    fListFeeddownLaCand->Add(v0); //fill TList with ass. particles, stemming from feeddown from Xi(bar) decays       	    
	  }
	}
	if(particletype == kAntiLambda){
	  mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
	  
	  if((motherType == -3312)||(motherType == -3322)){
	    fListFeeddownALaCand->Add(v0); //fill TList with ass. particles, stemming from feeddown from Xi(bar) decays	      	          
	  }
	}
      }
      
      //_only true primary particles survive the following checks_______________________________________________________________________________________________
      
      if(particletype == kK0){
	mclabelcheck = MCLabelCheck(v0, kK0, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
	if(mclabelcheck == kFALSE)continue;
      }
      if(particletype == kLambda){
	mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
	if(mclabelcheck == kFALSE)continue;
      }
      if(particletype == kAntiLambda){
	mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode);
	if(mclabelcheck == kFALSE)continue;
      }
      
      if(fPhysicalPrimary != 1)continue; //V0 candidate (K0s, Lambda or Antilambda) must be physical primary, this means there is no mother particle existing
      
    }
   

    list->Add(v0);

  }
  
  Int_t nPart=list->GetSize();
 
  return nPart;
} // end GetListOfV0s()

// -------------------------------------------------------------------------------------------------------

void AliAnalysisTaskJetChem::CalculateInvMass(AliAODv0* v0vtx, const Int_t particletype, Double_t& invM, Double_t& trackPt){ 

   //particletype:
   // * kaon = 1
   // * lambda = 2
   // * antilambda = 3

     invM    = 0;
     trackPt = 0;
   
     Double_t pp[3]={0,0,0}; //3-momentum positive charged track
     Double_t pm[3]={0,0,0}; //3-momentum negative charged track

     const Double_t massPi = 0.13957018; //better use PDG code at this point
     const Double_t massP  = 0.93827203;
    
     Double_t mass1=0;  
     Double_t mass2=0;
    
     TLorentzVector vector;  //lorentzvector V0 particle
     TLorentzVector fourmom1;//lorentzvector positive daughter
     TLorentzVector fourmom2;//lorentzvector negative daughter
  
   //--------------------------------------------------------------
   
   AliAODTrack *trackPos = (AliAODTrack *) (v0vtx->GetSecondaryVtx()->GetDaughter(0));//index 0 defined as positive charged track in AliESDFilter
  
   if( trackPos->Charge() == 1 ){
    
     pp[0]=v0vtx->MomPosX();
     pp[1]=v0vtx->MomPosY();
     pp[2]=v0vtx->MomPosZ();

     pm[0]=v0vtx->MomNegX();
     pm[1]=v0vtx->MomNegY();
     pm[2]=v0vtx->MomNegZ();
   }

   if( trackPos->Charge() == -1 ){ 
    
     pm[0]=v0vtx->MomPosX();
     pm[1]=v0vtx->MomPosY();
     pm[2]=v0vtx->MomPosZ();

     pp[0]=v0vtx->MomNegX();
     pp[1]=v0vtx->MomNegY();
     pp[2]=v0vtx->MomNegZ();
   }

   if (particletype == kK0){ // case K0s 
     mass1 = massPi;//positive particle
     mass2 = massPi;//negative particle
   } else if (particletype == kLambda){ // case Lambda
     mass1 = massP;//positive particle
     mass2 = massPi;//negative particle
   } else if (particletype == kAntiLambda){ //case AntiLambda
     mass1 = massPi;//positive particle
     mass2 = massP; //negative particle
   }
  
   fourmom1.SetXYZM(pp[0],pp[1],pp[2],mass1);//positive track
   fourmom2.SetXYZM(pm[0],pm[1],pm[2],mass2);//negative track
   vector=fourmom1 + fourmom2;
  
   invM    = vector.M(); 
   trackPt = vector.Pt();

   /*// don't apply AliAODv0 methods to get the inv. mass for the OnFly finder, since the daughter labels are sometimes switched!!!! For Offline V0 finder no problem

   if(particletype == kK0){
     std::cout << "invMK0s: " << invM <<std::endl;
     std::cout << "v0vtx->MassK0Short(): " << v0vtx->MassK0Short() << std::endl;
     std::cout << "    " <<std::endl;
     //invM = v0vtx->MassK0Short();
   }

   if(particletype == kLambda){
     std::cout << "invMLambda: " << invM <<std::endl;
     std::cout << "v0vtx->MassMassLambda(): " << v0vtx->MassLambda() << std::endl; 
     std::cout << "    " <<std::endl;
    //invM = v0vtx->MassLambda();
   }

   if(particletype == kAntiLambda){
     std::cout << "invMAntiLambda: " << invM <<std::endl;
     std::cout << "v0vtx->MassAntiLambda(): " << v0vtx->MassAntiLambda() << std::endl;
     std::cout << "    " <<std::endl;
     //invM = v0vtx->MassAntiLambda();
   }
   */

   return;
}


//_____________________________________________________________________________________
Int_t AliAnalysisTaskJetChem::GetListOfMCParticles(TList *outputlist, const Int_t particletype, AliAODEvent *mcaodevent) //(list to fill here e.g. fListMCgenK0s, particle species to search for)
{

  outputlist->Clear();

  TClonesArray *stack = 0x0;
  Double_t mcXv=0., mcYv=0., mcZv=0.;//MC vertex position
  Int_t ntrk =0;

  // get MC generated particles

  Int_t fPdgcodeCurrentPart = 0; //pdg code current particle
  Double_t fRapCurrentPart  = 0; //get rapidity
  Double_t fPtCurrentPart   = 0; //get transverse momentum
  Double_t fEtaCurrentPart = 0;  //get pseudorapidity 

  //variable for check: physical primary particle
  //Bool_t IsPhysicalPrimary = -1;
  //Int_t index = 0; //check number of injected particles
  //****************************
  // Start loop over MC particles
 
  TList *lst = mcaodevent->GetList();
  
  if(!lst){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }
  
  stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack) {
    Printf("ERROR: stack not available");
    return -1;
  }
  
  AliAODMCHeader *mcHdr=(AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHdr)return -1;

  mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ(); // position of the MC primary vertex

  
  ntrk=stack->GetEntriesFast();
  
  //if(TMath::Abs(mcZv)>10)return; //i also cut at the reconstructed particles - here i also want to cut for a second time on z vertex (?) -> could be possible bias because of resolution effects on edges of acceptance, also the case for pseudorapidity...


  for (Int_t iMc = 0; iMc < ntrk; iMc++) {                                             //loop over mc generated particles
    
   
    AliAODMCParticle *p0=(AliAODMCParticle*)stack->UncheckedAt(iMc);
    if (!p0) {
      //Printf("ERROR: particle with label %d not found in stack (mc loop)", iMc);
      continue;
    }
    fPdgcodeCurrentPart = p0->GetPdgCode();

    // Keep only K0s, Lambda and AntiLambda, Xi and Phi:
    //if ( (fPdgcodeCurrentPart != 310 ) && (fPdgcodeCurrentPart != 3122 ) && (fPdgcodeCurrentPart != -3122 ) && (fPdgcodeCurrentPart != 3312 ) && (fPdgcodeCurrentPart != -3312) && (fPdgcodeCurrentPart != -333) ) continue;
 
    
   
      //Rejection of Pythia injected particles with David Chinellatos method - not the latest method, better Method with TString from MC generator in IsInjected() function below!
      
    /*     if( (p0->GetStatus()==21) ||
	  ((p0->GetPdgCode() == 443) &&
	   (p0->GetMother() == -1) &&
	   (p0->GetDaughter(0) == (iMc))) ){ index++; } 
 
      if(p0->GetStatus()==21){std::cout<< "hello !!!!" <<std::endl;}

      std::cout<< "MC particle status:  " << p0->GetStatus() <<std::endl;
    */


      //if(index>=1){std::cout<< "MC particle status:  " << p0->GetStatus() <<std::endl;}//if first injected MC particle was found, the Status is printed out for this and every following MC particle
	

	//injected particles could be from GenBox (single high pt tracks) or jet related tracks, both generated from PYTHIA MC generator      
      
      //Check: MC particle mother
      
      //for feed-down checks
      /* //MC gen particles
        Int_t iMother = p0->GetMother(); //Motherparticle of V0 candidate (e.g. phi particle,..)
	   if(iMother >= 0){
	AliAODMCParticle *partM = (AliAODMCParticle*)stack->UncheckedAt(iMother);
	Int_t codeM = -1;
	if(partM) codeM = TMath::Abs(partM->GetPdgCode());

	
         3312    Xi-                       -3312    Xibar+          
         3322    Xi0                       -3322    Xibar0 
	

	if((codeM == 3312)||(codeM == 3322))// feeddown for Lambda coming from Xi- and Xi0



	   }
	*/
	   /*	//Check: MC gen. particle decays via 2-pion decay? -> only to be done for the rec. particles !! (-> branching ratio ~ 70 % for K0s -> pi+ pi-)
	
	Int_t daughter0Label = p0->GetDaughter(0);
	AliAODMCParticle *mcDaughter0 = (AliAODMCParticle *)stack->UncheckedAt(daughter0Label);
      if(daughter0Label >= 0)
	{daughter0Type = mcDaughter0->GetPdgCode();}
      
      Int_t daughter1Label = p0->GetDaughter(1);
      AliAODMCParticle *mcDaughter1 = (AliAODMCParticle *)stack->UncheckedAt(daughter1Label);
      
      if(daughter1Label >= 1)
	{daughter1Type = mcDaughter1->GetPdgCode();}	//requirement that daughters are pions is only done for the reconstructed V0s in GetListofV0s() below	
    }
	   */
        

      // Keep only K0s, Lambda and AntiLambda: 
      if ( (fPdgcodeCurrentPart != 310 ) && (fPdgcodeCurrentPart != 3122 ) && (fPdgcodeCurrentPart != -3122 )) continue;
      // Check: Is physical primary
      
      //do not use anymore: //IsPhysicalPrimary = p0->IsPhysicalPrimary();
      //if(!IsPhysicalPrimary)continue;

      Float_t fDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)

      // Get the distance between production point of the MC mother particle and the primary vertex
      
      Double_t dx = mcXv-p0->Xv();//mc primary vertex - mc gen. v0 vertex 
      Double_t dy = mcYv-p0->Yv();
      Double_t dz = mcZv-p0->Zv();

      Double_t fDistPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      Bool_t fPhysicalPrimary = (fDistPrimary < fDistPrimaryMax);
      
      if(!fPhysicalPrimary)continue;

      //if(fPhysicalPrimary){std::cout<<"hello**********************"<<std::endl;}
      
      /* std::cout<<"dx: "<<dx<<std::endl;
      std::cout<<"dy: "<<dy<<std::endl;
      std::cout<<"dz: "<<dz<<std::endl;

      std::cout<<"start: "<<std::endl;
      std::cout<<"mcXv: "<<mcXv<<std::endl;
      std::cout<<"mcYv: "<<mcYv<<std::endl;
      std::cout<<"mcZv: "<<mcZv<<std::endl;
            
      std::cout<<"p0->Xv(): "<<p0->Xv()<<std::endl;
      std::cout<<"p0->Yv(): "<<p0->Yv()<<std::endl;
      std::cout<<"p0->Zv(): "<<p0->Zv()<<std::endl;
      std::cout<<" "<<std::endl;
      std::cout<<"fDistPrimary"<<fDistPrimary<<std::endl;
      std::cout<<"fPhysicalPrimary"<<fPhysicalPrimary<<std::endl;
      */
 
      //Is close enough to primary vertex to be considered as primary-like?
      
      fRapCurrentPart   = MyRapidity(p0->E(),p0->Pz());
      fEtaCurrentPart   = p0->Eta();
      fPtCurrentPart    = p0->Pt();
            
      if (TMath::Abs(fEtaCurrentPart) < fCutEta){
	// if (TMath::Abs(fRapCurrentPart) > fCutRap)continue;	  //rap cut for crosschecks
	
	if(particletype == kK0){                                      //MC gen. K0s  
	  if (fPdgcodeCurrentPart==310){                       
	    outputlist->Add(p0);
	  }
	}  
    
    if(particletype == kLambda){                               	//MC gen. Lambdas
      if (fPdgcodeCurrentPart==3122)  {  
	outputlist->Add(p0);
      }
    }
    
    if(particletype == kAntiLambda){ 
      if (fPdgcodeCurrentPart==-3122) {                         //MC gen. Antilambdas
	outputlist->Add(p0);
      }
    }
  }
  
  }//end  loop over MC generated particle
  
  Int_t nMCPart=outputlist->GetSize();
  

  return nMCPart;
    
}

//---------------------------------------------------------------------------
/*
Bool_t AliAnalysisTaskJetChem::FillFeeddownMatrix(TList* fListFeeddownCand, Int_t particletype)
{

	// Define Feeddown matrix 
	Double_t lFeedDownMatrix [100][100]; 
	// FeedDownMatrix [Lambda Bin][Xi Bin];

	//Initialize entries of matrix:	
	for(Int_t ilb = 0; ilb<100; ilb++){
	  for(Int_t ixb = 0; ixb<100; ixb++){ 
	    lFeedDownMatrix[ilb][ixb]=0; //first lambda bins, xi bins
	  }
	}
}
*/
//----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetChem::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 
//----------------------------------------------------------------------------

// ________________________________________________________________________________________________________________________//function to get the MC gen. jet particles

void AliAnalysisTaskJetChem::GetTracksInCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, 
								const Double_t radius, Double_t& sumPt, const Double_t minPt, const Double_t maxPt, Bool_t& isBadPt)
{
  // fill list of tracks in cone around jet axis  

  sumPt = 0;
  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;   

  Double_t jetMom[3];
  if(!jet)return;
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

      if(maxPt>0 && track->Pt()>maxPt) isBadMaxPt = kTRUE;   // reject jets containing any track with pt larger than this value, use GetFFMaxTrackPt()
      if(minPt>0 && track->Pt()>minPt) isBadMinPt = kFALSE;  // reject jets with leading track with pt smaller than this value, use GetFFMinLTrackPt()

    }
  }

  isBadPt = kFALSE; 
  if(minPt>0 && isBadMinPt) isBadPt = kTRUE; //either the jet is bad because of too small leading track pt.. (probability to be purely combinatorial jet is too high to accept it)
  if(maxPt>0 && isBadMaxPt) isBadPt = kTRUE; //..or because of leading track with too high pt (could be fake track) 
   
  outputlist->Sort();
  
}

//____________________________________________________________________________________________________________________


void AliAnalysisTaskJetChem::GetTracksInPerpCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, 
								const Double_t radius, Double_t& sumPerpPt)
{
  // fill list of tracks in two cones around jet axis rotated in phi +/- 90 degrees
 
  Double_t jetMom[3];         //array for entries in TVector3
  Double_t perpjetplusMom[3];     //array for entries in TVector3
  Double_t perpjetnegMom[3];
 
  if(!jet)return;

  jet->PxPyPz(jetMom);   //get 3D jet momentum
  Double_t jetPerpPt = jet->Pt(); //original jet pt, invariant under rotations
  Double_t jetPhi = jet->Phi(); //original jet phi

  Double_t jetPerpposPhi = jetPhi + ((TMath::Pi())*0.5);//get new perp. jet axis phi clockwise
  Double_t jetPerpnegPhi = jetPhi - ((TMath::Pi())*0.5);//get new perp. jet axis phi counterclockwise

  TVector3 jet3mom(jetMom); //3-Vector for original rec. jet axis
 
  perpjetplusMom[0]=(TMath::Sin(jetPerpposPhi)*jetPerpPt); //x coordinate (sidewards - when looking in beam direction)
  perpjetplusMom[1]=(TMath::Cos(jetPerpposPhi)*jetPerpPt); //y coordinate (upwards - when looking in beam direction)
  perpjetplusMom[2]=jetMom[2];                          //z coordinate (along beam axis), invariant under azimuthal rotation
              
  perpjetnegMom[0]=(TMath::Sin(jetPerpnegPhi)*jetPerpPt); //x coordinate (sidewards - when looking in beam direction)
  perpjetnegMom[1]=(TMath::Cos(jetPerpnegPhi)*jetPerpPt); //y coordinate (upwards - when looking in beam direction)
  perpjetnegMom[2]=jetMom[2];                          //z coordinate (along beam axis), invariant under azimuthal rotation
    
         
  TVector3 perpjetplus3mom(perpjetplusMom);  //3-Vector for new perp. jet axis, clockwise rotated
  TVector3 perpjetneg3mom(perpjetnegMom);    //3-Vector for new perp. jet axis, counterclockwise rotated

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){ //collect V0 content in perp cone, rotated clockwise 

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack)); //inputlist is fListK0s, all reconstructed K0s in event
    if(!track){std::cout<<"K0s track not found!!!"<<std::endl; continue;}

    Double_t trackMom[3];//3-mom of V0 particle
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = perpjetplus3mom.DeltaR(track3mom);

    if(dR<radius){

      outputlist->Add(track); // output list is jetPerpConeK0list
      
      sumPerpPt += track->Pt();


    }
  }


  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){//collect V0 content in perp cone, rotated counterclockwise 

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack)); //inputlist is fListK0s, all reconstructed K0s in event
    if(!track){std::cout<<"K0s track not found!!!"<<std::endl; continue;}

    Double_t trackMom[3];//3-mom of V0 particle
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = perpjetneg3mom.DeltaR(track3mom);

    if(dR<radius){

      outputlist->Add(track); // output list is jetPerpConeK0list
      
      sumPerpPt += track->Pt();


    }
  }

  // pay attention: this list contains the double amount of V0s, found in both cones
  // before using it, devide spectra by 2!!!
  sumPerpPt = sumPerpPt*0.5; //correct to do this?

   
  outputlist->Sort();
  
}


// _______________________________________________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskJetChem::MCLabelCheck(AliAODv0* v0, Int_t particletype,const AliAODTrack* trackNeg, const AliAODTrack* trackPos, TList *listmc, Int_t& negDaughterpdg, Int_t& posDaughterpdg, Int_t& motherType, Int_t& v0Label, Double_t& MCPt, Bool_t& fPhysicalPrimary, Int_t& MCv0PDGCode){
                                
  TClonesArray *stackmc = 0x0;
  stackmc = (TClonesArray*)listmc->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
  if (!stackmc)
    {
      Printf("ERROR: stack not available");
      return kFALSE;
    }
  else
    {	   
      Int_t negAssLabel = TMath::Abs(trackNeg->GetLabel());                       //negative (reconstructed) charged track label in MC stack
      Int_t posAssLabel = TMath::Abs(trackPos->GetLabel());                       //positive (reconstructed) charged track label in MC stack
      
      //injected particle checks

      Double_t mcXv = 0;
      Double_t mcYv = 0;
      Double_t mcZv = 0;
      
      AliAODMCHeader *header=(AliAODMCHeader*)listmc->FindObject(AliAODMCHeader::StdBranchName());
      if(!header)return kFALSE;
     
      mcXv=header->GetVtxX(); mcYv=header->GetVtxY(); mcZv=header->GetVtxZ();

      Int_t trackinjected = IsTrackInjected(v0, header, stackmc); //requires AliAODTrack instead of AliVTrack

      if(trackinjected == 0){std::cout<<"HIJING track injected!!: "<<trackinjected<<std::endl;}

      //mc label checks
     
      if(negAssLabel>=0 && negAssLabel < stackmc->GetEntriesFast() && posAssLabel>=0 && posAssLabel < stackmc->GetEntriesFast()){//safety check if label has valid value of stack

	AliAODMCParticle *mcNegPart =(AliAODMCParticle*)stackmc->UncheckedAt(negAssLabel);//fetch the, with one MC truth track associated (reconstructed), negative charged track 
	v0Label = mcNegPart->GetMother();
	negDaughterpdg = mcNegPart->GetPdgCode();
	AliAODMCParticle *mcPosPart =(AliAODMCParticle*)stackmc->UncheckedAt(posAssLabel);//fetch the, with one MC truth track associated (reconstructed), positive charged track 
	Int_t v0PosLabel = mcPosPart->GetMother();                                        //get mother label of positive charged track label
	posDaughterpdg = mcPosPart->GetPdgCode();

	if(v0Label >= 0 && v0Label < stackmc->GetEntriesFast() && v0Label == v0PosLabel){//first v0 mc label check, then: check if both daughters are stemming from same particle
  
	  AliAODMCParticle *mcv0 = (AliAODMCParticle *)stackmc->UncheckedAt(v0Label);  //fetch MC ass. particle to v0 (mother of the both charged daughter tracks)
	 
	  //do not use anymore: 
	  //fPhysicalPrimary = mcv0->IsPhysicalPrimary(); 

	 
	  Float_t fDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)
	  
	  // Get the distance between production point of the MC mother particle and the primary vertex
	 
	  Double_t dx = mcXv-mcv0->Xv();//mc primary vertex - mc particle production vertex 
	  Double_t dy = mcYv-mcv0->Yv();
	  Double_t dz = mcZv-mcv0->Zv();
	  
	  Float_t fDistPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	  fPhysicalPrimary = kFALSE;//init

	  fPhysicalPrimary = (fDistPrimary < fDistPrimaryMax);

	  //if(fPhysicalPrimary == kTRUE){std::cout<<"hello*********!!!!!!!!!!!!! "<<std::endl;} 
	  
	  
	  
	  MCv0PDGCode = mcv0->GetPdgCode();
	  
	  //std::cout<<"MCv0PDGCode: "<<MCv0PDGCode<<std::endl;
	  
	  MCPt = mcv0->Pt();//for MC data, always use MC gen. pt for any pt distributions, also for the spectra, used for normalisation
	  //for feed-down checks later
	  
	  Int_t motherLabel = mcv0->GetMother();  //get mother particle label of v0 particle
	  // std::cout<<"motherLabel: "<<motherLabel<<std::endl;
	  
	  if(motherLabel >= 0 && v0Label < stackmc->GetEntriesFast())                 //do safety check for mother label
	    {
	      AliAODMCParticle *mcMother = (AliAODMCParticle *)stackmc->UncheckedAt(motherLabel);  //get mother particle
	      motherType = mcMother->GetPdgCode(); //get PDG code of mother 
	      
	      Double_t XiPt = 0.;
	      Double_t XibarPt = 0.;
	      
	      if(particletype == kLambda){
		if((motherType == 3312)||(motherType == 3322)){ //if v0 mother is Xi0 or Xi- fill MC gen. pt in FD La histogram
		  XiPt = mcMother->Pt();
		  fh1MCXiPt->Fill(XiPt);
		}     
	      }
	      if(particletype == kAntiLambda){
		if((motherType == -3312)||(motherType == -3322)){ //if v0 mother is Xibar0 or Xibar+ fill MC gen. pt in FD ALa histogram
		  XibarPt = mcMother->Pt();
		  fh1MCXibarPt->Fill(XibarPt);
		} 
		}  
	      
	    }
	  
	  //pdg code checks etc..
	  
	  if(particletype == kK0){
	       
	       if(TMath::Abs(posDaughterpdg) != 211){return kFALSE;}//one or both of the daughters are not a pion
	       if(TMath::Abs(negDaughterpdg) != 211){return kFALSE;}
	       
	       if(MCv0PDGCode != 310)  {fh1noAssociatedK0s->Fill(1.);return kFALSE;}
	     }
	     
	     if(particletype == kLambda){
	       if(MCv0PDGCode != 3122)return kFALSE;//if particle is not Antilambda, v0 is rejected
	       if(posDaughterpdg != 2212)return kFALSE;
	       if(negDaughterpdg != -211)return kFALSE;    //pdg code check for Lambda daughters
	       
	       //{if((motherType == 3312)||(motherType == 3322)){continue;}//if Xi0 and Xi- is motherparticle of Lambda, particle is rejected, pay attention, most possible Xi-, Xi0 and Omega- are not distributed physically and are much more abundant than expected by physics       //}
	     }
	     
	     if(particletype == kAntiLambda){
	       if(MCv0PDGCode != -3122)return kFALSE;
	       if(posDaughterpdg != 211)return kFALSE;
	       if(negDaughterpdg !=-2212)return kFALSE;    //pdg code check for Antilambda daughters
	       //if(fPhysicalPrimary == 1){fh1MCmotherALa->Fill(7.);}
	       
	       //{if((motherType == -3312)||(motherType == -3322)){continue;}//if bar{Xi0} and Xi+ is motherparticle of Antilambda, particle is rejected
	       //}	  
	     }
	     	
	     return kTRUE;                     //check was successful
	   }//end mc v0 label check
      }// end of stack label check
    }//end of else



  return kFALSE;                               //check wasn't successful
} 
//________________________________________________________________________________________________________________________________________________________


Bool_t AliAnalysisTaskJetChem::IsParticleMatching(const AliAODMCParticle* mcp0, const Int_t v0Label){

  const Int_t mcp0label = mcp0->GetLabel();
    
  if(v0Label == mcp0label)return kTRUE;
 
  return kFALSE;
}

//_______________________________________________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskJetChem::DaughterTrackCheck(AliAODv0* v0, Int_t& nnum, Int_t& pnum){
  

  if(v0->GetNDaughters() != 2) return kFALSE;//case v0 has more or less than 2 daughters, avoids seg. break at some AOD files                                   //reason?


  // safety check of input parameters
  if(v0 == NULL)
    {
      if(fDebug > 1){std::cout << std::endl
		<< "Warning in AliAnalysisTaskJetChem::DaughterTrackCheck:" << std::endl
			   << "v0 = " << v0 << std::endl;}

      return kFALSE;
    }
  else
    {
      //Daughters track check: its Luke Hanrattys method to check daughters charge

      nnum = 1; 
      pnum = 0;

         
      AliAODTrack *ntracktest =(AliAODTrack*)(v0->GetDaughter(nnum));
  
      if(ntracktest == NULL)
	{
	  if(fDebug > 1){std::cout << std::endl
		    << "Warning in AliAnalysisTaskJetChem::DaughterTrackCheck:" << std::endl
			 << "ntracktest = " << ntracktest << std::endl;}

	  return kFALSE;
	}

      if(ntracktest->Charge() > 0)
	{
	  nnum = 0; 
	  pnum = 1;
	}
  
      const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
      const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
  
      //Check if both tracks are available
      if (!trackPos || !trackNeg) {
	if(fDebug > 1) Printf("strange analysis::UserExec:: Error:Could not retrieve one of the daughter tracks\n");
	return kFALSE;
      }
  
    
      //remove like sign V0s
      if ( trackPos->Charge() == trackNeg->Charge() ){
	//if(fDebug>1) Printf("%s:%d found like-sign V0", (char*)__FILE__,__LINE__);
	return kFALSE;
      }  
  
      return kTRUE;
    }
}

//_______________________________________________________________________________________________________________________________________________________

Int_t AliAnalysisTaskJetChem::IsTrackInjected(AliAODv0 *v0, AliAODMCHeader *header, TClonesArray *arrayMC){//info in TString should be available from 2011 data productions on..
  
  if(!v0){std::cout << " !part " << std::endl;return 1;}
  if(!header){std::cout << " !header " << std::endl;return 1;}
  if(!arrayMC){std::cout << " !arrayMC " << std::endl;return 1;}

  Int_t lab=v0->GetLabel();
  if(lab<0) {return 1;} 
  TString bbb = GetGenerator(lab,header);
  TString empty="";
  
  // std::cout << " TString bbb: " << bbb << std::endl;

  //  std::cout << " FIRST CALL " << bbb << std::endl;
  
  while(bbb.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){return 1;}
    Int_t mother = mcpart->GetMother();
    lab=mother;
    bbb = GetGenerator(mother,header);
    std::cout << "Testing " << bbb << " "  << std::endl;
  }

  std::cout << " FINAL CALL " << bbb << std::endl;
  
  //std::transform(bbb.begin(), bbb.end(), bbb.begin(), ::tolower);   //convert TString bbb into lower case, to avoid that TString could be  written in lower or upper case

  if(bbb.Contains("ijing")){std::cout << " particle is injected!! " << std::endl; return 0;}//if TString returns something with "ijing" return this method with 0 -> select out all HIJING particles, all others return with "1"
  
 
  return 1;
}

//______________________________________________________________________
TString AliAnalysisTaskJetChem::GetGenerator(Int_t label, AliAODMCHeader* header){
   Int_t nsumpart=0;
   TList *lh=header->GetCocktailHeaders();
   Int_t nh=lh->GetEntries();
   for(Int_t i=0;i<nh;i++){
     AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
     TString genname=gh->GetName();
     Int_t npart=gh->NProduced();
     if(label>=nsumpart && label<(nsumpart+npart)) return genname;
     nsumpart+=npart;
   }
   TString empty="";
   return empty;
 }

//_________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetChem::SmearJetPt(Double_t jetPt, Int_t cent, Double_t jetRadius, Double_t ptmintrack, Double_t& jetPtSmear){	   
  
  TF1 *fsmear = new TF1("f1","[0]*exp(-1*(x-[1])*(x-[1])/(2*[2]*[2]))",-100.,100.);   //smearing according to gaussian function in between  +/- 10 GeV/c
  
  jetRadius = 0.4;
  ptmintrack = 0.15;
  cent = 10;
  //Int_t cl = 1;
  
  /*  if(cent>10) cl = 2; 
  if(cent>30) cl = 3;
  if(cent>50) cl = 4;
  */

  fsmear->SetParameters(1,0,11.19);//for 2010 PbPb jets, R=0.4, ptmintrack = 0.15 GeV/c, cent 00-10%, delta-pt width estimated via single track embedding
  //fsmear->SetParameters(1,0,3.28);//for 2010 PbPb jets, R=0.4, ptmintrack = 0.15 GeV/c, cent 50-60%, delta-pt width estimated via single track embedding
  
  //fsmear->SetParameters(1,0,4.472208);// for 2010 PbPb jets, R=0.2, ptmintrack = 0.15 GeV/c, cent 00-10%
  
  /* //delta-pt width for anti-kt jet finder:
     
  // jet cone R = 0.4
  if((cl == 1)&&(jetRadius == 0.4)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,10.178069);//(max.,mean,sigma) of gaussian, needs to be adjusted for every combination of jet cone size, centrality and min. pt constituents cut
  }  
  if((cl == 2)&&(jetRadius == 0.4)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,8.536195);
  }
  if((cl == 3)&&(jetRadius == 0.4)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,?);
  }
  if((cl == 4)&&(jetRadius == 0.4)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,5.229839);
  }
  
  // jet cone R = 0.3     
  if((cl == 1)&&(jetRadius == 0.3)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,7.145967);
  }
  if((cl == 2)&&(jetRadius == 0.3)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,5.844796);
  }
  if((cl == 3)&&(jetRadius == 0.3)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,?);
  }
  if((cl == 4)&&(jetRadius == 0.3)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,3.630751);
  }
  
  // jet cone R = 0.2
  if((cl == 1)&&(jetRadius == 0.2)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,4.472208);
  }
  if((cl == 2)&&(jetRadius == 0.2)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,3.543938);
  }
  if((cl == 3)&&(jetRadius == 0.2)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,?);
  }
  if((cl == 4)&&(jetRadius == 0.2)&&(ptmintrack == 0.15)){
  fsmear->SetParameters(1,0,1.037476);
  }
  
  */
  
  Double_t r = fsmear->GetRandom();
  jetPtSmear = jetPt + r;
  
  //  std::cout<<"jetPt: "<<jetPt<<std::endl;
  //  std::cout<<"jetPtSmear: "<<jetPtSmear<<std::endl;
  //  std::cout<<"r: "<<r<<std::endl;
  
  delete fsmear;
  return jetPtSmear;
}


// _______________________________________________________________________________________________________________________
AliAODJet* AliAnalysisTaskJetChem::GetMedianCluster()
{
  // fill tracks from bckgCluster branch, 
  // using cluster with median density (odd number of clusters) 
  // or picking randomly one of the two closest to median (even number)
  
  Int_t nBckgClusters = fBckgJetsRec->GetEntries(); // not 'recCuts': use all clusters in full eta range

  if(nBckgClusters<3) return 0; // need at least 3 clusters (skipping 2 highest)

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

    //Int_t medianIndex = indices[(Int_t) (0.5*(nBckgClusters-1))];
    Int_t medianIndex = indices[(Int_t) (0.5*(nBckgClusters+1))];

    medianCluster = (AliAODJet*)(fBckgJetsRec->At(medianIndex));
    
    Double_t clusterPt = medianCluster->Pt();
    Double_t area      = medianCluster->EffectiveAreaCharged();
    
    if(area>0) medianDensity = clusterPt/area;
  }
  else{

    //Int_t medianIndex1 = indices[(Int_t) (0.5*nBckgClusters-1)];
    //Int_t medianIndex2 = indices[(Int_t) (0.5*nBckgClusters)];

    Int_t medianIndex1 = indices[(Int_t) (0.5*nBckgClusters)];
    Int_t medianIndex2 = indices[(Int_t) (0.5*nBckgClusters+1)];

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
    
    medianCluster = ( (gRandom->Rndm()>0.5) ? medianCluster1 : medianCluster2 );  // select one randomly to avoid adding areas
  }
    
  delete[] bgrDensity;
  delete[] indices; 

  return medianCluster;
}    
