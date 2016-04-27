/*************************************************************************
 *                                                                       *
 *                                                                       *
 *      Task for Jet Chemistry Analysis in PWG-JE Jet Task Force Train   *
 *    Analysis of K0s, Lambda and Antilambda with and without Jetevents  *
 *                                                                       *
 *************************************************************************/

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
//Task for K0s, Lambda and Antilambda analysis in jets in Pb-Pb collisions
//Author: Alice Zimmermann (zimmermann@physi.uni-heidelberg.de)
  

/* $Id: */

#include "Riostream.h"
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
#include "TProfile.h"
#include "THnSparse.h"
#include <algorithm>
#include <string> 
#include "AliAnalysisHelperJetTasks.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h" 
#include "AliAODHeader.h" 
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
#include "TRandom3.h"

ClassImp(AliAnalysisTaskJetChem)

using std::cout;
using std::endl;


//____________________________________________________________________________
AliAnalysisTaskJetChem::AliAnalysisTaskJetChem()
   : AliAnalysisTaskFragmentationFunction()
 
   ,fRandom(0)
   ,fMatchMode(0)
   ,fIsNJEventEmb(0)
   ,fAnalysisMC(0)
   ,fDeltaVertexZ(0)
   ,fCutjetEta(0)
   ,fCuttrackNegNcls(0)
   ,fCuttrackPosNcls(0)
   ,fCutPostrackRap(0)
   ,fCutNegtrackRap(0)
   ,fCutRap(0)
   ,fCutPostrackEta(0)
   ,fCutNegtrackEta(0)
   ,fCutEta(0)
   ,fCutK0cosPointAngle(0)
   ,fCutLacosPointAngle(0)
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
   ,fCutFractionPtEmbedded(0)
   ,fCutDeltaREmbedded(0)
  //,fBranchEmbeddedJets("")  
   ,fK0Type(0)  
   ,fFilterMaskK0(0)
   ,jettracklist(0)
   ,jetConeK0list(0)
   ,jetConeLalist(0)
   ,jetConeALalist(0)
   ,jetConeK0Emblist(0)
   ,jetConeLaEmblist(0)
   ,jetConeALaEmblist(0) 
   ,jetConeK0EmbStlist(0)
   ,jetConeLaEmbStlist(0)
   ,jetConeALaEmbStlist(0)

   ,jetConeK0EmbMClist(0)
   ,jetConeLaEmbMClist(0)
   ,jetConeALaEmbMClist(0)

   ,jetPerpConeK0list(0)
   ,jetPerpRecCutslist(0)
   ,jetPerpConeK0Emblist(0)
   ,jetPerpConeLalist(0)
   ,jetPerpConeLaEmblist(0)
   ,jetPerpConeALalist(0)
   ,jetPerpConeALaEmblist(0)
   ,jetMedianConeK0list(0)
   ,jetMedianConeLalist(0)
   ,jetMedianConeALalist(0)
   ,jetMedianRecCutslist(0)
   ,fListK0sRC(0)
   ,fListLaRC(0)
   ,fListALaRC(0)
   ,fTracksRecCutsRC(0)
   ,fTracksRecBckgCuts(0) 
   ,fTracksPerpCone(0)
   ,fListK0s(0)
   ,fListK0sMC(0)
   ,fListK0sStandard(0)
   ,fPIDResponse(0)
   ,fV0QAK0(0)
   ,fFFHistosRecCutsK0Evt(0)      
  //,fFFHistosIMK0AllEvt(0)        
  //,fFFHistosIMK0Jet(0)           
  //,fFFHistosIMK0Cone(0)
   ,fLaType(0) 
   ,fFilterMaskLa(0)
   ,fListLa(0)
   ,fListLaMC(0)
   ,fListLaStandard(0)
  // ,fFFHistosIMLaAllEvt(0)        
  // ,fFFHistosIMLaJet(0)           
  //,fFFHistosIMLaCone(0)
   ,fALaType(0) 
   ,fFilterMaskALa(0)
   ,fListALa(0)
   ,fListALaMC(0)
   ,fListALaStandard(0)
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
   ,fUseExtraTracks(0)
   ,fUseExtraJetPt(0)
   ,fUseEmbeddedJetPt(0)

   ,fUseStandard(0)
  // ,fFFHistosIMALaAllEvt(0)        
  // ,fFFHistosIMALaJet(0)           
  // ,fFFHistosIMALaCone(0)
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
   ,fh1nEmbeddedJets(0)
   ,fh1nGenJets(0)
   ,fh1IndexEmbedded(0)
   ,fh1IndexEmbeddedMC(0)
   ,fh1PtEmbBeforeMatch(0)
   ,fh1PtEmbExtraOnly(0)
   ,fh1PtEmbReject(0)
   ,fh2PtEtaEmbReject(0)
   ,fh1PtEmbAfterMatch(0)
   ,fh1FractionPtEmbedded(0)
   ,fh1DeltaREmbedded(0)
   ,fh1FractionPtEmbeddedMC(0)
   ,fh2FractionPtVsEmbeddedJetPtMC(0)
   ,fh1DeltaREmbeddedMC(0)
   ,fh1JetPtEmbGenAfterMatch(0)
   ,fh2TracksPerpCone(0)
   ,fh1PerpCone(0)
  //,fh1V0JetPt(0)
   ,fh1V0PtCandidate(0)
   ,fh1IMK0Cone(0)
   ,fh1IMLaCone(0)
   ,fh1IMALaCone(0)
   ,fh1IMK0EmbCone(0)
   ,fh1IMLaEmbCone(0)
   ,fh1IMALaEmbCone(0)
   ,fh2FFJetTrackEta(0)
  //,fh1trackPosNCls(0)           
  //,fh1trackNegNCls(0)   
   ,fh1trackPosRap(0)            
   ,fh1trackNegRap(0)          
  //,fh1V0Rap(0)        
   ,fh1trackPosEta(0)            
   ,fh1trackNegEta(0)          
   ,fh1V0Eta(0)
  //,fh1V0totMom(0)
  ,fh1CosPointAngle(0)                       
   ,fh1DecayLengthV0(0)    
   ,fh2ProperLifetimeK0sVsPtBeforeCut(0)    
   ,fh2ProperLifetimeK0sVsPtAfterCut(0)
   ,fh1V0Radius(0)     
   ,fh1DcaV0Daughters(0)        
   ,fh1DcaPosToPrimVertex(0)   
   ,fh1DcaNegToPrimVertex(0)    
   ,fh2ArmenterosBeforeCuts(0)
   ,fh2ArmenterosAfterCuts(0)
   ,fh2BBLaPos(0)
   ,fh2BBLaNeg(0)
   ,fh1PosDaughterCharge(0)
   ,fh1NegDaughterCharge(0)
   ,fh1PtMCK0s(0)
   ,fh1PtMCLa(0)
   ,fh1PtMCALa(0)
   ,fh1EtaK0s(0)
   ,fh1EtaLa(0)
   ,fh1EtaALa(0)
   ,fh1RC(0)
   ,fh1RCBiasK0(0)
   ,fh1RCBiasLa(0)
   ,fh1RCBiasALa(0)
   ,fh1MCC(0)
   ,fh1OC(0)
   ,fh1NJ(0)
   ,fh1NJEmbEvt(0)
   ,fh1BckgJets(0)
   ,fh1BckgJetsPtBias(0)
   ,fhnInvMassEtaTrackPtK0s(0)
   ,fhnInvMassEtaTrackPtLa(0)
   ,fhnInvMassEtaTrackPtALa(0)
   ,fh1TrackMultCone(0)
   ,fh2TrackMultCone(0)
   ,fhnNJK0(0)
   ,fhnNJLa(0)
   ,fhnNJALa(0)
  //,fh2ChTracksNJ(0)
  //,fh2ChTracksRC(0)
     // ,fh2ChTracksOC(0)
  //,fh2ChTracksMCC(0) 
  //,fh2ChTracksPC(0)  
  ,fhnMCgenK0Cone(0)
  ,fhnMCgenLaCone(0)
  // ,fh2MCgenALaCone(0) 
  //,fh2MCEtagenK0Cone(0)
  //,fh2MCEtagenLaCone(0)
  // ,fh2MCEtagenALaCone(0)
  /*  ,fh2CorrHijingLaProton(0)
   ,fh2CorrInjectLaProton(0)
   ,fh2CorrHijingALaAProton(0)
   ,fh2CorrInjectALaAProton(0)*/
  ,fh1IMK0ConeSmear(0)
  ,fh1IMLaConeSmear(0)
  ,fh1IMALaConeSmear(0)
  ,fh2MC2K0Cone(0)
  ,fh2MC2LaCone(0)
  ,fh2MC2ALaCone(0)  
  /*,fh2MCEtaVsPtHijingLa(0)
   ,fh2MCEtaVsPtInjectLa(0)
   ,fh2MCEtaVsPtHijingALa(0)
   ,fh2MCEtaVsPtInjectALa(0)
   ,fhnrecMCHijingLaIncl(0)
   ,fhnrecMCHijingLaCone(0)
   ,fhnrecMCHijingALaIncl(0)
   ,fhnrecMCHijingALaCone(0)
   ,fhnrecMCInjectLaIncl(0)
   ,fhnrecMCInjectLaCone(0)
   ,fhnrecMCInjectALaIncl(0)
   ,fhnrecMCInjectALaCone(0)*/
   ,fhnMCrecK0Cone(0)   
   ,fhnMCrecLaCone(0)   
   ,fhnMCrecALaCone(0)
  /*,fhnMCrecK0ConeSmear(0) 
   ,fhnMCrecLaConeSmear(0)   
   ,fhnMCrecALaConeSmear(0)
   ,fhnK0sSecContinCone(0)
   ,fhnLaSecContinCone(0)
   ,fhnALaSecContinCone(0)*/
   ,fhnK0sIncl(0)
   ,fhnK0sCone(0)
   ,fhnK0sEmbCone(0)
   ,fhnK0sEmbConeRef(0)
   ,fhnK0sEmbConeStandard(0)
   ,fhnLaIncl(0)
   ,fhnLaCone(0)
   ,fhnLaEmbCone(0)
   ,fhnLaEmbConeRef(0)
   ,fhnLaEmbConeStandard(0)
   ,fhnALaIncl(0)
   ,fhnALaCone(0)
   ,fhnALaEmbCone(0)
   ,fhnALaEmbConeRef(0)
   ,fhnALaEmbConeStandard(0)
   ,fh2MCEmbK0sJetPt(0)
   ,fh2MCEmbLaJetPt(0)
   ,fh2MCEmbALaJetPt(0)
   ,fhnK0sPC(0)
   ,fhnK0sEmbPC(0)
   ,fhnLaPC(0)
   ,fhnLaEmbPC(0)
   ,fhnALaPC(0)
   ,fhnALaEmbPC(0)
   ,fhnK0sMCC(0)
   ,fhnLaMCC(0)
   ,fhnALaMCC(0)
   ,fhnK0sRC(0)
   ,fhnLaRC(0)
   ,fhnALaRC(0)
   ,fhnK0sRCBias(0)
   ,fhnLaRCBias(0)
   ,fhnALaRCBias(0)
   ,fhnK0sOC(0)
   ,fhnLaOC(0)
   ,fhnALaOC(0)
   ,fh1AreaExcluded(0)
   ,fh1MedianEta(0)
   ,fh1JetPtMedian(0)
   ,fh1MCMultiplicityPrimary(0)
   ,fh1MCMultiplicityTracks(0)
   ,fhnFeedDownLa(0)
   ,fhnFeedDownALa(0)
   ,fhnFeedDownLaCone(0)
   ,fhnFeedDownALaCone(0)
   ,fh2FeedDownXiLa(0)
   ,fh2FeedDownXiALa(0)
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
  //,fh1MCRapK0s(0) 
  //,fh1MCRapLambda(0)
  //,fh1MCRapAntiLambda(0)
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

  ,fRandom(0)
  ,fMatchMode(0)
  ,fIsNJEventEmb(0)
  ,fAnalysisMC(0)
  ,fDeltaVertexZ(0)
  ,fCutjetEta(0)
  ,fCuttrackNegNcls(0)
  ,fCuttrackPosNcls(0)
  ,fCutPostrackRap(0)
  ,fCutNegtrackRap(0)
  ,fCutRap(0)
  ,fCutPostrackEta(0)
  ,fCutNegtrackEta(0)
  ,fCutEta(0)
  ,fCutK0cosPointAngle(0)
  ,fCutLacosPointAngle(0)
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
  ,fCutFractionPtEmbedded(0)
  ,fCutDeltaREmbedded(0)
    //,fBranchEmbeddedJets("") 
  ,fK0Type(0)  
  ,fFilterMaskK0(0)
  ,jettracklist(0)
  ,jetConeK0list(0)
  ,jetConeLalist(0)
  ,jetConeALalist(0)
  ,jetConeK0Emblist(0)
  ,jetConeLaEmblist(0)
  ,jetConeALaEmblist(0)
  ,jetConeK0EmbStlist(0)
  ,jetConeLaEmbStlist(0)
  ,jetConeALaEmbStlist(0)

  ,jetConeK0EmbMClist(0)
  ,jetConeLaEmbMClist(0)
  ,jetConeALaEmbMClist(0)

  ,jetPerpConeK0list(0)
  ,jetPerpRecCutslist(0)
  ,jetPerpConeK0Emblist(0)
  ,jetPerpConeLalist(0)
  ,jetPerpConeLaEmblist(0)
  ,jetPerpConeALalist(0)
  ,jetPerpConeALaEmblist(0)
  ,jetMedianConeK0list(0)
  ,jetMedianConeLalist(0)
  ,jetMedianConeALalist(0)
  ,jetMedianRecCutslist(0)
  ,fListK0sRC(0)
  ,fListLaRC(0)
  ,fListALaRC(0)
  ,fTracksRecCutsRC(0)
  ,fTracksRecBckgCuts(0)
  ,fTracksPerpCone(0)
  ,fListK0s(0)
  ,fListK0sMC(0)
  ,fListK0sStandard(0)
  ,fPIDResponse(0)
  ,fV0QAK0(0)
  ,fFFHistosRecCutsK0Evt(0)      
    //,fFFHistosIMK0AllEvt(0)        
    //,fFFHistosIMK0Jet(0)           
    //,fFFHistosIMK0Cone(0)
  ,fLaType(0)  
  ,fFilterMaskLa(0)
  ,fListLa(0)
  ,fListLaMC(0)
  ,fListLaStandard(0)
    //,fFFHistosIMLaAllEvt(0)        
    //,fFFHistosIMLaJet(0)           
    //,fFFHistosIMLaCone(0)
  ,fALaType(0)  
  ,fFilterMaskALa(0)
  ,fListALa(0)
  ,fListALaMC(0)
  ,fListALaStandard(0)
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
  ,fUseExtraTracks(0)
  ,fUseExtraJetPt(0)
  ,fUseEmbeddedJetPt(0)

  ,fUseStandard(0) 
    //,fFFHistosIMALaAllEvt(0)        
    //,fFFHistosIMALaJet(0)           
    // ,fFFHistosIMALaCone(0)
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
  ,fh1nEmbeddedJets(0)
  ,fh1nGenJets(0)
  ,fh1IndexEmbedded(0)
  ,fh1IndexEmbeddedMC(0)
  ,fh1PtEmbBeforeMatch(0)
  ,fh1PtEmbExtraOnly(0)
  ,fh1PtEmbReject(0)
  ,fh2PtEtaEmbReject(0)
  ,fh1PtEmbAfterMatch(0)
  ,fh1FractionPtEmbedded(0)
  ,fh1DeltaREmbedded(0)
  ,fh1FractionPtEmbeddedMC(0)
  ,fh2FractionPtVsEmbeddedJetPtMC(0)
  ,fh1DeltaREmbeddedMC(0)
  ,fh1JetPtEmbGenAfterMatch(0)
  ,fh2TracksPerpCone(0)
  ,fh1PerpCone(0)
    //  ,fh1V0JetPt(0)
  ,fh1V0PtCandidate(0)
  ,fh1IMK0Cone(0)
  ,fh1IMLaCone(0)
  ,fh1IMALaCone(0)
  ,fh1IMK0EmbCone(0)
  ,fh1IMLaEmbCone(0)
  ,fh1IMALaEmbCone(0)
  ,fh2FFJetTrackEta(0)  
    //  ,fh1trackPosNCls(0)           
    // ,fh1trackNegNCls(0) 
  ,fh1trackPosRap(0)            
  ,fh1trackNegRap(0)          
    //,fh1V0Rap(0)          
  ,fh1trackPosEta(0)            
  ,fh1trackNegEta(0)          
  ,fh1V0Eta(0)  
    // ,fh1V0totMom(0)            
  ,fh1CosPointAngle(0)        
  ,fh1DecayLengthV0(0) 
  ,fh2ProperLifetimeK0sVsPtBeforeCut(0)  
  ,fh2ProperLifetimeK0sVsPtAfterCut(0)            
  ,fh1V0Radius(0)       
  ,fh1DcaV0Daughters(0)        
  ,fh1DcaPosToPrimVertex(0)   
  ,fh1DcaNegToPrimVertex(0)    
  ,fh2ArmenterosBeforeCuts(0)
  ,fh2ArmenterosAfterCuts(0)
  ,fh2BBLaPos(0)
  ,fh2BBLaNeg(0)
  ,fh1PosDaughterCharge(0)
  ,fh1NegDaughterCharge(0)
  ,fh1PtMCK0s(0)
  ,fh1PtMCLa(0)
  ,fh1PtMCALa(0)
  ,fh1EtaK0s(0)
  ,fh1EtaLa(0)
  ,fh1EtaALa(0)
  ,fh1RC(0)
  ,fh1RCBiasK0(0)
  ,fh1RCBiasLa(0)
  ,fh1RCBiasALa(0)
  ,fh1MCC(0)
  ,fh1OC(0)
  ,fh1NJ(0)
  ,fh1NJEmbEvt(0)
  ,fh1BckgJets(0)
  ,fh1BckgJetsPtBias(0)
  ,fhnInvMassEtaTrackPtK0s(0)
  ,fhnInvMassEtaTrackPtLa(0)
  ,fhnInvMassEtaTrackPtALa(0)
  ,fh1TrackMultCone(0)
  ,fh2TrackMultCone(0)
  ,fhnNJK0(0)
  ,fhnNJLa(0)
  ,fhnNJALa(0)
  //,fh2ChTracksNJ(0)
  //,fh2ChTracksRC(0)
    // ,fh2ChTracksOC(0)
  //,fh2ChTracksMCC(0)
  //,fh2ChTracksPC(0)
   ,fhnMCgenK0Cone(0)
   ,fhnMCgenLaCone(0)
    //,fh2MCgenALaCone(0)
  //,fh2MCEtagenK0Cone(0)
  //,fh2MCEtagenLaCone(0)
    //,fh2MCEtagenALaCone(0)
  /* ,fh2CorrHijingLaProton(0)
  ,fh2CorrInjectLaProton(0)
  ,fh2CorrHijingALaAProton(0)
  ,fh2CorrInjectALaAProton(0)*/
  ,fh1IMK0ConeSmear(0)
  ,fh1IMLaConeSmear(0)
  ,fh1IMALaConeSmear(0) 
  ,fh2MC2K0Cone(0)
  ,fh2MC2LaCone(0)
  ,fh2MC2ALaCone(0)
  /* ,fh2MCEtaVsPtHijingLa(0)
  ,fh2MCEtaVsPtInjectLa(0)
  ,fh2MCEtaVsPtHijingALa(0)
  ,fh2MCEtaVsPtInjectALa(0) 
  ,fhnrecMCHijingLaIncl(0)
  ,fhnrecMCHijingLaCone(0)
  ,fhnrecMCHijingALaIncl(0)
  ,fhnrecMCHijingALaCone(0)
  ,fhnrecMCInjectLaIncl(0)
  ,fhnrecMCInjectLaCone(0)
  ,fhnrecMCInjectALaIncl(0)
  ,fhnrecMCInjectALaCone(0)*/
  ,fhnMCrecK0Cone(0)
  ,fhnMCrecLaCone(0)
  ,fhnMCrecALaCone(0) 
  /*,fhnMCrecK0ConeSmear(0) 
  ,fhnMCrecLaConeSmear(0)   
  ,fhnMCrecALaConeSmear(0)
 ,fhnK0sSecContinCone(0)
  ,fhnLaSecContinCone(0)
  ,fhnALaSecContinCone(0)*/
  ,fhnK0sIncl(0)
  ,fhnK0sCone(0)
  ,fhnK0sEmbCone(0)
  ,fhnK0sEmbConeRef(0)
  ,fhnK0sEmbConeStandard(0)
  ,fhnLaIncl(0)
  ,fhnLaCone(0)
  ,fhnLaEmbCone(0)
  ,fhnLaEmbConeRef(0)
  ,fhnLaEmbConeStandard(0)
  ,fhnALaIncl(0)
  ,fhnALaCone(0)
  ,fhnALaEmbCone(0)
  ,fhnALaEmbConeRef(0)
  ,fhnALaEmbConeStandard(0)
  ,fh2MCEmbK0sJetPt(0)
  ,fh2MCEmbLaJetPt(0)
  ,fh2MCEmbALaJetPt(0)
  ,fhnK0sPC(0)
  ,fhnK0sEmbPC(0)
  ,fhnLaPC(0)
  ,fhnLaEmbPC(0)
  ,fhnALaPC(0)
  ,fhnALaEmbPC(0)
  ,fhnK0sMCC(0)
  ,fhnLaMCC(0)
  ,fhnALaMCC(0)
  ,fhnK0sRC(0)
  ,fhnLaRC(0)
  ,fhnALaRC(0)
  ,fhnK0sRCBias(0)
  ,fhnLaRCBias(0)
  ,fhnALaRCBias(0)
  ,fhnK0sOC(0)
  ,fhnLaOC(0)
  ,fhnALaOC(0)
  ,fh1AreaExcluded(0)
  ,fh1MedianEta(0)
  ,fh1JetPtMedian(0)
  ,fh1MCMultiplicityPrimary(0)
  ,fh1MCMultiplicityTracks(0)
  ,fhnFeedDownLa(0)
  ,fhnFeedDownALa(0)
  ,fhnFeedDownLaCone(0)
  ,fhnFeedDownALaCone(0)
  ,fh2FeedDownXiLa(0)
  ,fh2FeedDownXiALa(0)
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
    //,fh1MCRapK0s(0) 
    //,fh1MCRapLambda(0)
    //,fh1MCRapAntiLambda(0)
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
  
  ,fRandom(copy.fRandom)
  ,fMatchMode(copy.fMatchMode)
  ,fIsNJEventEmb(copy.fIsNJEventEmb)
  ,fAnalysisMC(copy.fAnalysisMC)
  ,fDeltaVertexZ(copy.fDeltaVertexZ)
  ,fCutjetEta(copy.fCutjetEta)
  ,fCuttrackNegNcls(copy.fCuttrackNegNcls)
  ,fCuttrackPosNcls(copy.fCuttrackPosNcls)
  ,fCutPostrackRap(copy.fCutPostrackRap)
  ,fCutNegtrackRap(copy.fCutNegtrackRap)
  ,fCutRap(copy.fCutRap)
  ,fCutPostrackEta(copy.fCutPostrackEta)
  ,fCutNegtrackEta(copy.fCutNegtrackEta)
  ,fCutEta(copy.fCutEta)
  ,fCutK0cosPointAngle(copy.fCutK0cosPointAngle)
  ,fCutLacosPointAngle(copy.fCutLacosPointAngle)
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
  ,fCutFractionPtEmbedded(copy.fCutFractionPtEmbedded)
  ,fCutDeltaREmbedded(copy.fCutDeltaREmbedded)
    //,fBranchEmbeddedJets(copy.fBranchEmbeddedJets) 
  ,fK0Type(copy.fK0Type)              
  ,fFilterMaskK0(copy.fFilterMaskK0)
  ,jettracklist(copy.jettracklist)
  ,jetConeK0list(copy.jetConeK0list)
  ,jetConeLalist(copy.jetConeLalist)
  ,jetConeALalist(copy.jetConeALalist)
  ,jetConeK0Emblist(copy.jetConeK0Emblist)
  ,jetConeLaEmblist(copy.jetConeLaEmblist)
  ,jetConeALaEmblist(copy.jetConeALaEmblist)
  ,jetConeK0EmbStlist(copy.jetConeK0EmbStlist)
  ,jetConeLaEmbStlist(copy.jetConeLaEmbStlist)
  ,jetConeALaEmbStlist(copy.jetConeALaEmbStlist)

  ,jetConeK0EmbMClist(copy.jetConeK0EmbMClist)
  ,jetConeLaEmbMClist(copy.jetConeLaEmbMClist)
  ,jetConeALaEmbMClist(copy.jetConeALaEmbMClist)

  ,jetPerpConeK0list(copy.jetPerpConeK0list)
  ,jetPerpRecCutslist(copy.jetPerpRecCutslist)
  ,jetPerpConeK0Emblist(copy.jetPerpConeK0Emblist)
  ,jetPerpConeLalist(copy.jetPerpConeLalist)
  ,jetPerpConeLaEmblist(copy.jetPerpConeLaEmblist)
  ,jetPerpConeALalist(copy.jetPerpConeALalist)
  ,jetPerpConeALaEmblist(copy.jetPerpConeALaEmblist)
  ,jetMedianConeK0list(copy.jetMedianConeK0list)
  ,jetMedianConeLalist(copy.jetMedianConeLalist)
  ,jetMedianConeALalist(copy.jetMedianConeALalist)
  ,jetMedianRecCutslist(copy.jetMedianRecCutslist)
  ,fListK0sRC(copy.fListK0sRC)
  ,fListLaRC(copy.fListLaRC)
  ,fListALaRC(copy.fListALaRC)
  ,fTracksRecCutsRC(copy.fTracksRecCutsRC) 
  ,fTracksRecBckgCuts(copy.fTracksRecBckgCuts) 
  ,fTracksPerpCone(copy.fTracksPerpCone)
  ,fListK0s(copy.fListK0s)
  ,fListK0sMC(copy.fListK0sMC)
  ,fListK0sStandard(copy.fListK0sStandard)
  ,fPIDResponse(copy.fPIDResponse)
  ,fV0QAK0(copy.fV0QAK0)
  ,fFFHistosRecCutsK0Evt(copy.fFFHistosRecCutsK0Evt)      
    //,fFFHistosIMK0AllEvt(copy.fFFHistosIMK0AllEvt)        
    //,fFFHistosIMK0Jet(copy.fFFHistosIMK0Jet)           
    //,fFFHistosIMK0Cone(copy.fFFHistosIMK0Cone)          
  ,fLaType(copy.fLaType)                  
  ,fFilterMaskLa(copy.fFilterMaskLa)
  ,fListLa(copy.fListLa)
  ,fListLaMC(copy.fListLaMC)
  ,fListLaStandard(copy.fListLaStandard)
    //,fFFHistosIMLaAllEvt(copy.fFFHistosIMLaAllEvt)        
    //,fFFHistosIMLaJet(copy.fFFHistosIMLaJet)           
    //,fFFHistosIMLaCone(copy.fFFHistosIMLaCone)          
  ,fALaType(copy.fALaType)                 
  ,fFilterMaskALa(copy.fFilterMaskALa)
  ,fListALa(copy.fListALa)
  ,fListALaMC(copy.fListALaMC)
  ,fListALaStandard(copy.fListALaStandard)
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
  ,fUseExtraTracks(copy.fUseExtraTracks)
  ,fUseExtraJetPt(copy.fUseExtraJetPt)
  ,fUseEmbeddedJetPt(copy.fUseEmbeddedJetPt)

  ,fUseStandard(copy.fUseStandard)
    //,fFFHistosIMALaAllEvt(copy.fFFHistosIMALaAllEvt)        
    //,fFFHistosIMALaJet(copy.fFFHistosIMALaJet)           
    //,fFFHistosIMALaCone(copy.fFFHistosIMALaCone)          
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
  ,fh1nEmbeddedJets(copy.fh1nEmbeddedJets)
  ,fh1nGenJets(copy.fh1nGenJets)
  ,fh1IndexEmbedded(copy.fh1IndexEmbedded)
  ,fh1IndexEmbeddedMC(copy.fh1IndexEmbeddedMC)
  ,fh1PtEmbBeforeMatch(copy.fh1PtEmbBeforeMatch)
  ,fh1PtEmbExtraOnly(copy.fh1PtEmbExtraOnly)
  ,fh1PtEmbReject(copy.fh1PtEmbReject)
  ,fh2PtEtaEmbReject(copy.fh2PtEtaEmbReject)
  ,fh1PtEmbAfterMatch(copy.fh1PtEmbAfterMatch)
  ,fh1FractionPtEmbedded(copy.fh1FractionPtEmbedded)
  ,fh1DeltaREmbedded(copy.fh1DeltaREmbedded)
  ,fh1FractionPtEmbeddedMC(copy.fh1FractionPtEmbeddedMC)
  ,fh2FractionPtVsEmbeddedJetPtMC(copy.fh2FractionPtVsEmbeddedJetPtMC)
  ,fh1DeltaREmbeddedMC(copy.fh1DeltaREmbeddedMC)
  ,fh1JetPtEmbGenAfterMatch(copy.fh1JetPtEmbGenAfterMatch)
  ,fh2TracksPerpCone(copy.fh2TracksPerpCone)
  ,fh1PerpCone(copy.fh1PerpCone)
    //,fh1V0JetPt(copy.fh1V0JetPt)
  ,fh1V0PtCandidate(copy.fh1V0PtCandidate)
  ,fh1IMK0Cone(copy.fh1IMK0Cone)
  ,fh1IMLaCone(copy.fh1IMLaCone)
  ,fh1IMALaCone(copy.fh1IMALaCone)
  ,fh1IMK0EmbCone(copy.fh1IMK0EmbCone)
  ,fh1IMLaEmbCone(copy.fh1IMLaEmbCone)
  ,fh1IMALaEmbCone(copy.fh1IMALaEmbCone)
  ,fh2FFJetTrackEta(copy.fh2FFJetTrackEta) 
    //,fh1trackPosNCls(copy.fh1trackPosNCls)           
    //,fh1trackNegNCls(copy.fh1trackNegNCls)
  ,fh1trackPosRap(copy.fh1trackPosRap)            
  ,fh1trackNegRap(copy.fh1trackNegRap)          
    //,fh1V0Rap(copy.fh1V0Rap)         
  ,fh1trackPosEta(copy.fh1trackPosEta)            
  ,fh1trackNegEta(copy.fh1trackNegEta)          
  ,fh1V0Eta(copy.fh1V0Eta)   
    //,fh1V0totMom(copy.fh1V0totMom)           
  ,fh1CosPointAngle(copy.fh1CosPointAngle)           
  ,fh1DecayLengthV0(copy.fh1DecayLengthV0)  
  ,fh2ProperLifetimeK0sVsPtBeforeCut(copy.fh2ProperLifetimeK0sVsPtBeforeCut)  
  ,fh2ProperLifetimeK0sVsPtAfterCut(copy.fh2ProperLifetimeK0sVsPtAfterCut)    
  ,fh1V0Radius(copy.fh1V0Radius)          
  ,fh1DcaV0Daughters(copy.fh1DcaV0Daughters)        
  ,fh1DcaPosToPrimVertex(copy.fh1DcaPosToPrimVertex)   
  ,fh1DcaNegToPrimVertex(copy.fh1DcaNegToPrimVertex)    
  ,fh2ArmenterosBeforeCuts(copy.fh2ArmenterosBeforeCuts)
  ,fh2ArmenterosAfterCuts(copy.fh2ArmenterosAfterCuts)
  ,fh2BBLaPos(copy.fh2BBLaPos)
  ,fh2BBLaNeg(copy.fh2BBLaPos)
  ,fh1PosDaughterCharge(copy.fh1PosDaughterCharge)
  ,fh1NegDaughterCharge(copy.fh1NegDaughterCharge)
  ,fh1PtMCK0s(copy.fh1PtMCK0s)
  ,fh1PtMCLa(copy.fh1PtMCLa)
  ,fh1PtMCALa(copy.fh1PtMCALa)
  ,fh1EtaK0s(copy.fh1EtaK0s)
  ,fh1EtaLa(copy.fh1EtaLa)
  ,fh1EtaALa(copy.fh1EtaALa)
  ,fh1RC(copy.fh1RC)
  ,fh1RCBiasK0(copy.fh1RCBiasK0)
  ,fh1RCBiasLa(copy.fh1RCBiasLa)
  ,fh1RCBiasALa(copy.fh1RCBiasALa)
  ,fh1MCC(copy.fh1MCC)
  ,fh1OC(copy.fh1OC)
  ,fh1NJ(copy.fh1NJ)
  ,fh1NJEmbEvt(copy.fh1NJEmbEvt)
  ,fh1BckgJets(copy.fh1BckgJets)
  ,fh1BckgJetsPtBias(copy.fh1BckgJetsPtBias)
  ,fhnInvMassEtaTrackPtK0s(copy.fhnInvMassEtaTrackPtK0s)
  ,fhnInvMassEtaTrackPtLa(copy.fhnInvMassEtaTrackPtLa)
  ,fhnInvMassEtaTrackPtALa(copy.fhnInvMassEtaTrackPtALa)
  ,fh1TrackMultCone(copy.fh1TrackMultCone)
  ,fh2TrackMultCone(copy.fh2TrackMultCone)
  ,fhnNJK0(copy.fhnNJK0)
  ,fhnNJLa(copy.fhnNJLa)
  ,fhnNJALa(copy.fhnNJALa)
  //,fh2ChTracksNJ(copy.fh2ChTracksNJ)
  //,fh2ChTracksRC(copy.fh2ChTracksRC)
    // ,fh2ChTracksOC(copy.fh2ChTracksOC)
  //,fh2ChTracksMCC(copy.fh2ChTracksMCC)
  //,fh2ChTracksPC(copy.fh2ChTracksPC)
   ,fhnMCgenK0Cone(copy.fhnMCgenK0Cone)
   ,fhnMCgenLaCone(copy.fhnMCgenLaCone)
    //,fh2MCgenALaCone(copy.fh2MCgenALaCone)
  //,fh2MCEtagenK0Cone(copy.fh2MCEtagenK0Cone)
  // ,fh2MCEtagenLaCone(copy.fh2MCEtagenLaCone)
    //,fh2MCEtagenALaCone(copy.fh2MCEtagenALaCone)
  /* ,fh2CorrHijingLaProton(copy.fh2CorrHijingLaProton)
  ,fh2CorrInjectLaProton(copy.fh2CorrInjectLaProton)
  ,fh2CorrHijingALaAProton(copy.fh2CorrHijingALaAProton)
  ,fh2CorrInjectALaAProton(copy.fh2CorrInjectALaAProton)*/
  ,fh1IMK0ConeSmear(copy.fh1IMK0ConeSmear)
  ,fh1IMLaConeSmear(copy.fh1IMLaConeSmear)
  ,fh1IMALaConeSmear(copy.fh1IMALaConeSmear)
  ,fh2MC2K0Cone(copy.fh2MC2K0Cone)
  ,fh2MC2LaCone(copy.fh2MC2LaCone)
  ,fh2MC2ALaCone(copy.fh2MC2ALaCone)
  /*
  ,fh2MCEtaVsPtHijingLa(copy.fh2MCEtaVsPtHijingLa)
  ,fh2MCEtaVsPtInjectLa(copy.fh2MCEtaVsPtInjectLa)
  ,fh2MCEtaVsPtHijingALa(copy.fh2MCEtaVsPtHijingALa)
  ,fh2MCEtaVsPtInjectALa(copy.fh2MCEtaVsPtInjectALa)
  ,fhnrecMCHijingLaIncl(copy.fhnrecMCHijingLaIncl)
  ,fhnrecMCHijingLaCone(copy.fhnrecMCHijingLaCone)
  ,fhnrecMCHijingALaIncl(copy.fhnrecMCHijingALaIncl)
  ,fhnrecMCHijingALaCone(copy.fhnrecMCHijingALaCone)
  ,fhnrecMCInjectLaIncl(copy.fhnrecMCInjectLaIncl)
  ,fhnrecMCInjectLaCone(copy.fhnrecMCInjectLaCone)
  ,fhnrecMCInjectALaIncl(copy.fhnrecMCInjectALaIncl)
  ,fhnrecMCInjectALaCone(copy.fhnrecMCInjectALaCone)   */
  ,fhnMCrecK0Cone(copy.fhnMCrecK0Cone)
  ,fhnMCrecLaCone(copy.fhnMCrecLaCone)
  ,fhnMCrecALaCone(copy.fhnMCrecALaCone) 
  /*,fhnMCrecK0ConeSmear(copy.fhnMCrecK0ConeSmear)
  ,fhnMCrecLaConeSmear(copy.fhnMCrecLaConeSmear)
  ,fhnMCrecALaConeSmear(copy.fhnMCrecALaConeSmear)
  ,fhnK0sSecContinCone(copy.fhnK0sSecContinCone)
  ,fhnLaSecContinCone(copy.fhnLaSecContinCone)
  ,fhnALaSecContinCone(copy.fhnALaSecContinCone)*/
  ,fhnK0sIncl(copy.fhnK0sIncl)
  ,fhnK0sCone(copy.fhnK0sCone)
  ,fhnK0sEmbCone(copy.fhnK0sEmbCone)
  ,fhnK0sEmbConeRef(copy.fhnK0sEmbConeRef)
  ,fhnK0sEmbConeStandard(copy.fhnK0sEmbConeStandard)
  ,fhnLaIncl(copy.fhnLaIncl)
  ,fhnLaCone(copy.fhnLaCone)
  ,fhnLaEmbCone(copy.fhnLaEmbCone)
  ,fhnLaEmbConeRef(copy.fhnLaEmbConeRef)
  ,fhnLaEmbConeStandard(copy.fhnLaEmbConeStandard)
  ,fhnALaIncl(copy.fhnALaIncl)
  ,fhnALaCone(copy.fhnALaCone)
  ,fhnALaEmbCone(copy.fhnALaEmbCone)
  ,fhnALaEmbConeRef(copy.fhnALaEmbConeRef)
  ,fhnALaEmbConeStandard(copy.fhnALaEmbConeStandard)

  ,fh2MCEmbK0sJetPt(copy.fh2MCEmbK0sJetPt)
  ,fh2MCEmbLaJetPt(copy.fh2MCEmbLaJetPt)
  ,fh2MCEmbALaJetPt(copy.fh2MCEmbALaJetPt)

  ,fhnK0sPC(copy.fhnK0sPC)
  ,fhnK0sEmbPC(copy.fhnK0sEmbPC)
  ,fhnLaPC(copy.fhnLaPC)
  ,fhnLaEmbPC(copy.fhnLaEmbPC)
  ,fhnALaPC(copy.fhnALaPC)
  ,fhnALaEmbPC(copy.fhnALaEmbPC)
  ,fhnK0sMCC(copy.fhnK0sMCC)
  ,fhnLaMCC(copy.fhnLaMCC)
  ,fhnALaMCC(copy.fhnALaMCC)
  ,fhnK0sRC(copy.fhnK0sRC)
  ,fhnLaRC(copy.fhnLaRC)
  ,fhnALaRC(copy.fhnALaRC)
  ,fhnK0sRCBias(copy.fhnK0sRCBias)
  ,fhnLaRCBias(copy.fhnLaRCBias)
  ,fhnALaRCBias(copy.fhnALaRCBias)
  ,fhnK0sOC(copy.fhnK0sOC)
  ,fhnLaOC(copy.fhnLaOC)
  ,fhnALaOC(copy.fhnALaOC)
  ,fh1AreaExcluded(copy.fh1AreaExcluded)
  ,fh1MedianEta(copy.fh1MedianEta)
  ,fh1JetPtMedian(copy.fh1JetPtMedian)
  ,fh1MCMultiplicityPrimary(copy.fh1MCMultiplicityPrimary)
  ,fh1MCMultiplicityTracks(copy.fh1MCMultiplicityTracks)
  ,fhnFeedDownLa(copy.fhnFeedDownLa)
  ,fhnFeedDownALa(copy.fhnFeedDownALa)
  ,fhnFeedDownLaCone(copy.fhnFeedDownLaCone)
  ,fhnFeedDownALaCone(copy.fhnFeedDownALaCone)
  ,fh2FeedDownXiLa(copy.fh2FeedDownXiLa)
  ,fh2FeedDownXiALa(copy.fh2FeedDownXiALa)  
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
    //,fh1MCRapK0s(copy.fh1MCRapK0s) 
    //,fh1MCRapLambda(copy.fh1MCRapLambda)
    //,fh1MCRapAntiLambda(copy.fh1MCRapAntiLambda)
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

    fRandom                         = o.fRandom;
    fMatchMode                      = o.fMatchMode;
    fIsNJEventEmb                   = o.fIsNJEventEmb; 
    fAnalysisMC                     = o.fAnalysisMC;
    fDeltaVertexZ                   = o.fDeltaVertexZ;
    fCutjetEta                      = o.fCutjetEta;
    fCuttrackNegNcls                = o.fCuttrackNegNcls;
    fCuttrackPosNcls                = o.fCuttrackPosNcls;
    fCutPostrackRap                 = o.fCutPostrackRap;
    fCutNegtrackRap                 = o.fCutNegtrackRap;  
    fCutRap                         = o.fCutRap;
    fCutPostrackEta                 = o.fCutPostrackEta;
    fCutNegtrackEta                 = o.fCutNegtrackEta;  
    fCutEta                         = o.fCutEta;
    fCutK0cosPointAngle             = o.fCutK0cosPointAngle;
    fCutLacosPointAngle             = o.fCutLacosPointAngle;
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
    fCutFractionPtEmbedded          = o.fCutFractionPtEmbedded;
    fCutDeltaREmbedded              = o.fCutDeltaREmbedded;
    //fBranchEmbeddedJets             = o.fBranchEmbeddedJets; 
    fK0Type                         = o.fK0Type;
    fFilterMaskK0                   = o.fFilterMaskK0;
    jettracklist                    = o.jettracklist;
    jetConeK0list                   = o.jetConeK0list;
    jetConeLalist                   = o.jetConeLalist;
    jetConeALalist                  = o.jetConeALalist;
    jetConeK0Emblist                = o.jetConeK0Emblist;
    jetConeLaEmblist                = o.jetConeLaEmblist;
    jetConeALaEmblist               = o.jetConeALaEmblist;
    jetConeK0EmbStlist              = o.jetConeK0EmbStlist;
    jetConeLaEmbStlist              = o.jetConeLaEmbStlist;
    jetConeALaEmbStlist             = o.jetConeALaEmbStlist;

    jetConeK0EmbMClist              = o.jetConeK0EmbMClist;
    jetConeLaEmbMClist              = o.jetConeLaEmbMClist;
    jetConeALaEmbMClist             = o.jetConeALaEmbMClist;

    jetPerpConeK0list               = o.jetPerpConeK0list;
    jetPerpRecCutslist              = o.jetPerpRecCutslist;
    jetPerpConeK0Emblist            = o.jetPerpConeK0Emblist;
    jetPerpConeLalist               = o.jetPerpConeLalist;
    jetPerpConeLaEmblist            = o.jetPerpConeLaEmblist;
    jetPerpConeALalist              = o.jetPerpConeALalist;
    jetPerpConeALaEmblist           = o.jetPerpConeALaEmblist;
    jetMedianConeK0list             = o.jetMedianConeK0list;
    jetMedianConeLalist             = o.jetMedianConeLalist;
    jetMedianConeALalist            = o.jetMedianConeALalist;
    jetMedianRecCutslist            = o.jetMedianRecCutslist;
    fListK0sRC                      = o.fListK0sRC;
    fListLaRC                       = o.fListLaRC;
    fListALaRC                      = o.fListALaRC;
    fTracksRecCutsRC                = o.fTracksRecCutsRC;
    fTracksRecBckgCuts              = o.fTracksRecBckgCuts;
    fTracksPerpCone                 = o.fTracksPerpCone;
    fListK0s                        = o.fListK0s;
    fListK0sMC                      = o.fListK0sMC;
    fListK0sStandard                = o.fListK0sStandard;
    fPIDResponse                    = o.fPIDResponse;
    fV0QAK0                         = o.fV0QAK0;
    fFFHistosRecCutsK0Evt           = o.fFFHistosRecCutsK0Evt;      
    //fFFHistosIMK0AllEvt             = o.fFFHistosIMK0AllEvt;        
    //fFFHistosIMK0Jet                = o.fFFHistosIMK0Jet;           
    //fFFHistosIMK0Cone               = o.fFFHistosIMK0Cone;          
    fLaType                         = o.fLaType;
    fFilterMaskLa                   = o.fFilterMaskLa;
    fListLa                         = o.fListLa;
    fListLaMC                       = o.fListLaMC;
    fListLaStandard                 = o.fListLaStandard;
    //fFFHistosIMLaAllEvt             = o.fFFHistosIMLaAllEvt;        
    //fFFHistosIMLaJet                = o.fFFHistosIMLaJet;           
    //fFFHistosIMLaCone               = o.fFFHistosIMLaCone;          
    fALaType                        = o.fALaType;
    fFilterMaskALa                  = o.fFilterMaskALa;
    fListALa                        = o.fListALa;
    fListALaMC                      = o.fListALaMC;
    fListALaStandard                = o.fListALaStandard;
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
    fUseExtraTracks                 = o.fUseExtraTracks;
    fUseExtraJetPt                  = o.fUseExtraJetPt;
    fUseEmbeddedJetPt               = o.fUseEmbeddedJetPt;

    fUseStandard                    = o.fUseStandard;
    // fFFHistosIMALaAllEvt            = o.fFFHistosIMALaAllEvt;        
    // fFFHistosIMALaJet               = o.fFFHistosIMALaJet;           
    // fFFHistosIMALaCone              = o.fFFHistosIMALaCone;          
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
    fh1nEmbeddedJets                = o.fh1nEmbeddedJets;
    fh1nGenJets                     = o.fh1nGenJets; 
    fh1IndexEmbedded                = o.fh1IndexEmbedded;
    fh1IndexEmbeddedMC              = o.fh1IndexEmbeddedMC;
    fh1PtEmbBeforeMatch             = o.fh1PtEmbBeforeMatch;
    fh1PtEmbExtraOnly               = o.fh1PtEmbExtraOnly;
    fh1PtEmbReject                  = o.fh1PtEmbReject;
    fh2PtEtaEmbReject               = o.fh2PtEtaEmbReject;
    fh1PtEmbAfterMatch              = o.fh1PtEmbAfterMatch;
    fh1FractionPtEmbedded           = o.fh1FractionPtEmbedded;
    fh1DeltaREmbedded               = o.fh1DeltaREmbedded;
    fh1FractionPtEmbeddedMC        = o.fh1FractionPtEmbeddedMC;
    fh2FractionPtVsEmbeddedJetPtMC = o.fh2FractionPtVsEmbeddedJetPtMC;
    fh1DeltaREmbeddedMC            = o.fh1DeltaREmbeddedMC;
    fh1JetPtEmbGenAfterMatch       = o.fh1JetPtEmbGenAfterMatch;
    fh2TracksPerpCone               = o.fh2TracksPerpCone;
    fh1PerpCone                     = o.fh1PerpCone;
    //fh1V0JetPt                     = o.fh1V0JetPt;
    fh1V0PtCandidate                = o.fh1V0PtCandidate;
    fh1IMK0Cone                     = o.fh1IMK0Cone;
    fh1IMLaCone                     = o.fh1IMLaCone;
    fh1IMALaCone                    = o.fh1IMALaCone;
    fh1IMK0EmbCone                  = o.fh1IMK0EmbCone;
    fh1IMLaEmbCone                  = o.fh1IMLaEmbCone;
    fh1IMALaEmbCone                 = o.fh1IMALaEmbCone;
    fh2FFJetTrackEta                = o.fh2FFJetTrackEta; 
    //fh1trackPosNCls                 = o.fh1trackPosNCls;           
    //fh1trackNegNCls                 = o.fh1trackNegNCls;    
    fh1trackPosRap                  = o.fh1trackPosRap;            
    fh1trackNegRap                  = o.fh1trackNegRap;        
    //fh1V0Rap                        = o.fh1V0Rap;        
    fh1trackPosEta                  = o.fh1trackPosEta;            
    fh1trackNegEta                  = o.fh1trackNegEta;        
    fh1V0Eta                        = o.fh1V0Eta;  
    // fh1V0totMom                     = o.fh1V0totMom;            
    fh1CosPointAngle                = o.fh1CosPointAngle;                      
    fh1DecayLengthV0                = o.fh1DecayLengthV0;  
    fh2ProperLifetimeK0sVsPtBeforeCut = o.fh2ProperLifetimeK0sVsPtBeforeCut;
    fh2ProperLifetimeK0sVsPtAfterCut= o.fh2ProperLifetimeK0sVsPtAfterCut; 
    fh1V0Radius                     = o.fh1V0Radius;         
    fh1DcaV0Daughters               = o.fh1DcaV0Daughters;        
    fh1DcaPosToPrimVertex           = o.fh1DcaPosToPrimVertex;   
    fh1DcaNegToPrimVertex           = o.fh1DcaNegToPrimVertex;    
    fh2ArmenterosBeforeCuts         = o.fh2ArmenterosBeforeCuts;
    fh2ArmenterosAfterCuts          = o.fh2ArmenterosAfterCuts;
    fh2BBLaPos                      = o.fh2BBLaPos;
    fh2BBLaNeg                      = o.fh2BBLaPos;
    fh1PosDaughterCharge            = o.fh1PosDaughterCharge;
    fh1NegDaughterCharge            = o.fh1NegDaughterCharge;
    fh1PtMCK0s                      = o.fh1PtMCK0s;
    fh1PtMCLa                       = o.fh1PtMCLa;
    fh1PtMCALa                      = o.fh1PtMCALa;
    fh1EtaK0s                       = o.fh1EtaK0s;
    fh1EtaLa                        = o.fh1EtaLa;
    fh1EtaALa                       = o.fh1EtaALa;
    fh1RC                           = o.fh1RC;
    fh1RCBiasK0                     = o.fh1RCBiasK0;
    fh1RCBiasLa                     = o.fh1RCBiasLa;
    fh1RCBiasALa                    = o.fh1RCBiasALa;
    fh1MCC                          = o.fh1MCC;
    fh1OC                           = o.fh1OC;
    fh1NJ                           = o.fh1NJ;
    fh1NJEmbEvt                     = o.fh1NJEmbEvt;
    fh1BckgJets                     = o.fh1BckgJets;
    fh1BckgJetsPtBias               = o.fh1BckgJetsPtBias;
    fhnInvMassEtaTrackPtK0s         = o.fhnInvMassEtaTrackPtK0s;
    fhnInvMassEtaTrackPtLa          = o.fhnInvMassEtaTrackPtLa;
    fhnInvMassEtaTrackPtALa         = o.fhnInvMassEtaTrackPtALa;
    fh1TrackMultCone                = o.fh1TrackMultCone;
    fh2TrackMultCone                = o.fh2TrackMultCone;
    fhnNJK0                         = o.fhnNJK0;
    fhnNJLa                         = o.fhnNJLa;
    fhnNJALa                        = o.fhnNJALa;
    //fh2ChTracksNJ                   = o.fh2ChTracksNJ;
    //fh2ChTracksRC                   = o.fh2ChTracksRC;
    //   fh2ChTracksOC                   = o.fh2ChTracksOC;
    //fh2ChTracksMCC                  = o.fh2ChTracksMCC;
    //fh2ChTracksPC                   = o.fh2ChTracksPC;
    fhnMCgenK0Cone                  = o.fhnMCgenK0Cone;
    fhnMCgenLaCone                  = o.fhnMCgenLaCone;
    //fh2MCgenALaCone                 = o.fh2MCgenALaCone; 
    //fh2MCEtagenK0Cone               = o.fh2MCEtagenK0Cone;
    //fh2MCEtagenLaCone               = o.fh2MCEtagenLaCone;
    //fh2MCEtagenALaCone              = o.fh2MCEtagenALaCone;
    fh1IMK0ConeSmear                = o.fh1IMK0ConeSmear;
    fh1IMLaConeSmear                = o.fh1IMLaConeSmear;
    fh1IMALaConeSmear               = o.fh1IMALaConeSmear;
    fh2MC2K0Cone                    = o.fh2MC2K0Cone;
    fh2MC2LaCone                    = o.fh2MC2LaCone;
    fh2MC2ALaCone                   = o.fh2MC2ALaCone;
    /*
    fh2MCEtaVsPtHijingLa            = o.fh2MCEtaVsPtHijingLa;
    fh2MCEtaVsPtInjectLa            = o.fh2MCEtaVsPtInjectLa;
    fh2MCEtaVsPtHijingALa           = o.fh2MCEtaVsPtHijingALa;
    fh2MCEtaVsPtInjectALa           = o.fh2MCEtaVsPtInjectALa;
    fhnrecMCHijingLaIncl            = o.fhnrecMCHijingLaIncl;
    fhnrecMCHijingLaCone            = o.fhnrecMCHijingLaCone;
    fhnrecMCHijingALaIncl           = o.fhnrecMCHijingALaIncl;
    fhnrecMCHijingALaCone           = o.fhnrecMCHijingALaCone;
    fhnrecMCInjectLaIncl            = o.fhnrecMCInjectLaIncl;
    fhnrecMCInjectLaCone            = o.fhnrecMCInjectLaCone;
    fhnrecMCInjectALaIncl           = o.fhnrecMCInjectALaIncl;
    fhnrecMCInjectALaCone           = o.fhnrecMCInjectALaCone;     */   
    fhnMCrecK0Cone                  = o.fhnMCrecK0Cone;
    fhnMCrecLaCone                  = o.fhnMCrecLaCone;
    fhnMCrecALaCone                 = o.fhnMCrecALaCone;
    /* fhnMCrecK0ConeSmear             = o.fhnMCrecK0ConeSmear;
    fhnMCrecLaConeSmear             = o.fhnMCrecLaConeSmear;
    fhnMCrecALaConeSmear            = o.fhnMCrecALaConeSmear;
    fhnK0sSecContinCone             = o.fhnK0sSecContinCone;
    fhnLaSecContinCone              = o.fhnLaSecContinCone;
    fhnALaSecContinCone             = o.fhnALaSecContinCone;*/
    fhnK0sIncl                      = o.fhnK0sIncl; 
    fhnK0sCone                      = o.fhnK0sCone;
    fhnK0sEmbCone                   = o.fhnK0sEmbCone;
    fhnK0sEmbConeRef                = o.fhnK0sEmbConeRef;
    fhnK0sEmbConeStandard           = o.fhnK0sEmbConeStandard;
    fhnLaIncl                       = o.fhnLaIncl;
    fhnLaCone                       = o.fhnLaCone;
    fhnLaEmbCone                    = o.fhnLaEmbCone;
    fhnLaEmbConeRef                 = o.fhnLaEmbConeRef;
    fhnLaEmbConeStandard            = o.fhnLaEmbConeStandard;
    fhnALaIncl                      = o.fhnALaIncl;
    fhnALaCone                      = o.fhnALaCone; 
    fhnALaEmbCone                   = o.fhnALaEmbCone;  
    fhnALaEmbConeRef                = o.fhnALaEmbConeRef;
    fhnALaEmbConeStandard           = o.fhnALaEmbConeStandard;

    fh2MCEmbK0sJetPt                = o.fh2MCEmbK0sJetPt;
    fh2MCEmbLaJetPt                 = o.fh2MCEmbLaJetPt;
    fh2MCEmbALaJetPt                = o.fh2MCEmbALaJetPt;

    fhnK0sPC                        = o.fhnK0sPC;
    fhnK0sEmbPC                     = o.fhnK0sEmbPC;
    fhnLaPC                         = o.fhnLaPC;
    fhnLaEmbPC                      = o.fhnLaEmbPC;
    fhnALaPC                        = o.fhnALaPC;
    fhnALaEmbPC                     = o.fhnALaEmbPC;
    fhnK0sRC                        = o.fhnK0sRC;
    fhnLaRC                         = o.fhnLaRC;
    fhnALaRC                        = o.fhnALaRC;
    fhnK0sRCBias                    = o.fhnK0sRCBias;
    fhnLaRCBias                     = o.fhnLaRCBias;
    fhnALaRCBias                    = o.fhnALaRCBias;
    fhnK0sOC                        = o.fhnK0sOC;
    fhnLaOC                         = o.fhnLaOC;
    fhnALaOC                        = o.fhnALaOC;
    fh1AreaExcluded                 = o.fh1AreaExcluded;
    fh1MedianEta                    = o.fh1MedianEta;
    fh1JetPtMedian                  = o.fh1JetPtMedian;
    fh1MCMultiplicityPrimary        = o.fh1MCMultiplicityPrimary;
    fh1MCMultiplicityTracks         = o.fh1MCMultiplicityTracks;
    fhnFeedDownLa                   = o.fhnFeedDownLa;
    fhnFeedDownALa                  = o.fhnFeedDownALa;
    fhnFeedDownLaCone               = o.fhnFeedDownLaCone;
    fhnFeedDownALaCone              = o.fhnFeedDownALaCone;
    fh2FeedDownXiLa                 = o.fh2FeedDownXiLa;
    fh2FeedDownXiALa                = o.fh2FeedDownXiALa;
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
    //fh1MCRapK0s                     = o.fh1MCRapK0s; 
    //fh1MCRapLambda                  = o.fh1MCRapLambda;
    //fh1MCRapAntiLambda              = o.fh1MCRapAntiLambda;
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

  if(jettracklist) delete jettracklist;
  if(jetConeK0list) delete jetConeK0list;
  if(jetConeLalist) delete jetConeLalist;
  if(jetConeALalist) delete jetConeALalist;
  if(jetConeK0Emblist) delete jetConeK0Emblist;
  if(jetConeLaEmblist) delete jetConeLaEmblist;
  if(jetConeALaEmblist) delete jetConeALaEmblist;
  if(jetConeK0EmbStlist) delete jetConeK0EmbStlist;
  if(jetConeLaEmbStlist) delete jetConeLaEmbStlist;
  if(jetConeALaEmbStlist) delete jetConeALaEmbStlist;

  if(jetConeK0EmbMClist) delete jetConeK0EmbMClist;
  if(jetConeLaEmbMClist) delete jetConeLaEmbMClist;
  if(jetConeALaEmbMClist) delete jetConeALaEmbMClist;

  if(jetPerpConeK0list) delete jetPerpConeK0list;
  if(jetPerpRecCutslist) delete jetPerpRecCutslist;
  if(jetPerpConeK0Emblist) delete jetPerpConeK0Emblist;
  if(jetPerpConeLalist) delete jetPerpConeLalist;
  if(jetPerpConeLaEmblist) delete jetPerpConeLaEmblist;
  if(jetPerpConeALalist) delete jetPerpConeALalist;
  if(jetPerpConeALaEmblist) delete jetPerpConeALaEmblist;
  if(jetMedianConeK0list) delete jetMedianConeK0list;
  if(jetMedianConeLalist) delete jetMedianConeLalist;
  if(jetMedianConeALalist) delete jetMedianConeALalist;
  if(jetMedianRecCutslist) delete jetMedianRecCutslist;
  if(fListK0sRC) delete fListK0sRC;
  if(fListLaRC) delete fListLaRC;
  if(fListALaRC) delete fListALaRC;
  if(fTracksRecCutsRC) delete fTracksRecCutsRC;
  if(fTracksRecBckgCuts) delete fTracksRecBckgCuts;
  if(fTracksPerpCone) delete fTracksPerpCone;
  if(fListK0s) delete fListK0s;
  if(fListK0sMC) delete fListK0sMC;
  if(fListLa) delete fListLa;
  if(fListLaMC) delete fListLaMC;
  if(fListALa) delete fListALa;
  if(fListALaMC) delete fListALaMC;
  if(fListK0sStandard) delete fListK0sStandard;
  if(fListLaStandard) delete fListLaStandard;
  if(fListALaStandard) delete fListALaStandard;
  if(fListFeeddownLaCand) delete fListFeeddownLaCand;
  if(fListFeeddownALaCand) delete fListFeeddownALaCand;
  if(jetConeFDLalist) delete jetConeFDLalist;
  if(jetConeFDALalist) delete jetConeFDALalist;   
  if(fListMCgenK0s) delete fListMCgenK0s;
  if(fListMCgenLa) delete fListMCgenLa;
  if(fListMCgenALa) delete fListMCgenALa;
  if(fListMCgenK0sCone) delete fListMCgenK0sCone;
  if(fListMCgenLaCone) delete fListMCgenLaCone;
  if(fListMCgenALaCone) delete fListMCgenALaCone;
  if(fRandom) delete fRandom;
  
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
  // fill FF, don't use TH3F anymore use THnSparse instead to save memory
 
  if(incrementJetPt) fh1JetPt->Fill(jetPt);    
  //fh3TrackPt->Fill(jetPt,invM,trackPt);//Fill(x,y,z)
  invM = 0;
  Double_t z = 0.;
  if(jetPt>0) z = trackPt / jetPt;
  // Double_t xi = 0;
  //if(z>0) xi = TMath::Log(1/z);
  
  //fh3Xi->Fill(jetPt,invM,xi);
  //fh3Z->Fill(jetPt,invM,z);
}

//___________________________________________________________________________________
void AliAnalysisTaskJetChem::AliFragFuncHistosInvMass::AddToOutput(TList* list) const
{
  // add histos to list

  list->Add(fh1JetPt);
  //list->Add(fh3TrackPt);
  //list->Add(fh3Xi);
  //list->Add(fh3Z);
}


//____________________________________________________
void AliAnalysisTaskJetChem::UserCreateOutputObjects()
{
  // create output objects
   
  fRandom = new TRandom3(0);
  fRandom->SetSeed(0);

  if(fDebug > 1) Printf("AliAnalysisTaskJetChem::UserCreateOutputObjects()");
 
  // create list of tracks and jets 
  jettracklist = new TList();
  jettracklist->SetOwner(kFALSE);
  jetConeK0list = new TList();
  jetConeK0list->SetOwner(kFALSE);
  jetConeLalist = new TList();
  jetConeLalist->SetOwner(kFALSE);
  jetConeALalist = new TList();
  jetConeALalist->SetOwner(kFALSE);
  jetConeK0Emblist = new TList();
  jetConeK0Emblist->SetOwner(kFALSE);
  jetConeLaEmblist = new TList();
  jetConeLaEmblist->SetOwner(kFALSE);
  jetConeALaEmblist = new TList();
  jetConeALaEmblist->SetOwner(kFALSE);
  jetConeK0EmbStlist = new TList();
  jetConeK0EmbStlist->SetOwner(kFALSE);
  jetConeLaEmbStlist = new TList();
  jetConeLaEmbStlist->SetOwner(kFALSE);
  jetConeALaEmbStlist = new TList();
  jetConeALaEmbStlist->SetOwner(kFALSE);

  jetConeK0EmbMClist = new TList();
  jetConeK0EmbMClist->SetOwner(kFALSE);
  jetConeLaEmbMClist = new TList();
  jetConeLaEmbMClist->SetOwner(kFALSE);
  jetConeALaEmbMClist = new TList();
  jetConeALaEmbMClist->SetOwner(kFALSE);

  jetPerpConeK0list = new TList();
  jetPerpConeK0list->SetOwner(kFALSE);
  jetPerpRecCutslist = new TList();
  jetPerpRecCutslist->SetOwner(kFALSE);
  jetPerpConeK0Emblist = new TList();
  jetPerpConeK0Emblist->SetOwner(kFALSE);
  jetPerpConeLalist = new TList(); 
  jetPerpConeLalist->SetOwner(kFALSE);
  jetPerpConeLaEmblist = new TList(); 
  jetPerpConeLaEmblist->SetOwner(kFALSE);
  jetPerpConeALalist = new TList();
  jetPerpConeALalist->SetOwner(kFALSE);
  jetPerpConeALaEmblist = new TList();
  jetPerpConeALaEmblist->SetOwner(kFALSE);
  jetMedianConeK0list = new TList();
  jetMedianConeK0list->SetOwner(kFALSE);
  jetMedianConeLalist = new TList();
  jetMedianConeLalist->SetOwner(kFALSE);
  jetMedianConeALalist = new TList();
  jetMedianConeALalist->SetOwner(kFALSE);
  jetMedianRecCutslist = new TList();
  jetMedianRecCutslist->SetOwner(kFALSE);
  fListK0sRC = new TList();
  fListK0sRC->SetOwner(kFALSE);
  fListLaRC = new TList();
  fListLaRC->SetOwner(kFALSE);
  fListALaRC = new TList();
  fListALaRC->SetOwner(kFALSE);
  fTracksRecCutsRC = new TList();
  fTracksRecCutsRC->SetOwner(kFALSE);
  fTracksRecBckgCuts = new TList();
  fTracksRecBckgCuts->SetOwner(kFALSE);
  fTracksRecCuts = new TList();
  fTracksRecCuts->SetOwner(kFALSE); //objects in TList wont be deleted when TList is deleted 
  fTracksGen = new TList();
  fTracksGen->SetOwner(kFALSE);
  fTracksPerpCone = new TList();
  fTracksPerpCone->SetOwner(kFALSE); 
  fJetsRecCuts = new TList();
  fJetsRecCuts->SetOwner(kFALSE);
  fJetsGen = new TList();
  fJetsGen->SetOwner(kFALSE);
  fJetsEmbedded = new TList();
  fJetsEmbedded->SetOwner(kFALSE);
  fBckgJetsRec = new TList();
  fBckgJetsRec->SetOwner(kFALSE);
  fListK0s = new TList(); 
  fListK0s->SetOwner(kFALSE);
  fListK0sMC = new TList(); 
  fListK0sMC->SetOwner(kFALSE);
  fListLa = new TList(); 
  fListLa->SetOwner(kFALSE);
  fListLaMC = new TList(); 
  fListLaMC->SetOwner(kFALSE);
  fListALa = new TList(); 
  fListALa->SetOwner(kFALSE);
  fListALaMC = new TList(); 
  fListALaMC->SetOwner(kFALSE);
  fListK0sStandard = new TList(); 
  fListK0sStandard->SetOwner(kFALSE);
  fListLaStandard = new TList(); 
  fListLaStandard->SetOwner(kFALSE);
  fListALaStandard = new TList(); 
  fListALaStandard->SetOwner(kFALSE);
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
  fListMCgenK0sCone = new TList();
  fListMCgenK0sCone->SetOwner(kFALSE);
  fListMCgenLaCone = new TList();
  fListMCgenLaCone->SetOwner(kFALSE);
  fListMCgenALaCone = new TList();
  fListMCgenALaCone->SetOwner(kFALSE);
  
  // Create histograms / output container
 
  fCommonHistList = new TList();
  fCommonHistList->SetOwner();
  
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
  fh1nRecJetsCuts            = new TH1F("fh1nRecJetsCuts","reconstructed jets per event",100,-0.5,99.5);
 
  // histograms JetChem task
 
  fh1EvtAllCent	                = new TH1F("fh1EvtAllCent","before centrality selection",100,0.,100.);
  fh1Evt                        = new TH1F("fh1Evt", "All events runned over", 3, 0.,1.);
  fh1EvtMult 	                = new TH1F("fh1EvtMult","multiplicity",240,0.,240.);
  fh1K0Mult 	                = new TH1F("fh1K0Mult","K0 multiplicity",100,0.,100.);//500. all
  fh1dPhiJetK0                  = new TH1F("fh1dPhiJetK0","",64,-1,5.4);
  fh1LaMult 	                = new TH1F("fh1LaMult","La multiplicity",100,0.,100.);
  fh1dPhiJetLa                  = new TH1F("fh1dPhiJetLa","",64,-1,5.4);
  fh1ALaMult 	                = new TH1F("fh1ALaMult","ALa multiplicity",100,0.,100.);
  fh1dPhiJetALa                 = new TH1F("fh1dPhiJetALa","",64,-1,5.4);
  fh1JetEta                     = new TH1F("fh1JetEta","#eta distribution of all jets",40,-2.,2.);
  fh1JetPhi                     = new TH1F("fh1JetPhi","#phi distribution of all jets",63,0.,6.3);
  fh2JetEtaPhi                  = new TH2F("fh2JetEtaPhi","#eta and #phi distribution of all jets",400,-2.,2.,63,0.,6.3);
  
  if(fDebug>2)std::cout<<"fBranchEmbeddedJets.Length(): "<<fBranchEmbeddedJets.Length()<<std::endl;
  if(fDebug>2)std::cout<<"fBranchGenJets.Length(): "<<fBranchGenJets.Length()<<std::endl;


  //embedding
  fh1nEmbeddedJets              = new TH1F("fh1nEmbeddedJets","Number of embedded jets",10,-0.5,9.5);
  fh1nGenJets                   = new TH1F("fh1nGenJets","generated jets per event",10,-0.5,9.5);
  fh2TracksPerpCone             = new TH2F("fh2TracksPerpCone","Charged tracks in 2 perp. cones;#it{p^{ch,jet}}_{T} (GeV/#it{c});#it{p^{ch}}_{T} (GeV/#it{c})",19,5.,100.,120,0.,12.);
  fh1IndexEmbedded              = new TH1F("fh1IndexEmbedded","",11,-1.,10.);

  fh1PtEmbBeforeMatch           = new TH1F("fh1PtEmbBeforeMatch","Pt spectrum of jets before JetMatching",19,5.,100.);
  fh1PtEmbExtraOnly             = new TH1F("fh1PtEmbExtraOnly","Pt spectrum of jets from ExtraOnly tracks (embedded truth)",19,5.,100.);
  
  fh1PtEmbReject                = new TH1F("fh1PtEmbReject","Pt spectrum of jets rejected by JetMatching cuts",19,5.,100.);
  fh2PtEtaEmbReject             = new TH2F("fh2PtEtaEmbReject","Pt #eta distribution of jets rejected by JetMatching cuts #eta; #it{p}_{T}",19,5.,100.,200,-1.,1.);
  fh1PtEmbAfterMatch            = new TH1F("fh1PtEmbAfterMatch","Pt spectrum of jets after JetMatching cuts and with leading constituent cut",19,5.,100.);
 
  fh1FractionPtEmbedded         = new TH1F("fh1FractionPtEmbedded","",110,0.,1.1);
  fh1DeltaREmbedded             = new TH1F("fh1DeltaREmbedded","",50,0.,0.5);
  
  fh1IndexEmbeddedMC            = new TH1F("fh1IndexEmbeddedMC","",11,-1.,10.);
  fh1FractionPtEmbeddedMC         = new TH1F("fh1FractionPtEmbeddedMC","",110,0,1.1);
  fh2FractionPtVsEmbeddedJetPtMC  = new TH2F("fh2FractionPtVsEmbeddedJetPtMC","",250,0,250,110,0,1.1);
  fh1DeltaREmbeddedMC             = new TH1F("fh1DeltaREmbeddedMC","",50,0,0.5);
  fh1JetPtEmbGenAfterMatch        = new TH1F("fh1JetPtEmbGenAfterMatch","Pt spectrum of jets after JetMatching cuts and with leading constituent cut",19,5.,100.);
  
  //


  fh1PerpCone                   = new TH1F("fh1PerpCone","Number of perp. cones for charged tracks in event",2.,0.5,1.5);
  fh1V0PtCandidate              = new TH1F("fh1V0PtCandidate","p_{T} distribution of all v0s candidates of PYTHIA",200,0.,200.);
  fh1IMK0Cone                   = new TH1F("fh1IMK0Cone","p_{T} distribution of all jets containing K0s candidates",19,5.,100.);
  fh1IMLaCone                   = new TH1F("fh1IMLaCone","p_{T} distribution of all jets containing #Lambda candidates",19,5.,100.);
  fh1IMALaCone                  = new TH1F("fh1IMALaCone","p_{T} distribution of all jets containing #bar{#Lambda} candidates",19,5.,100.);
  fh1IMK0EmbCone                = new TH1F("fh1IMK0EmbCone","p_{T} distribution of all embedded and selected jets containing K0s candidates",19,5.,100.);
  fh1IMLaEmbCone                = new TH1F("fh1IMLaEmbCone","p_{T} distribution of all embedded and selected jets containing #Lambda candidates",19,5.,100.);
  fh1IMALaEmbCone               = new TH1F("fh1IMALaEmbCone","p_{T} distribution of all embedded and selected jets containing #bar{#Lambda} candidates",19,5.,100.);
  fh2FFJetTrackEta              = new TH2F("fh2FFJetTrackEta","charged track eta distr. in jet cone",200,-1.,1.,40,0.,200.);  
  //fh1trackPosNCls             = new TH1F("fh1trackPosNCls","NTPC clusters positive daughters",10,0.,100.);
  //fh1trackNegNCls             = new TH1F("fh1trackNegNCls","NTPC clusters negative daughters",10,0.,100.);
  fh1trackPosEta                = new TH1F("fh1trackPosEta","eta positive daughters",100,-2.,2.);
  fh1trackNegEta                = new TH1F("fh1trackNegEta","eta negative daughters",100,-2.,2.);
  fh1V0Eta                      = new TH1F("fh1V0Eta","V0 eta",60,-1.5,1.5);
  //fh1V0totMom                 = new TH1F("fh1V0totMom","V0 tot mom",100,0.,20.); 
  fh1CosPointAngle              = new TH1F("fh1CosPointAngle", "Cosine of V0's pointing angle",50,0.99,1.0);
  fh1DecayLengthV0              = new TH1F("fh1DecayLengthV0", "V0s decay Length;decay length(cm)",1200,0.,120.);
  fh2ProperLifetimeK0sVsPtBeforeCut = new TH2F("fh2ProperLifetimeK0sVsPtBeforeCut"," K0s ProperLifetime vs Pt; p_{T} (GeV/#it{c})",150,0.,15.,250,0.,250.);
  fh2ProperLifetimeK0sVsPtAfterCut = new TH2F("fh2ProperLifetimeK0sVsPtAfterCut"," K0s ProperLifetime vs Pt; p_{T} (GeV/#it{c})",150,0.,15.,250,0.,250.);
  fh1V0Radius                   = new TH1F("fh1V0Radius", "V0s Radius;Radius(cm)",200,0.,40.);
  fh1DcaV0Daughters             = new TH1F("fh1DcaV0Daughters", "DCA between daughters;dca(cm)",200,0.,2.);
  fh1DcaPosToPrimVertex         = new TH1F("fh1DcaPosToPrimVertex", "Positive V0 daughter;dca(cm)",100,0.,10.);
  fh1DcaNegToPrimVertex         = new TH1F("fh1DcaNegToPrimVertex", "Negative V0 daughter;dca(cm)",100,0.,10.);
  fh2ArmenterosBeforeCuts       = new TH2F("fh2ArmenterosBeforeCuts","Armenteros Podolanski Plot for K0s Candidates;#alpha;(p^{arm})_{T}/(GeV/#it{c})",200,-1.2,1.2,350,0.,0.35);
  fh2ArmenterosAfterCuts        = new TH2F("fh2ArmenterosAfterCuts","Armenteros Podolanski Plot for K0s Candidates;#alpha;(p^{arm})_{T}/(GeV/#it{c});",200,-1.2,1.2,350,0.,0.35);
  fh2BBLaPos                    = new TH2F("fh2BBLaPos","PID of the positive daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",100,0,10,200,0,200);
  fh2BBLaNeg                    = new TH2F("fh2BBLaNeg","PID of the negative daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",100,0,10,200,0,200);
  fh1PosDaughterCharge          = new TH1F("fh1PosDaughterCharge","charge of V0 positive daughters; V0 daughters",3,-2.,2.);
  fh1NegDaughterCharge          = new TH1F("fh1NegDaughterCharge","charge of V0 negative daughters; V0 daughters",3,-2.,2.);
  fh1PtMCK0s                    = new TH1F("fh1PtMCK0s","Pt of MC rec K0s; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1PtMCLa                     = new TH1F("fh1PtMCLa","Pt of MC rec La; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1PtMCALa                    = new TH1F("fh1PtMCALa","Pt of MC rec ALa; #it{p}_{T} (GeV/#it{c})",200,0.,20.);
  fh1EtaK0s                     = new TH1F("fh1EtaK0s","K^{0}_{s} entries ;#eta",200,-1.,1.);
  fh1EtaLa                      = new TH1F("fh1EtaLa","#Lambda entries ;#eta",200,-1.,1.);
  fh1EtaALa                     = new TH1F("fh1EtaALa","#bar{#Lambda} entries ;#eta",200,-1.,1.);

  //histos for normalisation of MCC, RC, OC and NJ

  fh1RC                         = new TH1F("fh1RC"," # random cones used",1,0.5,1.5);
  fh1RCBiasK0                   = new TH1F("fh1RCBiasK0"," # random cones with K0s trigger particle",1,0.5,1.5);
  fh1RCBiasLa                   = new TH1F("fh1RCBiasLa"," # random cones with La trigger particle",1,0.5,1.5);
  fh1RCBiasALa                  = new TH1F("fh1RCBiasALa"," # random cones with ALa trigger particle",1,0.5,1.5);
  fh1MCC                        = new TH1F("fh1MCC","# median cluster cones used",1,0.5,1.5);
  fh1OC                         = new TH1F("fh1OC","# outside cones used, number of jet events",1,0.5,1.5);
  fh1NJ                         = new TH1F("fh1NJ","# non-jet events used",1,0.5,1.5);
  fh1NJEmbEvt                   = new TH1F("fh1NJEmbEvt","# Embedding non-jet events used",1,0.5,1.5);
  fh1BckgJets                   = new TH1F("fh1BckgJets","Jet pT distribution of bckg jets (anti-kt data jets) used in Embedding study",19,5.,100.);
  fh1BckgJetsPtBias             = new TH1F("fh1BckgJetsPtBias","Jet pT distribution after JetPtBias of bckg jets (anti-kt data jets) used in Embedding study",19,5.,100.);


  Int_t binsInvMassEtaTrackPtK0s[3] = {200, 200, 120};//eta,invM,trackPt
  Double_t xminInvMassEtaTrackPtK0s[3] = {-1.,0.3,0.};
  Double_t xmaxInvMassEtaTrackPtK0s[3] = {1.,0.7,12.};

  fhnInvMassEtaTrackPtK0s       = new THnSparseF("fhnInvMassEtaTrackPtK0s","#eta; K0s invM (GeV/{#it{c}}^{2}); #it{p}_{T} (GeV/#it{c})",3,binsInvMassEtaTrackPtK0s,xminInvMassEtaTrackPtK0s,xmaxInvMassEtaTrackPtK0s);

  Int_t binsInvMassEtaTrackPtLa[3] = {200, 200, 120};//eta,invM,trackPt
  Double_t xminInvMassEtaTrackPtLa[3] = {-1.,1.05,0.};
  Double_t xmaxInvMassEtaTrackPtLa[3] = {1.,1.25,12.};

  fhnInvMassEtaTrackPtLa       = new THnSparseF("fhnInvMassEtaTrackPtLa","#eta; #Lambda invM (GeV/{#it{c}}^{2}); #it{p}_{T} (GeV/#it{c})",3,binsInvMassEtaTrackPtLa,xminInvMassEtaTrackPtLa,xmaxInvMassEtaTrackPtLa);

  Int_t binsInvMassEtaTrackPtALa[3] = {200, 200, 120};//eta,invM,trackPt
  Double_t xminInvMassEtaTrackPtALa[3] = {-1.,1.05,0.};
  Double_t xmaxInvMassEtaTrackPtALa[3] = {1.,1.25,12.};

  fhnInvMassEtaTrackPtALa       = new THnSparseF("fhnInvMassEtaTrackPtALa","#eta; #bar{#Lambda} invM (GeV/{#it{c}}^{2}); #it{p}_{T} (GeV/#it{c})",3,binsInvMassEtaTrackPtALa,xminInvMassEtaTrackPtALa,xmaxInvMassEtaTrackPtALa);

  Int_t binsK0sPC[4] = {19, 200, 120, 200};
  Double_t xminK0sPC[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sPC[4] = {100.,0.7, 12., 1.};
  fhnK0sPC                      = new THnSparseF("fhnK0sPC","two perp cones;jet pT; K0s invM; particle pT; particle #eta",4,binsK0sPC,xminK0sPC,xmaxK0sPC);

  Int_t binsLaPC[4] = {19, 200, 120, 200};
  Double_t xminLaPC[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaPC[4] = {100.,1.25, 12., 1.};
  fhnLaPC                       = new THnSparseF("fhnLaPC","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",4,binsLaPC,xminLaPC,xmaxLaPC);

  Int_t binsALaPC[4] = {19, 200, 120, 200};
  Double_t xminALaPC[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaPC[4] = {100.,1.25, 12., 1.};
  fhnALaPC                      = new THnSparseF("fhnALaPC","two perp cones;jet pT; #bar#Lambda invM; particle pT; particle #eta",4,binsALaPC,xminALaPC,xmaxALaPC);

  Int_t binsK0sEmbPC[4] = {19, 200, 120, 200};
  Double_t xminK0sEmbPC[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sEmbPC[4] = {100.,0.7, 12., 1.};
  fhnK0sEmbPC                      = new THnSparseF("fhnK0sEmbPC","two perp cones;jet pT; K0s invM; particle pT; particle #eta",4,binsK0sEmbPC,xminK0sEmbPC,xmaxK0sEmbPC);

  Int_t binsLaEmbPC[4] = {19, 200, 120, 200};
  Double_t xminLaEmbPC[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaEmbPC[4] = {100.,1.25, 12., 1.};
  fhnLaEmbPC                       = new THnSparseF("fhnLaEmbPC","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",4,binsLaEmbPC,xminLaEmbPC,xmaxLaEmbPC);

  Int_t binsALaEmbPC[4] = {19, 200, 120, 200};
  Double_t xminALaEmbPC[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaEmbPC[4] = {100.,1.25, 12., 1.};
  fhnALaEmbPC                      = new THnSparseF("fhnALaEmbPC","two perp cones;jet pT; #bar#Lambda invM; particle pT; particle #eta",4,binsALaEmbPC,xminALaEmbPC,xmaxALaEmbPC);
  
  Int_t binsK0sMCC[3] = {200, 120, 200};
  Double_t xminK0sMCC[3] = {0.3, 0., -1.};
  Double_t xmaxK0sMCC[3] = {0.7, 12., 1.};
  fhnK0sMCC                     = new THnSparseF("fhnK0sMCC","two perp cones;jet pT; K0s invM; particle pT; particle #eta",3,binsK0sMCC,xminK0sMCC,xmaxK0sMCC);

  Int_t binsLaMCC[3] = {200, 120, 200};
  Double_t xminLaMCC[3] = {1.05, 0., -1.};
  Double_t xmaxLaMCC[3] = {1.25, 12., 1.};
  fhnLaMCC                      = new THnSparseF("fhnLaMCC","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",3,binsLaMCC,xminLaMCC,xmaxLaMCC);

  Int_t binsALaMCC[3] = {200, 120, 200};
  Double_t xminALaMCC[3] = {1.05, 0., -1.};
  Double_t xmaxALaMCC[3] = {1.25, 12., 1.};
  fhnALaMCC                = new THnSparseF("fhnALaMCC","two perp cones;jet pT; #bara#Lambda invM; particle pT; particle #eta",3,binsALaMCC,xminALaMCC,xmaxALaMCC);

  Int_t binsK0sRC[3] = {200, 120, 200};
  Double_t xminK0sRC[3] = {0.3, 0., -1.};
  Double_t xmaxK0sRC[3] = {0.7, 12., 1.};
  fhnK0sRC                 = new THnSparseF("fhnK0sRC","two perp cones;jet pT; K0s invM; particle pT; particle #eta",3,binsK0sRC,xminK0sRC,xmaxK0sRC);

  Int_t binsLaRC[3] = {200, 120, 200};
  Double_t xminLaRC[3] = {1.05, 0., -1.};
  Double_t xmaxLaRC[3] = {1.25, 12., 1.};
  fhnLaRC                  = new THnSparseF("fhnLaRC","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",3,binsLaRC,xminLaRC,xmaxLaRC);

  Int_t binsALaRC[3] = {200, 120, 200};
  Double_t xminALaRC[3] = {1.05, 0., -1.};
  Double_t xmaxALaRC[3] = {1.25, 12., 1.};
  fhnALaRC                 = new THnSparseF("fhnALaRC","two perp cones;jet pT; #bara#Lambda invM; particle pT; particle #eta",3,binsALaRC,xminALaRC,xmaxALaRC);

  Int_t binsK0sRCBias[3] = {200, 120, 200};
  Double_t xminK0sRCBias[3] = {0.3, 0., -1.};
  Double_t xmaxK0sRCBias[3] = {0.7, 12., 1.};
  fhnK0sRCBias             = new THnSparseF("fhnK0sRCBias","two perp cones;jet pT; K0s invM; particle pT; particle #eta",3,binsK0sRCBias,xminK0sRCBias,xmaxK0sRCBias);

  Int_t binsLaRCBias[3] = {200, 120, 200};
  Double_t xminLaRCBias[3] = {1.05, 0., -1.};
  Double_t xmaxLaRCBias[3] = {1.25, 12., 1.};
  fhnLaRCBias              = new THnSparseF("fhnLaRCBias","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",3,binsLaRCBias,xminLaRCBias,xmaxLaRCBias);

  Int_t binsALaRCBias[3] = {200, 120, 200};
  Double_t xminALaRCBias[3] = {1.05, 0., -1.};
  Double_t xmaxALaRCBias[3] = {1.25, 12., 1.};
  fhnALaRCBias             = new THnSparseF("fhnALaRCBias","two perp cones;jet pT; #bara#Lambda invM; particle pT; particle #eta",3,binsALaRCBias,xminALaRCBias,xmaxALaRCBias);

  Int_t binsK0sOC[3] = {200, 120, 200};
  Double_t xminK0sOC[3] = {0.3, 0., -1.};
  Double_t xmaxK0sOC[3] = {0.7, 12., 1.};
  fhnK0sOC                     = new THnSparseF("fhnK0sOC","two perp cones;jet pT; K0s invM; particle pT; particle #eta",3,binsK0sOC,xminK0sOC,xmaxK0sOC);

  Int_t binsLaOC[3] = {200, 120, 200};
  Double_t xminLaOC[3] = {1.05, 0., -1.};
  Double_t xmaxLaOC[3] = {1.25, 12., 1.};
  fhnLaOC                      = new THnSparseF("fhnLaOC","two perp cones;jet pT; #Lambda invM; particle pT; particle #eta",3,binsLaOC,xminLaOC,xmaxLaOC);

  Int_t binsALaOC[3] = {200, 120, 200};
  Double_t xminALaOC[3] = {1.05, 0., -1.};
  Double_t xmaxALaOC[3] = {1.25, 12., 1.};

  fhnALaOC                      = new THnSparseF("fhnALaOC","two perp cones;jet pT; #bara#Lambda invM; particle pT; particle #eta",3,binsALaOC,xminALaOC,xmaxALaOC);

  fh1AreaExcluded               = new TH1F("fh1AreaExcluded","area excluded for selected jets in event acceptance",50,0.,1.);

  fh1MedianEta                  = new TH1F("fh1MedianEta","Median cluster axis ;#eta",200,-1.,1.);
  fh1JetPtMedian                = new TH1F("fh1JetPtMedian"," (selected) jet it{p}_{T} distribution for MCC method; #GeV/it{c}",19,5.,100.);

  fh1TrackMultCone              = new TH1F("fh1TrackMultCone","track multiplicity in jet cone; number of tracks",20,0.,50.);

  fh2TrackMultCone              = new TH2F("fh2TrackMultCone","track multiplicity in jet cone vs. jet momentum; number of tracks; jet it{p}_{T} (GeV/it{c})",50,0.,50.,19,5.,100.);

  Int_t binsNJK0[3] = {200, 120, 200};
  Double_t xminNJK0[3] = {0.3, 0., -1.};
  Double_t xmaxNJK0[3] = {0.7, 12., 1.};
  fhnNJK0                       = new THnSparseF("fhnNJK0","K0s candidates in events wo selected jets;",3,binsNJK0,xminNJK0,xmaxNJK0);

  Int_t binsNJLa[3] = {200, 120, 200};
  Double_t xminNJLa[3] = {1.05, 0., -1.};
  Double_t xmaxNJLa[3] = {1.25, 12., 1.};
  fhnNJLa                      = new THnSparseF("fhnNJLa","La candidates in events wo selected jets; ",3,binsNJLa,xminNJLa,xmaxNJLa);

  Int_t binsNJALa[3] = {200, 120, 200};
  Double_t xminNJALa[3] = {1.05, 0., -1.};
  Double_t xmaxNJALa[3] = {1.25, 12., 1.};
  fhnNJALa                     = new THnSparseF("fhnNJALa","ALa candidates in events wo selected jets; ",3,binsNJALa,xminNJALa,xmaxNJALa);

  //fh2ChTracksNJ                = new TH2F("fh2ChTracksNJ","charged tracks in non-jet events; it{p}_{T} (GeV/it{c};#eta)",120,0.,12.,200,-1.,1.);

  //fh2ChTracksRC                = new TH2F("fh2ChTracksRC","charged tracks in random cones; it{p}_{T} (GeV/it{c};#eta)",120,0.,12.,200,-1.,1.);
  //fh2ChTracksOC                = new TH2F("fh2ChTracksOC","charged tracks outside cones; it{p}_{T} (GeV/it{c};#eta)",120,0.,12.,200,-1.,1.);

  //fh2ChTracksMCC               = new TH2F("fh2ChTracksMCC","charged tracks in non-jet events; it{p}_{T} (GeV/it{c};#eta)",120,0.,12.,200,-1.,1.);

  //fh2ChTracksPC                = new TH2F("fh2ChTracksPC","charged tracks in perpendicular cones; it{p}_{T} (GeV/it{c};#eta)",120,0.,12.,200,-1.,1.);

  fFFHistosRecCuts   	        = new AliFragFuncHistos("RecCuts", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  fV0QAK0                       = new AliFragFuncQATrackHistos("V0QAK0",fQATrackNBinsPt, fQATrackPtMin, fQATrackPtMax, 
							    fQATrackNBinsEta, fQATrackEtaMin, fQATrackEtaMax,
							    fQATrackNBinsPhi, fQATrackPhiMin, fQATrackPhiMax, 
							    fQATrackHighPtThreshold);
  
  fFFHistosRecCutsK0Evt         = new AliFragFuncHistos("RecCutsK0Evt", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						     fFFNBinsPt, fFFPtMin, fFFPtMax, 
						     fFFNBinsXi, fFFXiMin, fFFXiMax,  
						     fFFNBinsZ , fFFZMin , fFFZMax);
  
  fFFHistosGen   	      = new AliFragFuncHistos("Gen", fFFNBinsJetPt, fFFJetPtMin, fFFJetPtMax, 
						      fFFNBinsPt, fFFPtMin, fFFPtMax, 
						      fFFNBinsXi, fFFXiMin, fFFXiMax,  
						      fFFNBinsZ , fFFZMin , fFFZMax);
 
  //***************
  // MC histograms
  //***************

  //fh2MCgenK0Cone                = new TH2F("fh2MCgenK0Cone", "MC gen {K^{0}}^{s} #it{p}_{T}  in cone around rec jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",19,5.,100.,200,0.,20.);
  //fh2MCgenLaCone                = new TH2F("fh2MCgenLaCone", "MC gen #Lambda #it{p}_{T} in cone around rec jet axis versus jet #it{p}_{T} ; jet #it{p}_{T}",19,5.,100.,200,0.,20.);
  //fh2MCgenALaCone               = new TH2F("fh2MCgenALaCone", "MC gen #Antilambda #it{p}_{T} in cone around rec jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",19,5.,100.,200,0.,20.);

  Int_t binsMCgenK0Cone[4] = {19, 120, 200};
  Double_t xminMCgenK0Cone[4] = {5., 0., -1.};
  Double_t xmaxMCgenK0Cone[4] = {100., 12., 1.};
  fhnMCgenK0Cone                = new THnSparseF("fhnMCgenK0Cone", "MC gen {K^{0}}^{s} #it{p}_{T}  in cone around jet axis",3,binsMCgenK0Cone,xminMCgenK0Cone,xmaxMCgenK0Cone);  


  Int_t binsMCgenLaCone[4] = {19, 120, 200};
  Double_t xminMCgenLaCone[4] = {5., 0., -1.};
  Double_t xmaxMCgenLaCone[4] = {100., 12., 1.};
  fhnMCgenLaCone                = new THnSparseF("fhnMCgenLaCone", "MC gen #Lambda #it{p}_{T}  in cone around jet axis",3,binsMCgenLaCone,xminMCgenLaCone,xmaxMCgenLaCone);  


  //fh2MCgenK0Cone->GetYaxis()->SetTitle("MC gen K^{0}}^{s} #it{p}_{T}");
  //fh2MCgenLaCone->GetYaxis()->SetTitle("MC gen #Lambda #it{p}_{T}");
  //fh2MCgenALaCone->GetYaxis()->SetTitle("MC gen #Antilambda #it{p}_{T}");

  //fh2MCEtagenK0Cone             = new TH2F("fh2MCEtagenK0Cone","MC gen {K^{0}}^{s} #it{p}_{T} #eta distribution in jet cone;#eta",19,5.,100.,200,-1.,1.);
  //fh2MCEtagenLaCone             = new TH2F("fh2MCEtagenLaCone","MC gen #Lambda #it{p}_{T} #eta distribution in jet cone;#eta",19,5.,100.,200,-1.,1.);
  //fh2MCEtagenALaCone            = new TH2F("fh2MCEtagenALaCone","MC gen #Antilambda #it{p}_{T} #eta distribution in jet cone;#eta",19,5.,100.,200,-1.,1.);
  
  fh1IMK0ConeSmear                = new TH1F("fh1IMK0ConeSmear","Smeared jet pt study for K0s-in-cone-jets; smeared jet #it{p}_{T}", 19,5.,100.);
  fh1IMLaConeSmear                = new TH1F("fh1IMLaConeSmear","Smeared jet pt study for La-in-cone-jets; smeared jet #it{p}_{T}", 19,5.,100.);
  fh1IMALaConeSmear               = new TH1F("fh1IMALaConeSmear","Smeared jet pt study for ALa-in-cone-jets; smeared jet #it{p}_{T}", 19,5.,100.);

  fh2MC2K0Cone                    = new TH2F("fh2MC2K0Cone", "MC true {K^{0}}^{s} #it{p}_{T} in cone around jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",19,5.,100.,120,0.,12.);
  fh2MC2LaCone                    = new TH2F("fh2MC2LaCone", "MC true {#Lambda #it{p}_{T} in cone around jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",19,5.,100.,120,0.,12.);
  fh2MC2ALaCone                    = new TH2F("fh2MC2ALaCone", "MC true {#bar{#Lambda} #it{p}_{T} in cone around jet axis versus jet #it{p}_{T}; jet #it{p}_{T}",19,5.,100.,120,0.,12.);

  /*
  fh2CorrHijingLaProton           = new TH2F("fh2CorrHijingLaProton","#Lambda - proton pT correlation, Hijing;#it{p^{#Lambda}}_{T} (GeV/c);#it{p^{proton}}_{T} (GeV/c)",120,0.,12.,120,0.,12.);        
  fh2CorrInjectLaProton           = new TH2F("fh2CorrInjectLaProton","#Lambda - proton pT correlation, Injected;#it{p^{#Lambda}}_{T} (GeV/c);#it{p^{proton}}_{T} (GeV/c)",120,0.,12.,120,0.,12.);
  fh2CorrHijingALaAProton         = new TH2F("fh2CorrHijingALaAProton","#bar{#Lambda} - proton pT correlation, Hijing;#it{p^{#Lambda}}_{T} (GeV/c);#it{p^{#bar{proton}}}_{T} (GeV/c)",120,0.,12.,120,0.,12.);        
  fh2CorrInjectALaAProton         = new TH2F("fh2CorrInjectALaAProton","#bar{#Lambda} - proton pT correlation, Injected;#it{p^{#Lambda}}_{T} (GeV/c);#it{p^{#bar{proton}}}_{T} (GeV/c)",120,0.,12.,120,0.,12.);
  //12 new histograms: Cone, Incl, Lambda, Antilambda, Hijing, Injected:
  
  fh2MCEtaVsPtHijingLa              = new TH2F("fh2MCEtaVsPtHijingLa","MC Hijing gen. #Lambda #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);
  fh2MCEtaVsPtInjectLa              = new TH2F("fh2MCEtaVsPtInjectLa","MC injected gen. #Lambda  #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);
  fh2MCEtaVsPtHijingALa             = new TH2F("fh2MCEtaVsPtHijingALa","MC gen. Hijing  #bar{#Lambda} #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);
  fh2MCEtaVsPtInjectALa             = new TH2F("fh2MCEtaVsPtInjectALa","MC gen. injected #bar{#Lambda}  #eta; #it{p}_{T}",200,0.,20.,200,-1.,1.);

  Int_t binsrecMCHijingLaIncl[3] = {200, 120, 200};
  Double_t xminrecMCHijingLaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxrecMCHijingLaIncl[3] = {1.25, 12., 1.};
  fhnrecMCHijingLaIncl          = new THnSparseF("fhnrecMCHijingLaIncl","La inv. mass; particle pT; particle #eta",3,binsrecMCHijingLaIncl,xminrecMCHijingLaIncl,xmaxrecMCHijingLaIncl);

  Int_t binsrecMCHijingLaCone[4] = {19, 200, 120, 200};
  Double_t xminrecMCHijingLaCone[4] = {5., 1.05, 0., -1.};
  Double_t xmaxrecMCHijingLaCone[4] = {100., 1.25, 12., 1.};
  fhnrecMCHijingLaCone          = new THnSparseF("fhnrecMCHijingLaCone","La inv. mass; particle pT; particle #eta",4,binsrecMCHijingLaCone,xminrecMCHijingLaCone,xmaxrecMCHijingLaCone);

  Int_t binsrecMCHijingALaIncl[3] = {200, 120, 200};
  Double_t xminrecMCHijingALaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxrecMCHijingALaIncl[3] = {1.25, 12., 1.};
  fhnrecMCHijingALaIncl         = new THnSparseF("fhnrecMCHijingALaIncl","ALa inv. mass; particle pT; particle #eta",3,binsrecMCHijingALaIncl,xminrecMCHijingALaIncl,xmaxrecMCHijingALaIncl);

  Int_t binsrecMCHijingALaCone[4] = {19, 200, 120, 200};
  Double_t xminrecMCHijingALaCone[4] = {5., 1.05, 0., -1.};
  Double_t xmaxrecMCHijingALaCone[4] = {100., 1.25, 12., 1.};
  fhnrecMCHijingALaCone         = new THnSparseF("fhnrecMCHijingALaCone","ALa inv. mass; particle pT; particle #eta",4,binsrecMCHijingALaCone,xminrecMCHijingALaCone,xmaxrecMCHijingALaCone);

  Int_t binsrecMCInjectLaIncl[3] = {200, 120, 200};
  Double_t xminrecMCInjectLaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxrecMCInjectLaIncl[3] = {1.25, 12., 1.};
  fhnrecMCInjectLaIncl          = new THnSparseF("fhnrecMCInjectLaIncl","La inv. mass; particle pT; particle #eta",3,binsrecMCInjectLaIncl,xminrecMCInjectLaIncl,xmaxrecMCInjectLaIncl);

  Int_t binsrecMCInjectLaCone[4] = {19, 200, 120, 200};
  Double_t xminrecMCInjectLaCone[4] = {5., 1.05, 0., -1.};
  Double_t xmaxrecMCInjectLaCone[4] = {100., 1.25, 12., 1.};
  fhnrecMCInjectLaCone          = new THnSparseF("fhnrecMCInjectLaCone","La jet pT;inv. mass; particle pT; particle #eta",4,binsrecMCInjectLaCone,xminrecMCInjectLaCone,xmaxrecMCInjectLaCone);

  Int_t binsrecMCInjectALaIncl[3] = {200, 120, 200};
  Double_t xminrecMCInjectALaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxrecMCInjectALaIncl[3] = {1.25, 12., 1.};
  fhnrecMCInjectALaIncl         = new THnSparseF("fhnrecMCInjectALaIncl","ALa inv. mass; particle pT; particle #eta",3,binsrecMCInjectALaIncl,xminrecMCInjectALaIncl,xmaxrecMCInjectALaIncl);

  Int_t binsrecMCInjectALaCone[4] = {19, 200, 120, 200};
  Double_t xminrecMCInjectALaCone[4] = {5., 1.05, 0., -1.};
  Double_t xmaxrecMCInjectALaCone[4] = {100., 1.25, 12., 1.};
  fhnrecMCInjectALaCone         = new THnSparseF("fhnrecMCInjectALaCone","ALa inv. mass; particle pT; particle #eta",4,binsrecMCInjectALaCone,xminrecMCInjectALaCone,xmaxrecMCInjectALaCone);
  */

  Int_t binsMCrecK0Cone[4] = {19, 200, 120, 200};
  Double_t xminMCrecK0Cone[4] = {5.,0.3, 0., -1.};
  Double_t xmaxMCrecK0Cone[4] = {100.,0.7, 12., 1.};
  fhnMCrecK0Cone                = new THnSparseF("fhnMCrecK0Cone", "MC rec {K^{0}}^{s} #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecK0Cone,xminMCrecK0Cone,xmaxMCrecK0Cone);  

  Int_t binsMCrecLaCone[4] = {19, 200, 120, 200};
  Double_t xminMCrecLaCone[4] = {5.,0.3, 0., -1.};
  Double_t xmaxMCrecLaCone[4] = {100.,0.7, 12., 1.};
  fhnMCrecLaCone                = new THnSparseF("fhnMCrecLaCone", "MC rec {#Lambda #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecLaCone,xminMCrecLaCone,xmaxMCrecLaCone); 

  Int_t binsMCrecALaCone[4] = {19, 200, 120, 200};
  Double_t xminMCrecALaCone[4] = {5.,0.3, 0., -1.};
  Double_t xmaxMCrecALaCone[4] = {100.,0.7, 12., 1.};
  fhnMCrecALaCone                = new THnSparseF("fhnMCrecALaCone", "MC rec {#bar{#Lambda} #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecALaCone,xminMCrecALaCone,xmaxMCrecALaCone); 
  /*
  Int_t binsMCrecK0ConeSmear[4] = {19, 200, 120, 200};
  Double_t xminMCrecK0ConeSmear[4] = {5.,0.3, 0., -1.};
  Double_t xmaxMCrecK0ConeSmear[4] = {100.,0.7, 12., 1.};
  fhnMCrecK0ConeSmear                = new THnSparseF("fhnMCrecK0ConeSmear", "MC rec {K^{0}}^{s} #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecK0ConeSmear,xminMCrecK0ConeSmear,xmaxMCrecK0ConeSmear);  

  Int_t binsMCrecLaConeSmear[4] = {19, 200, 120, 200};
  Double_t xminMCrecLaConeSmear[4] = {5.,1.05, 0., -1.};
  Double_t xmaxMCrecLaConeSmear[4] = {100.,1.25, 12., 1.};
  fhnMCrecLaConeSmear                = new THnSparseF("fhnMCrecLaConeSmear", "MC rec {#Lambda #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecLaConeSmear,xminMCrecLaConeSmear,xmaxMCrecLaConeSmear); 

  Int_t binsMCrecALaConeSmear[4] = {19, 200, 120, 200};
  Double_t xminMCrecALaConeSmear[4] = {5.,1.05, 0., -1.};
  Double_t xmaxMCrecALaConeSmear[4] = {100.,1.25, 12., 1.};
  fhnMCrecALaConeSmear                = new THnSparseF("fhnMCrecALaConeSmear", "MC rec {#bar{#Lambda} #it{p}_{T}  in cone around jet axis matching MC gen particle; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",4,binsMCrecALaConeSmear,xminMCrecALaConeSmear,xmaxMCrecALaConeSmear); */
  /*
  Int_t binsK0sSecContinCone[3] = {19, 120, 200};
  Double_t xminK0sSecContinCone[3] = {5.,0., -1.};
  Double_t xmaxK0sSecContinCone[3] = {100.,12., 1.};
  fhnK0sSecContinCone                = new THnSparseF("fhnK0sSecContinCone", "Secondary contamination {K^{0}}^{s} #it{p}_{T}  in cone around jet axis; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",3,binsK0sSecContinCone,xminK0sSecContinCone,xmaxK0sSecContinCone);  

  Int_t binsLaSecContinCone[3] = {19, 120, 200};
  Double_t xminLaSecContinCone[3] = {5.,0., -1.};
  Double_t xmaxLaSecContinCone[3] = {100.,12., 1.};
  fhnLaSecContinCone                = new THnSparseF("fhnLaSecContinCone", "Secondary contamination {#Lambda #it{p}_{T}  in cone around jet axis; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",3,binsLaSecContinCone,xminLaSecContinCone,xmaxLaSecContinCone);

  Int_t binsALaSecContinCone[3] = {19, 120, 200};
  Double_t xminALaSecContinCone[3] = {5.,0., -1.};
  Double_t xmaxALaSecContinCone[3] = {100.,12., 1.};
  fhnALaSecContinCone                = new THnSparseF("fhnALaSecContinCone", "Secondary contamination {#bar{#Lambda} #it{p}_{T}  in cone around jet axis; jet #it{p}_{T}; inv mass (GeV/#it{c}^{2};#it{p}_{T}",3,binsALaSecContinCone,xminALaSecContinCone,xmaxALaSecContinCone);
  */
  Int_t binsK0sIncl[3] = {200, 120, 200};
  Double_t xminK0sIncl[3] = {0.3, 0., -1.};
  Double_t xmaxK0sIncl[3] = {0.7, 12., 1.};
  fhnK0sIncl                    = new THnSparseF("fhnK0sIncl","Inclusive K0s",3,binsK0sIncl,xminK0sIncl,xmaxK0sIncl);

  Int_t binsK0sCone[4] = {19, 200, 120, 200};
  Double_t xminK0sCone[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sCone[4] = {100.,0.7, 12., 1.};
  fhnK0sCone                    = new THnSparseF("fhnK0sCone","K0s in jet cone",4,binsK0sCone,xminK0sCone,xmaxK0sCone);

  Int_t binsK0sEmbCone[4] = {19, 200, 120, 200};
  Double_t xminK0sEmbCone[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sEmbCone[4] = {100.,0.7, 12., 1.};
  fhnK0sEmbCone                    = new THnSparseF("fhnK0sEmbCone","Embedded K0s in jet cone",4,binsK0sEmbCone,xminK0sEmbCone,xmaxK0sEmbCone);


  Int_t binsK0sEmbConeRef[4] = {19, 200, 120, 200};
  Double_t xminK0sEmbConeRef[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sEmbConeRef[4] = {100.,0.7, 12., 1.};
  fhnK0sEmbConeRef             = new THnSparseF("fhnK0sEmbConeRef","K0s Embedded reference in jet cone",4,binsK0sEmbConeRef,xminK0sEmbConeRef,xmaxK0sEmbConeRef);

  Int_t binsK0sEmbConeStandard[4] = {19, 200, 120, 200};
  Double_t xminK0sEmbConeStandard[4] = {5.,0.3, 0., -1.};
  Double_t xmaxK0sEmbConeStandard[4] = {100.,0.7, 12., 1.};
  fhnK0sEmbConeStandard             = new THnSparseF("fhnK0sEmbConeStandard","Standard K0s in matched jet cone",4,binsK0sEmbConeStandard,xminK0sEmbConeStandard,xmaxK0sEmbConeStandard);


  Int_t binsLaIncl[3] = {200, 120, 200};
  Double_t xminLaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxLaIncl[3] = {1.25, 12., 1.};
  fhnLaIncl                    = new THnSparseF("fhnLaIncl","Inclusive #Lambda",3,binsLaIncl,xminLaIncl,xmaxLaIncl);

  Int_t binsLaCone[4] = {19, 200, 120, 200};
  Double_t xminLaCone[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaCone[4] = {100.,1.25, 12., 1.};
  fhnLaCone                    = new THnSparseF("fhnLaCone","#Lambda in jet cone",4,binsLaCone,xminLaCone,xmaxLaCone);

  Int_t binsLaEmbCone[4] = {19, 200, 120, 200};
  Double_t xminLaEmbCone[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaEmbCone[4] = {100.,1.25, 12., 1.};
  fhnLaEmbCone                  = new THnSparseF("fhnLaEmbCone","Embedded #Lambda in jet cone",4,binsLaEmbCone,xminLaEmbCone,xmaxLaEmbCone);


  Int_t binsLaEmbConeRef[4] = {19, 200, 120, 200};
  Double_t xminLaEmbConeRef[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaEmbConeRef[4] = {100.,1.25, 12., 1.};
  fhnLaEmbConeRef                  = new THnSparseF("fhnLaEmbConeRef","#Lambda Embedded reference in jet cone",4,binsLaEmbConeRef,xminLaEmbConeRef,xmaxLaEmbConeRef);

  Int_t binsLaEmbConeStandard[4] = {19, 200, 120, 200};
  Double_t xminLaEmbConeStandard[4] = {5.,1.05, 0., -1.};
  Double_t xmaxLaEmbConeStandard[4] = {100.,1.25, 12., 1.};
  fhnLaEmbConeStandard             = new THnSparseF("fhnLaEmbConeStandard","Standard #Lambda in matched jet cone",4,binsLaEmbConeStandard,xminLaEmbConeStandard,xmaxLaEmbConeStandard);


  Int_t binsALaIncl[3] = {200, 120, 200};
  Double_t xminALaIncl[3] = {1.05, 0., -1.};
  Double_t xmaxALaIncl[3] = {1.25, 12., 1.};
  fhnALaIncl                    = new THnSparseF("fhnALaIncl","Inclusive #bar{#Lambda}",3,binsALaIncl,xminALaIncl,xmaxALaIncl);

  Int_t binsALaCone[4] = {19, 200, 120, 200};
  Double_t xminALaCone[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaCone[4] = {100.,1.25, 12., 1.};
  fhnALaCone                    = new THnSparseF("fhnALaCone","#bar{#Lambda} in jet cone",4,binsALaCone,xminALaCone,xmaxALaCone);
 
  Int_t binsALaEmbCone[4] = {19, 200, 120, 200};
  Double_t xminALaEmbCone[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaEmbCone[4] = {100.,1.25, 12., 1.};
  fhnALaEmbCone                 = new THnSparseF("fhnALaEmbCone","Embedded #bar{#Lambda} in jet cone",4,binsALaEmbCone,xminALaEmbCone,xmaxALaEmbCone);

  Int_t binsALaEmbConeRef[4] = {19, 200, 120, 200};
  Double_t xminALaEmbConeRef[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaEmbConeRef[4] = {100.,1.25, 12., 1.};
  fhnALaEmbConeRef                 = new THnSparseF("fhnALaEmbConeRef","#bar{#Lambda} Embedded reference in jet cone",4,binsALaEmbConeRef,xminALaEmbConeRef,xmaxALaEmbConeRef);

  Int_t binsALaEmbConeStandard[4] = {19, 200, 120, 200};
  Double_t xminALaEmbConeStandard[4] = {5.,1.05, 0., -1.};
  Double_t xmaxALaEmbConeStandard[4] = {100.,1.25, 12., 1.};
  fhnALaEmbConeStandard             = new THnSparseF("fhnALaEmbConeStandard","Standard #Lambda in matched jet cone",4,binsALaEmbConeStandard,xminALaEmbConeStandard,xmaxALaEmbConeStandard);


  fh2MCEmbK0sJetPt              = new TH2F("fh2MCEmbK0sJetPt","PYTHIA gen. K^{0}_{s} ;#it{p^{ch,jet}}_{T}; MC gen.#it{p}_{T}",19,5.,100.,120,0.,12.);
  fh2MCEmbLaJetPt               = new TH2F("fh2MCEmbLaJetPt"," PYTHIA gen. #Lambda ;#it{p^{ch,jet}}_{T}; MC gen.#it{p}_{T}",19,5.,100.,120,0.,12.);
  fh2MCEmbALaJetPt              = new TH2F("fh2MCEmbALaJetPt","PYTHIA gen. #bar{#Lambda};#it{p^{ch,jet}}_{T}; MC gen.#it{p}_{T}",19,5.,100.,120,0.,12.);


  fh1MCMultiplicityPrimary      = new TH1F("fh1MCMultiplicityPrimary", "MC Primary Particles;NPrimary;Count", 201, -0.5, 200.5);
  fh1MCMultiplicityTracks       = new TH1F("h1MCMultiplicityTracks", "MC Tracks;Ntracks;Count", 201, -0.5, 200.5);


  Int_t binsFeedDownLa[3] = {19, 200, 120};
  Double_t xminFeedDownLa[3] = {5.,1.05, 0.};
  Double_t xmaxFeedDownLa[3] = {100.,1.25, 12.};
  fhnFeedDownLa                 = new THnSparseF("fhnFeedDownLa","#Lambda stemming from feeddown from Xi(0/-)",3,binsFeedDownLa,xminFeedDownLa,xmaxFeedDownLa);

  Int_t binsFeedDownALa[3] = {19, 200, 120};
  Double_t xminFeedDownALa[3] = {5.,1.05, 0.};
  Double_t xmaxFeedDownALa[3] = {100.,1.25, 12.};
  fhnFeedDownALa                 = new THnSparseF("fhnFeedDownALa","#bar#Lambda stemming from feeddown from Xibar(0/+)",3,binsFeedDownALa,xminFeedDownALa,xmaxFeedDownALa);

  Int_t binsFeedDownLaCone[3] = {19, 200, 120};
  Double_t xminFeedDownLaCone[3] = {5.,1.05, 0.};
  Double_t xmaxFeedDownLaCone[3] = {100.,1.25, 12.};
  fhnFeedDownLaCone             = new THnSparseF("fhnFeedDownLaCone","#Lambda stemming from feeddown from Xi(0/-) in jet cone",3,binsFeedDownLaCone,xminFeedDownLaCone,xmaxFeedDownLaCone);

  Int_t binsFeedDownALaCone[3] = {19, 200, 120};
  Double_t xminFeedDownALaCone[3] = {5.,1.05, 0.};
  Double_t xmaxFeedDownALaCone[3] = {100.,1.25, 12.};
  fhnFeedDownALaCone             = new THnSparseF("fhnFeedDownALaCone","#bar#Lambda stemming from feeddown from Xibar(0/+) in jet cone",3,binsFeedDownALaCone,xminFeedDownALaCone,xmaxFeedDownALaCone);
 
  fh2FeedDownXiLa               = new TH2F("fh2FeedDownXiLa","MC gen. #Xi #it{p}_{T}; MC gen. #Lambda #it{p}_{T}",120,0.,12.,120,0.,12.);

  fh2FeedDownXiALa              = new TH2F("fh2FeedDownXiALa","MC gen. #bar{#Xi} #it{p}_{T}; MC gen.#bar{#Lambda} #it{p}_{T}",120,0.,12.,120,0.,12.);

  fh1MCProdRadiusK0s            = new TH1F("fh1MCProdRadiusK0s","MC gen. MC K0s prod radius",100,0.,100.);
  fh1MCProdRadiusLambda         = new TH1F("fh1MCProdRadiusLambda","MC gen. MC La prod radius",100,0.,100.);
  fh1MCProdRadiusAntiLambda     = new TH1F("fh1MCProdRadiusAntiLambda","MC gen. MC ALa prod radius",100,0.,100.);

  // Pt and inv mass distributions

  fh1MCPtV0s                    = new TH1F("fh1MCPtV0s", "MC gen. V^{0} in rap range;#it{p}_{T} (GeV/#it{c})",120,0,12.);// 0.1 GeV/c steps
  fh1MCPtK0s                    = new TH1F("fh1MCPtK0s", "MC gen. K^{0}_{s} in eta range;#it{p}_{T} (GeV/#it{c})",120,0.,12.);
  fh1MCPtLambda                 = new TH1F("fh1MCPtLambda", "MC gen. #Lambda in rap range;#it{p}_{T} (GeV/#it{c})",120,0.,12.);
  fh1MCPtAntiLambda             = new TH1F("fh1MCPtAntiLambda", "MC gen. #AntiLambda in rap range;#it{p}_{T} (GeV/#it{c})",120,0.,12.);
  fh1MCXiPt                     = new TH1F("fh1MCXiPt", "MC gen. #Xi^{-/o};#it{p}_{T} (GeV/#it{c})",120,0.,12.);
  fh1MCXibarPt                  = new TH1F("fh1MCXibarPt", "MC gen. #bar{#Xi}^{+/o};#it{p}_{T} (GeV/#it{c})",120,0.,12.);
  fh2MCEtaVsPtK0s               = new TH2F("fh2MCEtaVsPtK0s","MC gen. K^{0}_{s} #eta; #it{p}_{T}",120,0.,12.,200,-1.,1.);
  fh2MCEtaVsPtLa                = new TH2F("fh2MCEtaVsPtLa","MC gen. #Lambda #eta; #it{p}_{T}",120,0.,12.,200,-1.,1.);
  fh2MCEtaVsPtALa               = new TH2F("fh2MCEtaVsPtALa","MC gen. #bar{#Lambda}  #eta; #it{p}_{T}",120,0.,12.,200,-1.,1.);

  // Rapidity
  //fh1MCRapK0s                   = new TH1F("fh1MCRapK0s", "MC gen. K0s;rap with cut",200,-10,10); 
  //fh1MCRapLambda                = new TH1F("fh1MCRapLambda", "MC gen. #Lambda;rap",200,-10,10);
  //fh1MCRapAntiLambda            = new TH1F("fh1MCRapAntiLambda", "MC gen. #bar{#Lambda};rap",200,-10,10);
  fh1MCEtaAllK0s                = new TH1F("fh1MCEtaAllK0s", "MC gen. K0s;#eta",200,-1.,1.); 
  fh1MCEtaK0s                   = new TH1F("fh1MCEtaK0s", "MC gen. K0s;#eta with cut",200,-1.,1.); 
  fh1MCEtaLambda                = new TH1F("fh1MCEtaLambda", "MC gen. #Lambda;#eta",200,-1.,1.);
  fh1MCEtaAntiLambda            = new TH1F("fh1MCEtaAntiLambda", "MC gen. #bar{#Lambda};#eta",200,-1.,1.);

  fV0QAK0->DefineHistos();
  fFFHistosRecCuts->DefineHistos();
  fFFHistosRecCutsK0Evt->DefineHistos();
  fFFHistosGen->DefineHistos();

  /* fFFHistosIMK0AllEvt->DefineHistos();
  fFFHistosIMK0Jet->DefineHistos();
  fFFHistosIMK0Cone->DefineHistos();
  fFFHistosIMLaAllEvt->DefineHistos();
  fFFHistosIMLaJet->DefineHistos();
  fFFHistosIMLaCone->DefineHistos();
  fFFHistosIMALaAllEvt->DefineHistos();
  fFFHistosIMALaJet->DefineHistos();
  fFFHistosIMALaCone->DefineHistos();
  */

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
    
    if(fBranchEmbeddedJets.Length()){
    fCommonHistList->Add(fh1nEmbeddedJets);
    fCommonHistList->Add(fh1IndexEmbedded);
    //fCommonHistList->Add(fh1PtEmbExtraOnly);
    fCommonHistList->Add(fh1PtEmbBeforeMatch);
    fCommonHistList->Add(fh1PtEmbReject);
    fCommonHistList->Add(fh2PtEtaEmbReject);
    fCommonHistList->Add(fh1PtEmbAfterMatch);
    fCommonHistList->Add(fh1FractionPtEmbedded);
    fCommonHistList->Add(fh1DeltaREmbedded);
    
    if(fBranchGenJets.Length()&&(fMatchMode == 2)){
    fCommonHistList->Add(fh1IndexEmbeddedMC);
    fCommonHistList->Add(fh1FractionPtEmbeddedMC);
    fCommonHistList->Add(fh2FractionPtVsEmbeddedJetPtMC); 
    fCommonHistList->Add(fh1DeltaREmbeddedMC);
    fCommonHistList->Add(fh1JetPtEmbGenAfterMatch);
    }
    }
    
    fCommonHistList->Add(fh2TracksPerpCone);
    fCommonHistList->Add(fh1PerpCone);
    //fCommonHistList->Add(fh1V0JetPt); 
    fCommonHistList->Add(fh1V0PtCandidate);
    fCommonHistList->Add(fh1IMK0Cone);
    fCommonHistList->Add(fh1IMLaCone);
    fCommonHistList->Add(fh1IMALaCone);
    fCommonHistList->Add(fh1IMK0EmbCone);
    fCommonHistList->Add(fh1IMLaEmbCone);
    fCommonHistList->Add(fh1IMALaEmbCone);
    fCommonHistList->Add(fh2FFJetTrackEta);   
    // fCommonHistList->Add(fh1trackPosNCls);           
    //fCommonHistList->Add(fh1trackNegNCls);          
    fCommonHistList->Add(fh1trackPosEta);            
    fCommonHistList->Add(fh1trackNegEta);          
    fCommonHistList->Add(fh1V0Eta); 
    // fCommonHistList->Add(fh1V0totMom);        
    fCommonHistList->Add(fh1CosPointAngle);                      
    fCommonHistList->Add(fh1DecayLengthV0); 
    fCommonHistList->Add(fh2ProperLifetimeK0sVsPtBeforeCut);
    fCommonHistList->Add(fh2ProperLifetimeK0sVsPtAfterCut);
    fCommonHistList->Add(fh1V0Radius);     
    fCommonHistList->Add(fh1DcaV0Daughters);        
    fCommonHistList->Add(fh1DcaPosToPrimVertex);   
    fCommonHistList->Add(fh1DcaNegToPrimVertex);    
    fCommonHistList->Add(fh2ArmenterosBeforeCuts);
    fCommonHistList->Add(fh2ArmenterosAfterCuts);
    fCommonHistList->Add(fh2BBLaPos);
    fCommonHistList->Add(fh2BBLaNeg);
    fCommonHistList->Add(fh1PosDaughterCharge);
    fCommonHistList->Add(fh1NegDaughterCharge);
    fCommonHistList->Add(fh1PtMCK0s);
    fCommonHistList->Add(fh1PtMCLa);
    fCommonHistList->Add(fh1PtMCALa);
    fCommonHistList->Add(fh1EtaK0s);
    fCommonHistList->Add(fh1EtaLa);
    fCommonHistList->Add(fh1EtaALa);
    fCommonHistList->Add(fh1RC);   
    fCommonHistList->Add(fh1RCBiasK0);
    fCommonHistList->Add(fh1RCBiasLa);
    fCommonHistList->Add(fh1RCBiasALa);                       
    fCommonHistList->Add(fh1MCC);                          
    fCommonHistList->Add(fh1OC);                          
    fCommonHistList->Add(fh1NJ);
    fCommonHistList->Add(fh1NJEmbEvt);
    fCommonHistList->Add(fh1BckgJets);
    fCommonHistList->Add(fh1BckgJetsPtBias);
    fCommonHistList->Add(fhnInvMassEtaTrackPtK0s);
    fCommonHistList->Add(fhnInvMassEtaTrackPtLa);
    fCommonHistList->Add(fhnInvMassEtaTrackPtALa);
    fCommonHistList->Add(fh1TrackMultCone);
    fCommonHistList->Add(fh2TrackMultCone);
    fCommonHistList->Add(fhnNJK0);
    fCommonHistList->Add(fhnNJLa);
    fCommonHistList->Add(fhnNJALa);
    //fCommonHistList->Add(fh2ChTracksNJ);
    //fCommonHistList->Add(fh2ChTracksRC);
    //    fCommonHistList->Add(fh2ChTracksOC);
    //fCommonHistList->Add(fh2ChTracksMCC);
    //fCommonHistList->Add(fh2ChTracksPC);
    fCommonHistList->Add(fhnMCgenK0Cone);
    fCommonHistList->Add(fhnMCgenLaCone);
    //fCommonHistList->Add(fh2MCgenALaCone);
    //fCommonHistList->Add(fh2MCEtagenK0Cone);
    //fCommonHistList->Add(fh2MCEtagenLaCone);
    //fCommonHistList->Add(fh2MCEtagenALaCone);
    /* fCommonHistList->Add(fh2CorrHijingLaProton);
    fCommonHistList->Add(fh2CorrInjectLaProton);
    fCommonHistList->Add(fh2CorrHijingALaAProton);
    fCommonHistList->Add(fh2CorrInjectALaAProton);
    fCommonHistList->Add(fh2MCEtaVsPtHijingLa);
    fCommonHistList->Add(fh2MCEtaVsPtInjectLa);
    fCommonHistList->Add(fh2MCEtaVsPtHijingALa);
    fCommonHistList->Add(fh2MCEtaVsPtInjectALa);
    */
    fCommonHistList->Add(fh1IMK0ConeSmear);
    fCommonHistList->Add(fh1IMLaConeSmear);
    fCommonHistList->Add(fh1IMALaConeSmear);
    fCommonHistList->Add(fh2MC2K0Cone);
    fCommonHistList->Add(fh2MC2LaCone);
    fCommonHistList->Add(fh2MC2ALaCone);
    /* fCommonHistList->Add(fhnrecMCHijingLaIncl);
    fCommonHistList->Add(fhnrecMCHijingLaCone);
    fCommonHistList->Add(fhnrecMCHijingALaIncl);
    fCommonHistList->Add(fhnrecMCHijingALaCone);
    fCommonHistList->Add(fhnrecMCInjectLaIncl);
    fCommonHistList->Add(fhnrecMCInjectLaCone);
    fCommonHistList->Add(fhnrecMCInjectALaIncl);
    fCommonHistList->Add(fhnrecMCInjectALaCone);*/
    /*fCommonHistList->Add(fhnMCrecK0Cone);
    fCommonHistList->Add(fhnMCrecLaCone);
    fCommonHistList->Add(fhnMCrecALaCone); 
    fCommonHistList->Add(fhnMCrecK0ConeSmear);
    fCommonHistList->Add(fhnMCrecLaConeSmear);
    fCommonHistList->Add(fhnMCrecALaConeSmear); */
    /*  fCommonHistList->Add(fhnK0sSecContinCone);
    fCommonHistList->Add(fhnLaSecContinCone);
    fCommonHistList->Add(fhnALaSecContinCone);*/
    fCommonHistList->Add(fhnK0sIncl);
    fCommonHistList->Add(fhnK0sCone);
    fCommonHistList->Add(fhnK0sEmbCone);

    if(fBranchEmbeddedJets.Length()){
    if(fUseExtraTracks)fCommonHistList->Add(fhnK0sEmbConeRef);
    if((fUseExtraTracks == 1) && (fUseStandard == kTRUE)){fCommonHistList->Add(fhnK0sEmbConeStandard);}
    }

    fCommonHistList->Add(fhnLaIncl);
    fCommonHistList->Add(fhnLaCone);
    fCommonHistList->Add(fhnLaEmbCone);

    if(fBranchEmbeddedJets.Length()){
     if(fUseExtraTracks)fCommonHistList->Add(fhnLaEmbConeRef);
     if((fUseExtraTracks == 1) && (fUseStandard == kTRUE)){fCommonHistList->Add(fhnLaEmbConeStandard);}
    }

    fCommonHistList->Add(fhnALaIncl);
    fCommonHistList->Add(fhnALaCone);
    fCommonHistList->Add(fhnALaEmbCone);

    if(fBranchEmbeddedJets.Length()){
      if(fUseExtraTracks)fCommonHistList->Add(fhnALaEmbConeRef);
      if((fUseExtraTracks == 1) && (fUseStandard == kTRUE)){fCommonHistList->Add(fhnALaEmbConeStandard);}
    }
    if((fUseExtraTracks == 1) && (fMatchMode == 2)){
    fCommonHistList->Add(fh2MCEmbK0sJetPt); 
    fCommonHistList->Add(fh2MCEmbLaJetPt); 
    fCommonHistList->Add(fh2MCEmbALaJetPt); 
    }
    fCommonHistList->Add(fhnK0sPC);
    fCommonHistList->Add(fhnK0sEmbPC);
    fCommonHistList->Add(fhnLaPC);
    fCommonHistList->Add(fhnLaEmbPC);
    fCommonHistList->Add(fhnALaPC);
    fCommonHistList->Add(fhnALaEmbPC);
    fCommonHistList->Add(fhnK0sMCC);
    fCommonHistList->Add(fhnLaMCC);
    fCommonHistList->Add(fhnALaMCC);
    fCommonHistList->Add(fhnK0sRC);
    fCommonHistList->Add(fhnLaRC);
    fCommonHistList->Add(fhnALaRC);
    fCommonHistList->Add(fhnK0sRCBias);
    fCommonHistList->Add(fhnLaRCBias);
    fCommonHistList->Add(fhnALaRCBias);
    fCommonHistList->Add(fhnK0sOC);
    fCommonHistList->Add(fhnLaOC);
    fCommonHistList->Add(fhnALaOC);
    fCommonHistList->Add(fh1AreaExcluded); 
    fCommonHistList->Add(fh1MedianEta);
    fCommonHistList->Add(fh1JetPtMedian);
    fCommonHistList->Add(fh1MCMultiplicityPrimary);       
    fCommonHistList->Add(fh1MCMultiplicityTracks);       
    fCommonHistList->Add(fhnFeedDownLa);
    fCommonHistList->Add(fhnFeedDownALa);
    fCommonHistList->Add(fhnFeedDownLaCone);
    fCommonHistList->Add(fhnFeedDownALaCone);
    fCommonHistList->Add(fh2FeedDownXiLa);
    fCommonHistList->Add(fh2FeedDownXiALa);
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
    //fCommonHistList->Add(fh1MCRapK0s);
    //fCommonHistList->Add(fh1MCRapLambda);
    //fCommonHistList->Add(fh1MCRapAntiLambda);   
    fCommonHistList->Add(fh1MCEtaAllK0s);
    fCommonHistList->Add(fh1MCEtaK0s);
    fCommonHistList->Add(fh1MCEtaLambda);
    fCommonHistList->Add(fh1MCEtaAntiLambda);         
    fCommonHistList->Add(fh1nGenJets);

    fV0QAK0->AddToOutput(fCommonHistList);
    fFFHistosRecCuts->AddToOutput(fCommonHistList);
    fFFHistosRecCutsK0Evt->AddToOutput(fCommonHistList);

    if(fBranchGenJets.Length() && (fMatchMode == 2)){
    fFFHistosGen->AddToOutput(fCommonHistList);
    }
    
    // fFFHistosIMK0AllEvt->AddToOutput(fCommonHistList);
    // fFFHistosIMK0Jet->AddToOutput(fCommonHistList);
    // fFFHistosIMK0Cone->AddToOutput(fCommonHistList);
    // fFFHistosIMLaAllEvt->AddToOutput(fCommonHistList);
    // fFFHistosIMLaJet->AddToOutput(fCommonHistList);
    // fFFHistosIMLaCone->AddToOutput(fCommonHistList);
    // fFFHistosIMALaAllEvt->AddToOutput(fCommonHistList);
    // fFFHistosIMALaJet->AddToOutput(fCommonHistList);
    // fFFHistosIMALaCone->AddToOutput(fCommonHistList);
    
 
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
 PostData(1, fCommonHistList); 
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
  
  if((fDebug > 1)&&(!inputHandler)){std::cout<<"AliAnalysisTaskJetChem::AliInputEventHandler does not exist!! "<<std::endl;}

  //for AliPIDResponse:
  fPIDResponse = inputHandler->GetPIDResponse();

  if (!fPIDResponse){if(fDebug > 1) Printf("AliAnalysisTaskJetChem::UserExec(): fPIDResponse does not exist!"); return;}

  if(fDebug > 1){std::cout<<"inputHandler->IsEventSelected(): "<<inputHandler->IsEventSelected()<<std::endl;}
  if(fDebug > 1){std::cout<<"fEvtSelectionMask: "<<fEvtSelectionMask<<std::endl;}
  
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

  // event selection  *****************************************
  //remark: for centrality binning 0-10, 10-30, 30-50, 50-80


  Double_t centPercent = -1;
  Int_t cl = 0;

  if(fEventClass>0){// in case of PbPb, for pp cl = 0
    
    if(handler && handler->InheritsFrom("AliAODInputHandler")){ 
      
      centPercent = dynamic_cast<AliAODHeader*>(fAOD->GetHeader())->GetCentrality();
      cl = 1;
     
      fh1EvtAllCent->Fill(centPercent);
      
      //fEventClass set externally
      if(fEventClass >= 11){//to analyse also 5% central events (PWG-LF choice)

   	if(centPercent < 0) cl = -1;
	if(centPercent >= 0) cl = 11;//exception for analysis of 5% event centrality (PWG-LF choice)
	if(centPercent > 5)  cl = 12; 
	if(centPercent > 10) cl = 13;
	if(centPercent > 20) cl = 14;
	if(centPercent > 40) cl = 15; 
	if(centPercent > 60) cl = 16;
	if(centPercent > 80) cl = 17;
	if(centPercent > 90) cl = 18;

      }
      
      if(fEventClass < 11){//standard centrality estimator used in PWGJE analyses
	
	if(centPercent < 0) cl = -1;
	if(centPercent >= 0) cl = 1;
	if(centPercent > 10) cl = 2; 
	if(centPercent > 30) cl = 3;
	if(centPercent > 50) cl = 4;
	if(centPercent > 80) cl = 5; //takes centralities higher than my upper edge of 80%, not to be used
	
      }
    }//end of AliAOD header

    else {//in case of ESDs in Input

      cl = AliAnalysisHelperJetTasks::EventClass();
      
      if(fESD) centPercent = fESD->GetCentrality()->GetCentralityPercentile("V0M"); //ESD JetServices Task has the centrality binning 0-10,10-30,30-50,50-80
      fh1EvtAllCent->Fill(centPercent);

      
      if(fEventClass >= 11){//to analyse also 5% central events (PWG-LF choice)
	
   	if(centPercent < 0) cl = -1;
	if(centPercent >= 0) cl = 11;//exception for analysis of 5% event centrality (PWG-LF choice)
	if(centPercent > 5)  cl = 12; 
	if(centPercent > 10) cl = 13;
	if(centPercent > 20) cl = 14;
	if(centPercent > 40) cl = 15; 
	if(centPercent > 60) cl = 16;
	if(centPercent > 80) cl = 17;
	if(centPercent > 90) cl = 18;
	
      }
      
      if(fEventClass < 11){//standard centrality estimator used in PWGJE analyses
	
	if(centPercent < 0) cl = -1;
	if(centPercent >= 0) cl = 1;
	if(centPercent > 10) cl = 2; 
	if(centPercent > 30) cl = 3;
	if(centPercent > 50) cl = 4;
	if(centPercent > 80) cl = 5; //takes centralities higher than my upper edge of 80%, not to be used
	
      }
      
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
  
  //std::cout<<"Hallo 1!"<<std::endl;

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

  Int_t nJCuts = GetListOfJets(fJetsRecCuts, kJetsRecAcceptance);//fetch list with jets that survived all jet cuts: fJetsRecCuts, still missing is leading constituent cut

  if(fDebug>2)std::cout<<" nJCuts: "<<nJCuts<<std::endl;

  if(fDebug>2)std::cout<<" fBranchRecJets: "<<fBranchRecJets<<std::endl;
  
  Int_t nRecJetsCuts = 0;                                        //number of reconstructed jets after jet cuts
  if(nJCuts>=0) nRecJetsCuts = fJetsRecCuts->GetEntries(); 
  if(fDebug>2)Printf("%s:%d Selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  if(nRecJetsCuts != nJCuts) Printf("%s:%d Mismatch selected Rec jets after cuts: %d %d",(char*)__FILE__,__LINE__,nJCuts,nRecJetsCuts);
  fh1nRecJetsCuts->Fill(nRecJetsCuts);
  
 
  Int_t nGenJets = 0;
  
  Int_t nEmbeddedJets =  0; 
  Int_t nJGen = 0; 
  TArrayI iEmbeddedMatchIndex; 
  TArrayI iRecMatchIndex;
  TArrayF fRecMatchPtFraction; 
  
  TArrayI iEmbeddedMatchIndexMC; 
  TArrayI iGenMatchIndex; 
  TArrayF fGenMatchPtFraction; 
 
  //fetch all jets used for embedding mode
  if(!(fUseExtraTracks == 0)){ 

    if(fMatchMode == 2){
      nJGen  = GetListOfJets(fJetsGen, kJetsGenAcceptance);//fill list of embedded jets, generator level
      
      if(nJGen>=0) nGenJets = fJetsGen->GetEntries();
      if(fDebug>2)Printf("%s:%d Selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
      
      if(nJGen != nGenJets) Printf("%s:%d Mismatch selected Gen jets: %d %d",(char*)__FILE__,__LINE__,nJGen,nGenJets);
      fh1nGenJets->Fill(nGenJets);
      
    }
  
    if(fDebug>2)std::cout<<"fBranchEmbeddedJets for matching: "<<fBranchEmbeddedJets<<std::endl;


    if(fBranchEmbeddedJets.Length()){ 

    Int_t nJEmbedded = GetListOfJets(fJetsEmbedded, kJetsEmbedded);//fill list of embedded jets, detector level
    if(nJEmbedded>=0) nEmbeddedJets = fJetsEmbedded->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Embedded jets: %d %d",(char*)__FILE__,__LINE__,nJEmbedded,nEmbeddedJets);
    if(nJEmbedded != nEmbeddedJets) Printf("%s:%d Mismatch Selected Embedded Jets: %d %d",(char*)__FILE__,__LINE__,nJEmbedded,nEmbeddedJets);
    fh1nEmbeddedJets->Fill(nEmbeddedJets);
    
    Float_t maxDist = 0.3;//starting value, later the real DeltaR cut will be applied

    iEmbeddedMatchIndex.Set(nEmbeddedJets); 
    iRecMatchIndex.Set(nRecJetsCuts);
    fRecMatchPtFraction.Set(nRecJetsCuts);
    
    iEmbeddedMatchIndex.Reset(-1);
    iRecMatchIndex.Reset(-1);
    fRecMatchPtFraction.Reset(0);
   
    //get closest jets between embedded detector level jets and reconstructed detector level in extra tracks branch

    if(fDebug>2)std::cout<<"GetClosestJets(): "<<std::endl;
    if(fDebug>2)std::cout<<"nRecJetsCuts: "<<nRecJetsCuts<<std::endl;
    if(fDebug>2)std::cout<<"nEmbeddedJets: "<<nEmbeddedJets<<std::endl;


    AliAnalysisHelperJetTasks::GetClosestJets(fJetsEmbedded, nEmbeddedJets, 
					      fJetsRecCuts, nRecJetsCuts, 
                                              iRecMatchIndex,iEmbeddedMatchIndex,
					      0,maxDist);

  
    //std::cout<<"Hallo!"<<std::endl;

    //test:
    //for(Int_t i=0; i<nEmbeddedJets; i++){
    //AliAODJet* embJet = (AliAODJet*) fJetsEmbedded->At(i);
    //if(!embJet) continue; 
      
    //std::cout<<" embJet pt: "<<embJet->Pt()<<" embJet eta: "<<embJet->Eta()<<std::endl;
   
    //}



    // embedded pt fracion
    for(Int_t i=0; i<nRecJetsCuts; i++){
      AliAODJet* recJet = (AliAODJet*) fJetsRecCuts->At(i);
      if(!recJet) continue; 

      Int_t indexEmbedded = iRecMatchIndex[i];
      
      // std::cout<<"indexEmbedded: "<<indexEmbedded<<" recJet pt: "<<recJet->Pt()<<" recJet eta: "<<recJet->Eta()<<std::endl;
      
      //std::cout<<"nEmbeddedJets: "<<nEmbeddedJets<<std::endl;
      
      if(indexEmbedded>-1){
	AliAODJet* embeddedJet = (AliAODJet*) fJetsEmbedded->At(indexEmbedded);
	fRecMatchPtFraction[i] = AliAnalysisHelperJetTasks::GetFractionOfJet(recJet, embeddedJet, 1); // mode 1 / 2nd arg is denominator of fraction

	//std::cout<<"fRecMatchPtFraction["<<i<<"]: "<<fRecMatchPtFraction[i]<<std::endl;


      }
    }
    
    // match embedded rec. level and particle level jets   
    if(nGenJets > 0){

      iEmbeddedMatchIndexMC.Set(nEmbeddedJets); 
      iGenMatchIndex.Set(nGenJets); 
      fGenMatchPtFraction.Set(nGenJets); 
    
      iEmbeddedMatchIndexMC.Reset(-1);
      iGenMatchIndex.Reset(-1);
      fGenMatchPtFraction.Reset(0);
    
      if(fDebug > 2)std::cout<<" nGenJets "<<nGenJets<<" nEmbeddedJets "<<nEmbeddedJets<<std::endl;

      //get closest jets between embedded detector level and particle level PYTHIA jets

      AliAnalysisHelperJetTasks::GetClosestJets(fJetsEmbedded, nEmbeddedJets, 
						fJetsGen, nGenJets, 
						iGenMatchIndex,iEmbeddedMatchIndexMC,
						0,maxDist);

      // embedded pt fraction
      for(Int_t i=0; i<nGenJets; i++){
	AliAODJet* genJet = (AliAODJet*) fJetsGen->At(i);
	if(!genJet) continue; 
	Int_t indexEmbedded = iGenMatchIndex[i];
	if(indexEmbedded>-1){
	  AliAODJet* embeddedJet = (AliAODJet*) fJetsEmbedded->At(indexEmbedded);
	  fGenMatchPtFraction[i] = AliAnalysisHelperJetTasks::GetFractionOfJet(genJet, embeddedJet, 2); // mode 1 / 2nd arg is denominator of fraction
	}
      }
    }
  }

  }//end if embedding mode

  //--------------------------------------
  //____ fetch background clusters ___________(for data analysis these are the kt04 background clusters, for Embedding study these are the standard antikt jets of the data event (in which the PYTHIA jets are embedded into))___________

  if(fBranchRecBckgClusters.Length() != 0){
    
    Int_t nBJ = GetListOfBckgJets(fBckgJetsRec, kJetsRecAcceptance);
    
    if(fDebug>2)std::cout<<" fBranchRecBckgClusters: "<<fBranchRecBckgClusters<<std::endl;
    
    Int_t nRecBckgJets = 0;
    if(nBJ>=0) nRecBckgJets = fBckgJetsRec->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);
    if(nBJ != nRecBckgJets) Printf("%s:%d Mismatch Selected Rec background jets: %d %d",(char*)__FILE__,__LINE__,nBJ,nRecBckgJets);
    
    if(!(fUseExtraTracks == 0)){//for Embedding study only
      fTracksRecBckgCuts->Clear();
      Int_t nTBckgCuts = GetListOfTracks(fTracksRecBckgCuts, kTrackAODCuts);//all standard tracks of data event, no embedded tracks

      if(fDebug>2)Printf("%s:%d selected reconstructed standard Bckg tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTBckgCuts,fTracksRecBckgCuts->GetEntries());

      fIsNJEventEmb = kFALSE;//initialize boolean flag
      
      Int_t nRemainingBckgJets = nRecBckgJets; //init value for NJ events for embedding study 
      for(Int_t ij=0; ij<nRecBckgJets; ++ij){ // loop over all jets in event 
	
	AliAODJet* jetBckg = (AliAODJet*) (fBckgJetsRec->At(ij));
	Double_t sumBckgPt      = 0.;
	Bool_t isBadBckgJet     = kFALSE;
	
	//apply some further jet cuts since they are not carried out in GetListOfBckgJets()
	if( jetBckg->Pt() < fJetPtCut ) isBadBckgJet=kTRUE;
	if( jetBckg->EffectiveAreaCharged() < fJetMinArea ) isBadBckgJet=kTRUE;

	Double_t jetPt = jetBckg->Pt();
	if(isBadBckgJet == kFALSE)fh1BckgJets->Fill(jetPt);//all cuts on jets except LeadingTrackPt cut are applied here

	if((GetFFRadius()<=0)&&(isBadBckgJet == kFALSE)){
	  GetJetTracksTrackrefs(jettracklist, jetBckg, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadBckgJet);// list of jet tracks from trackrefs, normally not used here
	} else {
	  if(isBadBckgJet == kFALSE){GetJetTracksPointing(fTracksRecBckgCuts, jettracklist, jetBckg, GetFFRadius(), sumBckgPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadBckgJet);  // fill list of tracks in cone around jet axis with cone Radius (= 0.4 standard)
	  }
	}

	//APPLICATION OF REMAINING JET CUTS (leading track pt bias etc..):
	if(isBadBckgJet == kTRUE){
	  nRemainingBckgJets = nRemainingBckgJets-1;//remove one jet from nRemainingJets (was initialized with nRecJetsCuts) continue;//all bad jets are rejected
	  
	  if(nRemainingBckgJets == 0){fIsNJEventEmb = kTRUE;fh1NJEmbEvt->Fill(1.);}//set switch for Embedding into NJ events
	}

	if(isBadBckgJet == kFALSE)fh1BckgJetsPtBias->Fill(jetPt);//good jets are filled here
	if(fDebug>3)std::cout<<"isBadBckgJet: "<<isBadBckgJet<<std::endl;

	jettracklist->Clear();
	fTracksRecBckgCuts->Clear();
      }
    }
    
  }
  //------------------------------------
  //____ fetch reconstructed particles __________________________________________________________
  
  Int_t nTCuts;
  
  if(fUseExtraTracks ==  1) nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODExtraCuts);//all tracks in PYTHIA embedded events
    

  else if(fUseExtraTracks == -1) nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODExtraonlyCuts);// only tracks from PYTHIA embedding
  else                           nTCuts = GetListOfTracks(fTracksRecCuts, kTrackAODCuts);//all standard tracks of event, no embedded tracks
  

  if(fDebug>2)Printf("%s:%d selected reconstructed tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());
  if(fTracksRecCuts->GetEntries() != nTCuts) 
    Printf("%s:%d Mismatch selected reconstructed tracks after cuts: %d %d",(char*)__FILE__,__LINE__,nTCuts,fTracksRecCuts->GetEntries());

  if(nTCuts>=0) fh1EvtMult->Fill(fTracksRecCuts->GetEntries());

  Int_t nGenPart = 0;
  Int_t nTGen = 0;

  if(fMatchMode == 2){
    fTracksGen->Clear();
    nTGen = GetListOfTracks(fTracksGen,kTrackAODMCExtraonlyChargedCuts);
    if(nTGen>=0) nGenPart = fTracksGen->GetEntries();
    if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
    if(nGenPart != nTGen) Printf("%s:%d Mismatch selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nTGen,nGenPart);
  }

  //fetch V0 candidates

  Int_t nK0s = 0;
  Int_t nK0sMC = 0;
  Int_t nK0sStandard = 0;  


  if(fUseExtraTracks == 1)    { nK0s = GetListOfV0s(fListK0s,fK0Type,kK0,kTrackAODExtraCuts,myPrimaryVertex,fAOD);//all V0s in event with K0s assumption, all V0 cuts are applied int his function
    if(fUseStandard){nK0sStandard = GetListOfV0s(fListK0sStandard,fK0Type,kK0,kTrackAODCuts,myPrimaryVertex,fAOD);}//fill standard tracks for UE V0 subtraction with PYTHIA Embedding
    if(fMatchMode == 2){nK0sMC= GetListOfV0s(fListK0sMC,fK0Type,kK0,kV0AODMCExtraonlyCuts,myPrimaryVertex,fAOD);}//fill MC K0s in PYTHIA Embedding
  }

  if(fUseExtraTracks == -1)     nK0s = GetListOfV0s(fListK0s,fK0Type,kK0,kTrackAODExtraonlyCuts,myPrimaryVertex,fAOD);// only v0s from PYTHIA embedding
  if(fUseExtraTracks == 0)      nK0s = GetListOfV0s(fListK0s,fK0Type,kK0,kTrackAODCuts,myPrimaryVertex,fAOD);//all standard v0s of event, no embedded tracks
  
  if(fDebug>5){std::cout<<"fK0Type: "<<fK0Type<<" kK0: "<<kK0<<" myPrimaryVertex: "<<myPrimaryVertex<<" fAOD:  "<<fAOD<<std::endl;} 

  /*
  if(fUseExtraTracks == 1) std::cout<< "extra nK0s: "<<nK0s<<std::endl;
  if(fUseExtraTracks == -1) std::cout<< "extraonly nK0s: "<<nK0s<<std::endl;
  if(fUseExtraTracks == 0) std::cout<< "standard nK0s: "<<nK0s<<std::endl;
  */

  if(fDebug>2)Printf("%s:%d Selected Rec K0s candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  if(nK0s != fListK0s->GetEntries()) Printf("%s:%d Mismatch selected K0s: %d %d",(char*)__FILE__,__LINE__,nK0s,fListK0s->GetEntries());
  if(fMatchMode == 2){
  if(fDebug>2)Printf("%s:%d Selected PYTHIA MC gen K0s candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nK0sMC,fListK0sMC->GetEntries());
  }
  fh1K0Mult->Fill(fListK0s->GetEntries());

  
  Int_t nLa = 0;
  Int_t nLaMC = 0;
  Int_t nLaStandard = 0;
  
  if(fUseExtraTracks ==  1)      { nLa = GetListOfV0s(fListLa,fLaType,kLambda,kTrackAODExtraCuts,myPrimaryVertex,fAOD);//all V0s in event with K0s assumption
    if(fUseStandard == kTRUE){nLaStandard = GetListOfV0s(fListLaStandard,fLaType,kLambda,kTrackAODCuts,myPrimaryVertex,fAOD);}
    if(fMatchMode == 2){nLaMC= GetListOfV0s(fListLaMC,fLaType,kLambda,kV0AODMCExtraonlyCuts,myPrimaryVertex,fAOD);}
  }
  
  if(fUseExtraTracks == -1)      nLa = GetListOfV0s(fListLa,fLaType,kLambda,kTrackAODExtraonlyCuts,myPrimaryVertex,fAOD);// only v0s from PYTHIA embedding
  if(fUseExtraTracks ==  0)      nLa = GetListOfV0s(fListLa,fLaType,kLambda,kTrackAODCuts,myPrimaryVertex,fAOD);//all standard tracks of event, no embedded tracks
  
  if(fDebug>2)Printf("%s:%d Selected Rec La candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nLa,fListLa->GetEntries());
  if(nLa != fListLa->GetEntries()) Printf("%s:%d Mismatch selected La: %d %d",(char*)__FILE__,__LINE__,nLa,fListLa->GetEntries());
  fh1LaMult->Fill(fListLa->GetEntries());
  if(fMatchMode == 2){
  if(fDebug>2)Printf("%s:%d Selected PYTHIA MC gen La candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nLaMC,fListLaMC->GetEntries());
  }
  Int_t nALa = 0;
  Int_t nALaMC = 0;
  Int_t nALaStandard = 0;

  if(fUseExtraTracks ==  1)    { nALa = GetListOfV0s(fListALa,fALaType,kAntiLambda,kTrackAODExtraCuts,myPrimaryVertex,fAOD);//all V0s in event with K0s assumption
   if(fUseStandard == kTRUE){nALaStandard = GetListOfV0s(fListALaStandard,fALaType,kAntiLambda,kTrackAODCuts,myPrimaryVertex,fAOD);}
   if(fMatchMode == 2){nALaMC= GetListOfV0s(fListALaMC,fALaType,kAntiLambda,kV0AODMCExtraonlyCuts,myPrimaryVertex,fAOD);}
  }
  if(fUseExtraTracks == -1)    nALa = GetListOfV0s(fListALa,fALaType,kAntiLambda,kTrackAODExtraonlyCuts,myPrimaryVertex,fAOD);// only v0s from PYTHIA embedding
  if(fUseExtraTracks ==  0)    nALa = GetListOfV0s(fListALa,fALaType,kAntiLambda,kTrackAODCuts,myPrimaryVertex,fAOD);//all standard tracks of event, no embedded tracks

  if(fDebug>2)Printf("%s:%d Selected Rec ALa candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nALa,fListALa->GetEntries());
  if(nALa != fListALa->GetEntries()) Printf("%s:%d Mismatch selected ALa: %d %d",(char*)__FILE__,__LINE__,nALa,fListALa->GetEntries());
  if(fMatchMode == 2){
  if(fDebug>2)Printf("%s:%d Selected PYTHIA MC gen ALa candidates after cuts: %d %d",(char*)__FILE__,__LINE__,nALaMC,fListALaMC->GetEntries()); 
  }
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
      
      //Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      
      fh1MCEtaK0s->Fill(fEtaCurrentPart); 
      //fh1MCRapK0s->Fill(fRapCurrentPart);
      fh1MCPtK0s->Fill(fPtCurrentPart);	  
      
      fh2MCEtaVsPtK0s->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered

    }//end of the loop
        
    Int_t nMCgenLa = GetListOfMCParticles(fListMCgenLa,kLambda,fAOD); //fill TList with MC generated primary true Lambdas (list to fill, particletype, mc aod event)
    if(nMCgenLa != fListMCgenLa->GetEntries()) Printf("%s:%d Mismatch selected MCgenLa: %d %d",(char*)__FILE__,__LINE__,nMCgenLa,fListMCgenLa->GetEntries());

    TList *mclist = fAOD->GetList();	
    TClonesArray *stackMC = 0x0;
    stackMC = (TClonesArray*)mclist->FindObject(AliAODMCParticle::StdBranchName());
    if (!stackMC) {
      Printf("ERROR: AliAnalysisTaskJetChem.cxx: loop over MC gen. particles: stackMC not available!");
    }
    
    AliAODMCHeader *mcHdr=(AliAODMCHeader*)mclist->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHdr)Printf("ERROR: AliAnalysisTaskJetChem.cxx: loop over MC gen. particles: mcHdr not available!");
	  
    for(Int_t it=0; it<fListMCgenLa->GetSize(); ++it){ // loop MC generated La, filling histograms
      
      AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLa->At(it));
      if(!mcp0) continue;
	  
      //MC gen Lambdas  
      
      //Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      TString generatorName;
      
      fh1MCEtaLambda->Fill(fEtaCurrentPart); 
      //fh1MCRapLambda->Fill(fRapCurrentPart);
      fh1MCPtLambda->Fill(fPtCurrentPart);	  
      fh2MCEtaVsPtLa->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered

      //Int_t mcp0label = mcp0->GetLabel();
      //      Bool_t istrackInject = IsTrackInjected(mcp0label, mcHdr, stackMC, generatorName);  
    
      //std::cout<<"generatorName: "<<generatorName<<std::endl;

      /*
      if(generatorName == "Hijing"){
	fh2MCEtaVsPtHijingLa->Fill(fPtCurrentPart,fEtaCurrentPart);
      }
      
      if(istrackInject == kTRUE){
	fh2MCEtaVsPtInjectLa->Fill(fPtCurrentPart,fEtaCurrentPart);
	}  */
      
    }//end of the loop
      
      
    Int_t nMCgenALa = GetListOfMCParticles(fListMCgenALa,kAntiLambda,fAOD); //fill TList with MC generated primary true Antilambdas (list to fill, particletype, mc aod event)
    if(nMCgenALa != fListMCgenALa->GetEntries()) Printf("%s:%d Mismatch selected MCgenALa: %d %d",(char*)__FILE__,__LINE__,nMCgenALa,fListMCgenALa->GetEntries());
  
   	
    for(Int_t it=0; it<fListMCgenALa->GetSize(); ++it){ // loop MC generated ALa, filling histograms
      
      AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenALa->At(it));
      if(!mcp0) continue;
      
      //MC gen Antilambdas                  
      
      // Double_t fRapCurrentPart   = MyRapidity(mcp0->E(),mcp0->Pz());
      Double_t fEtaCurrentPart   = mcp0->Eta();
      Double_t fPtCurrentPart    = mcp0->Pt();
      TString generatorName;

      fh1MCEtaAntiLambda->Fill(fEtaCurrentPart); 
      //fh1MCRapAntiLambda->Fill(fRapCurrentPart);
      fh1MCPtAntiLambda->Fill(fPtCurrentPart);	  
      fh2MCEtaVsPtALa->Fill(fPtCurrentPart,fEtaCurrentPart);                  //eta cut, physical primary selection and decay mode considered


      // Int_t mcp0label = mcp0->GetLabel();
      //      Bool_t istrackInject = IsTrackInjected(mcp0label, mcHdr, stackMC, generatorName);  
    
      //std::cout<<"generatorName: "<<generatorName<<std::endl;

      /*
      if(generatorName == "Hijing"){
	fh2MCEtaVsPtHijingALa->Fill(fPtCurrentPart,fEtaCurrentPart);
      }
      
      if(istrackInject == kTRUE){
	fh2MCEtaVsPtInjectALa->Fill(fPtCurrentPart,fEtaCurrentPart);
	}  */

	
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
  Double_t dRadiusExcludeCone = 2*GetFFRadius(); //2 times jet radius
 
  //##########################################
 //------------------------------------------
  for(Int_t it=0; it<fListK0s->GetSize(); ++it){ //loop over all K0s candidates in PbPb event
        
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
    if(!v0) continue;

    // VO's main characteristics to check the reconstruction cuts
    
    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
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
    
 
    //OUTSIDE CONES:########

    //V0 particles
 
    Double_t fEta = v0->PseudoRapV0();
    Bool_t bIsInCone = kFALSE;//init boolean, is not in any cone (OC)
    Int_t nRemainingJets = nRecJetsCuts; //init value    Int_t nRemainingJets = nRecJetsCuts; //init value 

    for(Int_t ij=0; ij<nRecJetsCuts; ++ij){ // loop over all jets in event 
      
      AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));
      jettracklist->Clear();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
 
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);// list of jet tracks from trackrefs
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);  // fill list of tracks in cone around jet axis with cone Radius (= 0.4 standard)
      }
      
      //leading track pt bias on jets inside this small jet loop
      //APPLICATION OF REMAINING JET CUTS (leading track pt bias etc..):
      if(isBadJet){
	nRemainingJets = nRemainingJets-1;//remove one jet from nRemainingJets (was initialized with nRecJetsCuts) //all bad jets are rejected
	continue;
      }

      //if jet is selected, then check whether V0 is part of the jet cone:
      if(IsParticleInCone(jet, v0, dRadiusExcludeCone) == kTRUE) {bIsInCone = kTRUE;}//is v0 of PbPb event part of the (selected) jet cone?
      
      jettracklist->Clear();
    }
    
    //################OCs##########################################standard method#####################################
    if(fUseExtraTracks == 0){//standard PbPb analysis
      if((bIsInCone==kFALSE)&&(nRemainingJets > 0)){//K0s is not part of any selected jet in event, but its a jet event
	Double_t vK0sOC[3] = {invMK0s,trackPt,fEta};
	fhnK0sOC->Fill(vK0sOC);      
      }
    }
    
    

     //#################################################################################################################
   
    //end of outside cone K0s
   
    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);

    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();
    
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);

    
    fV0QAK0->FillTrackQA(v0->Eta(), TVector2::Phi_0_2pi(v0->Phi()), v0->Pt()); 
    //fFFHistosIMK0AllEvt->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaK0s->Fill(fEta);

    Double_t vK0sIncl[3] = {invMK0s,trackPt,fEta}; //fill all K0s in event into THnSparse of 3 dimensions
    fhnK0sIncl->Fill(vK0sIncl);


    if(fAnalysisMC){
      TString generatorName;
      Bool_t isinjected;
      TList *listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kK0, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
      //if(fPhysicalPrimary == kFALSE)continue;
      //std::cout<<"mclabelcheck: "<<mclabelcheck<<std::endl;
      //std::cout<<"IsPhysicalPrimary: "<<fPhysicalPrimary<<std::endl;

      if(mclabelcheck == kFALSE)continue;
      
      Double_t vInvMassEtaTrackPtK0s[3] = {fEta,invMK0s,trackPt};
      fhnInvMassEtaTrackPtK0s->Fill(vInvMassEtaTrackPtK0s);//includes also feeddown particles, mainly phi particles whose decay products are considered here as primary


      fh1PtMCK0s->Fill(MCPt);
    }
 

    fh1V0Eta->Fill(fEta);
    //fh1V0totMom->Fill(fV0TotalMomentum);
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
    //  Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
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
    
    Double_t fEta = v0->PseudoRapV0();
    Bool_t bIsInCone = kFALSE;//init boolean, is not in any cone (OC)
    Int_t nRemainingJets = nRecJetsCuts; //init value 

    CalculateInvMass(v0, kLambda, invMLa, trackPt);//function to calculate invMass with TLorentzVector class
    
    
    for(Int_t ij=0; ij<nRecJetsCuts; ++ij){ // loop over all jets in event 
      
      AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));
      jettracklist->Clear();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
 
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);// list of jet tracks from trackrefs
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);  // fill list of tracks in cone around jet axis with cone Radius (= 0.4 standard) + Check of leading constituent condition (= minPtL) carried out by method implemented in class AliAnalysisTaskFragmentationFunction
      }
      

      //leading track pt bias on jets inside this small jet loop
   
      if(isBadJet){
	nRemainingJets = nRemainingJets-1;//remove one jet from nRemainingJets (was initialized with nRecJetsCuts) continue;//all bad jets are rejected
	continue;
      }

      if(IsParticleInCone(jet, v0, dRadiusExcludeCone) == kTRUE) {bIsInCone = kTRUE;}
     
      jettracklist->Clear();  
    }  //end jet loop  
    
//at this point, the v0 candidate received a flag (bIsInCone) whether it belongs to any selected jet in event

 //standard OCs method:

    if(fUseExtraTracks == 0){//standard tracks
      
      if((bIsInCone == kFALSE)&&(nRemainingJets > 0)){//success! Lambda doesn't belong to any selected jet in event
	Double_t vLaOC[3] = {invMLa, trackPt,fEta};
	fhnLaOC->Fill(vLaOC); 
      }
    } 


    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();  
    
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
    
    //fFFHistosIMLaAllEvt->FillFF(trackPt, invMLa, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaLa->Fill(fEta);

    Double_t vLaIncl[3] = {invMLa,trackPt,fEta};
    fhnLaIncl->Fill(vLaIncl);

    if(fAnalysisMC){  
   
      TString generatorName;
      Bool_t isinjected;
      TList* listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
      if(mclabelcheck == kFALSE)continue; 
      //if(fPhysicalPrimary == kFALSE)continue;
      /*
      if(generatorName == "Hijing"){
	Double_t vrecMCHijingLaIncl[3] = {invMLa,trackPt,fEta};
	fhnrecMCHijingLaIncl->Fill(vrecMCHijingLaIncl);

	Double_t protonPt = trackPos->Pt();
      	fh2CorrHijingLaProton->Fill(trackPt,protonPt);
      }

      if(isinjected == kTRUE){
	Double_t vrecMCInjectLaIncl[3] = {invMLa,trackPt,fEta};
	fhnrecMCInjectLaIncl->Fill(vrecMCInjectLaIncl);

	Double_t protonPt = trackPos->Pt();
      	fh2CorrInjectLaProton->Fill(trackPt,protonPt);
	}*/

      Double_t vInvMassEtaTrackPtLa[3] = {fEta,invMLa,trackPt};
      fhnInvMassEtaTrackPtLa->Fill(vInvMassEtaTrackPtLa);//includes also feed-down particles
      fh1PtMCLa->Fill(MCPt);
      

      fh1PtMCLa->Fill(MCPt);
  

    }
    fh1V0Eta->Fill(fEta);
    //fh1V0totMom->Fill(fV0TotalMomentum);
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

    Double_t fEta = v0->PseudoRapV0();
    Bool_t bIsInCone = kFALSE;//init boolean for OC     
    Int_t nRemainingJets = nRecJetsCuts; //init value 
       

    CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
      
    for(Int_t ij=0; ij<nRecJetsCuts; ++ij){ // loop over all jets in event 
      
      AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));
      jettracklist->Clear();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
 

      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);// list of jet tracks from trackrefs
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);  // fill list of tracks in cone around jet axis with cone Radius (= 0.4 standard)
      }
      
      //leading track pt bias on jets inside this small jet loop
  
      if(isBadJet){
	nRemainingJets = nRemainingJets-1;//remove one jet from nRemainingJets (was initialized with nRecJetsCuts) continue;//all bad jets are rejected
	continue;
      }


      if(IsParticleInCone(jet, v0, dRadiusExcludeCone) == kTRUE){
	bIsInCone = kTRUE;	
      }

      jettracklist->Clear();
    }
 

    if(fUseExtraTracks == 0){//standard PbPb analysis
      if((bIsInCone == kFALSE)&&(nRemainingJets > 0)){//success!
	Double_t vALaOC[3] = {invMALa, trackPt,fEta};
	fhnALaOC->Fill(vALaOC); 
      }
    }

    
    Double_t fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
    lV0Position[0]= v0->DecayVertexV0X();  
    lV0Position[1]= v0->DecayVertexV0Y();  
    lV0Position[2]= v0->DecayVertexV0Z();  
    Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
    fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
       
    //fFFHistosIMALaAllEvt->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
    //fh1trackPosNCls->Fill(trackPosNcls);
    //fh1trackNegNCls->Fill(trackNegNcls);
    fh1EtaALa->Fill(fEta);

    Double_t vALaIncl[3] = {invMALa,trackPt,fEta};
    fhnALaIncl->Fill(vALaIncl);

    if(fAnalysisMC){
      TString generatorName;
      Bool_t isinjected;
      TList* listmc = fAOD->GetList();
      Bool_t mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
      if(mclabelcheck == kFALSE)continue; 
      //if(fPhysicalPrimary == kFALSE)continue;//take also feeddown particles into account
      /*
      if(generatorName == "Hijing"){
	Double_t vrecMCHijingALaIncl[3] = {invMALa,trackPt,fEta};
	fhnrecMCHijingALaIncl->Fill(vrecMCHijingALaIncl);

	Double_t aprotonPt = trackNeg->Pt();
       	fh2CorrHijingALaAProton->Fill(trackPt,aprotonPt);
      }


      if(isinjected == kTRUE){
	Double_t vrecMCInjectALaIncl[3] = {invMALa,trackPt,fEta};
	fhnrecMCInjectALaIncl->Fill(vrecMCInjectALaIncl);

	Double_t aprotonPt = trackNeg->Pt();
       	fh2CorrInjectALaAProton->Fill(trackPt,aprotonPt);

	}*/


      Double_t vInvMassEtaTrackPtALa[3] = {fEta,invMALa,trackPt};
      fhnInvMassEtaTrackPtALa->Fill(vInvMassEtaTrackPtALa);
      fh1PtMCALa->Fill(MCPt);

    }
    fh1V0Eta->Fill(fEta);
    //fh1V0totMom->Fill(fV0TotalMomentum);
    fh1CosPointAngle->Fill(fV0cosPointAngle);
    fh1DecayLengthV0->Fill(fV0DecayLength);
    fh1V0Radius->Fill(fV0Radius);
    fh1DcaV0Daughters->Fill(fDcaV0Daughters);
    fh1DcaPosToPrimVertex->Fill(fDcaPosToPrimVertex);
    fh1DcaNegToPrimVertex->Fill(fDcaNegToPrimVertex);
    fh1trackPosEta->Fill(PosEta);
    fh1trackNegEta->Fill(NegEta);
  }
  
  //_____no jets events______________________________________________________________________________________________________________________________________

  if(nRecJetsCuts == 0){//no jet events, before the remaining jet cuts are applied, the second part for the non-jet events comes inside the jet loop
        
    fh1NJ->Fill(1.);//for normalisation by number of NJ events for events in which no rec. jets are found right from the beginning and before even the leading track bias is applied
    
    if(fDebug>6) { std::cout<<"################## nRecJetsCuts == 0 ###################"<<std::endl;
      //std::cout<<"fListK0s->GetSize() in NJ event: "<<fListK0s->GetSize()<<std::endl;
    }
    
    for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 
      
      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
      if(!v0) continue;
      
      Double_t invMK0s =0;
      Double_t trackPt=0;
      CalculateInvMass(v0, kK0, invMK0s, trackPt);
      Double_t fEta = v0->Eta();

      Double_t vNJK0[3] = {invMK0s,trackPt,fEta}; //fill all K0s in events wo selected jets
      fhnNJK0->Fill(vNJK0);
      
    }
    
    for(Int_t it=0; it<fListLa->GetSize(); ++it){ // loop all La 
      
      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLa->At(it));
      if(!v0) continue;
      
      Double_t invMLa =0;
      Double_t trackPt=0;	
      CalculateInvMass(v0, kLambda, invMLa, trackPt);
      Double_t fEta = v0->Eta();

      Double_t vNJLa[3] = {invMLa,trackPt,fEta}; //fill all K0s in events wo selected jets
      fhnNJLa->Fill(vNJLa);

    } 
    
    for(Int_t it=0; it<fListALa->GetSize(); ++it){ // loop all ALa 
      
      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALa->At(it));
      if(!v0) continue;
      
      Double_t invMALa =0;
      Double_t trackPt=0;	
      CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);

      Double_t fEta = v0->Eta();

      Double_t vNJALa[3] = {invMALa,trackPt,fEta}; //fill all K0s in events wo selected jets
      fhnNJALa->Fill(vNJALa);
   
      
    } 
    
  }//no jet events
  
  //____ fill all jet related histos  ________________________________________________________________________________________________________________________
  //##########################################################################################################################################################
  //##########################JET LOOP########################################################################################################################
  //##########################################################################################################################################################
  //##########################################################################################################################################################

  Int_t nSelJets = nRecJetsCuts; //init value
  Bool_t IsOCEvt = kFALSE; //init for this outside cones normalisation histo (total number of OC events)
  Bool_t IsRCEvt = kFALSE; //init for that the random cone is placed only once per event
  Bool_t IsMCCEvt = kFALSE; //init for that the median cluster cone is placed only once per event

  //fill jet histos in general
  for(Int_t ij=0; ij<nRecJetsCuts; ++ij){                               // ij is an index running over the list of the reconstructed jets after most of the cuts, but not yet the leading track bias, all jets in event are looped
    
    AliAODJet* jet = (AliAODJet*) (fJetsRecCuts->At(ij));
    
    Double_t jetPt  = jet->Pt();
    Double_t jetEta = jet->Eta();
    Double_t jetPhi = jet->Phi();
    
    //if(ij==0){ // loop over leading jets for ij = 0, for ij>= 0 look into all jets
    
    if(ij>=0){//all jets in event
      
      jettracklist->Clear();
      Double_t sumPt      = 0.;
      Bool_t isBadJet     = kFALSE;
      Int_t njetTracks    = 0;
      
      if(GetFFRadius()<=0){
 	GetJetTracksTrackrefs(jettracklist, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);// list of jet tracks from trackrefs
      } else {
 	GetJetTracksPointing(fTracksRecCuts, jettracklist, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);  // fill list of charged hybrid tracks in cone around jet axis with cone Radius (= 0.4 standard), check of leading track bias condition
      }
      
      //not applied at the moment:
      if((GetFFMinNTracks()>0) && (jettracklist->GetSize() <= GetFFMinNTracks())) isBadJet = kTRUE; // reject jets with less tracks than fFFMinNTracks
      
      //APPLICATION OF REMAINING JET CUTS (leading track pt bias etc..) + NJ events
      
      if(isBadJet) {
	
	nSelJets = nSelJets-1;//remove one jet from nSelJets (was initialized with nRecJetsCuts)
	
	if(nSelJets == 0){//case that event doesn't contain no selected jets at all and there are no jets remaining to be looped over
	  
	  fh1NJ->Fill(1.);//for normalisation by number of NJ events
     	  
	  for(Int_t it=0; it<fListK0s->GetSize(); ++it){ // loop all K0s 
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0s->At(it));
	    if(!v0) continue;
	    
	    Double_t invMK0s =0;
	    Double_t trackPt=0;
	    CalculateInvMass(v0, kK0, invMK0s, trackPt);
	    Double_t fEta = v0->Eta();
	    
	    Double_t vNJK0[3] = {invMK0s,trackPt,fEta}; //fill all K0s in events wo selected jets
	    fhnNJK0->Fill(vNJK0);
	    
	  }
	  
	  for(Int_t it=0; it<fListLa->GetSize(); ++it){ // loop all La 
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLa->At(it));
	    if(!v0) continue;
	    
	    Double_t invMLa =0;
	    Double_t trackPt=0;	
	    CalculateInvMass(v0, kLambda, invMLa, trackPt);
	    Double_t fEta = v0->Eta();
	    
	    Double_t vNJLa[3] = {invMLa,trackPt,fEta}; //fill all K0s in events wo selected jets
	    fhnNJLa->Fill(vNJLa);
	    
	   }
	  
	  for(Int_t it=0; it<fListALa->GetSize(); ++it){ // loop all ALa 
      
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALa->At(it));
	    if(!v0) continue;
	    
	    Double_t invMALa =0;
	    Double_t trackPt=0;	
	    CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);
	    
	    Double_t fEta = v0->Eta();
	    
	    Double_t vNJALa[3] = {invMALa,trackPt,fEta}; //fill all K0s in events wo selected jets
	    fhnNJALa->Fill(vNJALa);
	    
	    
	  } 
	  
	  //Analysing charged tracks in non-jet events: for estimation of syst. uncertainty of UE subtraction methods
	  /*
	  for(Int_t in=0; in<fTracksRecCuts->GetEntries(); ++in){
	    
	    AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksRecCuts->At(in)); //inputlist is fListK0s, all reconstructed K0s in event
	    if(!track){std::cout<<"AliAnalysisTaskJetChem(): In NJ event, charged track not found!!!"<<std::endl; continue;}
	    
	    Double_t trackPt = track->Pt();
	    Double_t trackEta = track->Eta();
	    if(TMath::Abs(trackEta) <= fCutEta){
	    fh2ChTracksNJ->Fill(trackPt,trackEta);   	    
	    }
	  }*/
	  
	}


	continue;//rejection of current jet
      } // rejects jets in which no track has a track pt higher than 5 GeV/c (see AddTask macro)
      
      if(IsOCEvt == kFALSE){IsOCEvt = kTRUE;fh1OC->Fill(1.);}//the first found jet triggers an OC event and is filled only once into normalisation histo      
      //Float_t fJetAreaMin = 0.6*TMath::Pi()*GetFFRadius()*GetFFRadius(); // minimum jet area cut, already applied in JetListOfJets() in FF Task

      //if(fDebug > 2)  {if (jet->EffectiveAreaCharged() < fJetAreaMin) {std::cout<<" fCutjetArea cut removed a jet!!!!! Should not have to be done again!!"<<std::endl;}}// cut on jet area, already done by jet selection in FF task
      
      Double_t dAreaExcluded = TMath::Pi()*dRadiusExcludeCone*dRadiusExcludeCone; // area of the cone
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone,fCutjetEta-jet->Eta()); // positive eta overhang
      dAreaExcluded -= AreaCircSegment(dRadiusExcludeCone,fCutjetEta+jet->Eta()); // negative eta overhang
      fh1AreaExcluded->Fill(dAreaExcluded);//histo contains all areas that are jet related and have to be excluded concerning OC UE pt spectrum normalisation by area

      fh1JetEta->Fill(jetEta);        
      fh1JetPhi->Fill(jetPhi);                
      fh2JetEtaPhi->Fill(jetEta,jetPhi);  
  
      // printf("pT = %f, eta = %f, phi = %f, leadtr pt = %f\n, ",jetPt,jetEta,jetphi,leadtrack);

    
      //###################PYTHIA JET EMBEDDING#####################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
     


      if(!(fUseExtraTracks == 0)&&(fIsNJEventEmb == kTRUE)){//this following big block is used only for Embedding study
    
	
	TList* mclist = fAOD->GetList();
	if (!mclist){std::cout<<"mclist does not exist for Embedding study: "<<std::endl; continue;}
	
	Double_t ptFractionEmbedded = 0; 
	Double_t deltaREmbedded     = 0;
	AliAODJet* embeddedJet      = 0; 

	if(fDebug>2)std::cout<<"fBranchEmbeddedJets before index embedded: "<<fBranchEmbeddedJets<<std::endl;
	
	if(fBranchEmbeddedJets.Length()){ // find embedded jet
	  
	  Int_t indexEmbedded         = iRecMatchIndex[ij];
	  ptFractionEmbedded = fRecMatchPtFraction[ij]; 
	  
	  if(fDebug>2)std::cout<<"index embedded: "<<indexEmbedded<<std::endl;
	  if(fDebug>2)std::cout<<"ptFractionEmbedded: "<<ptFractionEmbedded<<std::endl;

  
	  fh1IndexEmbedded->Fill(indexEmbedded);
	  fh1FractionPtEmbedded->Fill(ptFractionEmbedded);
	  
	  if(indexEmbedded>-1){ 
	    
	    embeddedJet = dynamic_cast<AliAODJet*>(fJetsEmbedded->At(indexEmbedded));//fetch embedded jet
	    if(!embeddedJet) continue;

	    //std::cout<<"pointer to embeddedJet: "<<embeddedJet<<std::endl;	
    
	    deltaREmbedded   = jet->DeltaR((AliVParticle*) (embeddedJet)); 
	    if(fDebug>2)std::cout<<"deltaREmbedded: "<<deltaREmbedded<<std::endl;	  
  
	    fh1DeltaREmbedded->Fill(deltaREmbedded);

	  }
	}
	
	
	if(!embeddedJet)continue;
	
	
	Double_t JetPtEmb = embeddedJet->Pt();
	fh1PtEmbBeforeMatch->Fill(JetPtEmb);
	
	
	if((ptFractionEmbedded < fCutFractionPtEmbedded) || (deltaREmbedded > fCutDeltaREmbedded)){
	  Double_t JetPtRej = embeddedJet->Pt();
	  Double_t JetEtaRej = embeddedJet->Eta();

	  fh1PtEmbReject->Fill(JetPtRej);
	  fh2PtEtaEmbReject->Fill(JetPtRej,JetEtaRej);
	}
	
	
	if(ptFractionEmbedded >= fCutFractionPtEmbedded && deltaREmbedded <= fCutDeltaREmbedded) // end: cut embedded ratio
	  {
	    if(fMatchMode == 1){
	      FillEmbeddedHistos(embeddedJet, jet, nK0s, nLa, nALa, jettracklist);//fetch V0s for matched jets and fill embedding histos, 'jet' is matched jet here
	    }
	  }
	//################################end V0 embedding part
	//################################
	
      }//end of fTracksExtra != 0 check, end of embedding part
      
      
	//#####################End of embedding study in MatchMode 1################################################################################


      //############################################################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
      //############################################################################################################################################
    
	//std::cout<<"fUseExtraTracks: "<<fUseExtraTracks<<std::endl;
	
	if(fUseExtraTracks == 0){//no embedded jets are used, normally the case
	  
	
	  for(Int_t it=0; it<jettracklist->GetSize(); ++it){//loop over all charged tracks in jet
	    
	    AliVParticle* trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));//all tracks in jet cone	
	    if(!trackVP)continue;
	  
	    Float_t trackPt = trackVP->Pt();//transversal momentum of jet particle
	    Float_t trackEta = trackVP->Eta();
	    
	    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    
	    fFFHistosRecCuts->FillFF(trackPt, jetPt, incrementJetPt);//histo with tracks/jets after cut selection, for all events
	    if(nK0s>0) fFFHistosRecCutsK0Evt->FillFF(trackPt, jetPt, incrementJetPt);//only for K0s events
	    fh2FFJetTrackEta->Fill(trackEta,jetPt);
	    
	  }
	}
	

	njetTracks = jettracklist->GetSize();
	
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
	//	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	
	//	if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	Double_t invMK0s =0;
	Double_t trackPt=0;	
	CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	
	//	fFFHistosIMK0Jet->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
 	

	if(dPhiJetK0<fh1dPhiJetK0->GetXaxis()->GetXmin()) dPhiJetK0 += 2*TMath::Pi();
	fh1dPhiJetK0->Fill(dPhiJetK0);
	
      }

      // if(fListK0s->GetSize() == 0){ // no K0: increment jet pt spectrum 
	
      //	Bool_t incrementJetPt = kTRUE;
	//	fFFHistosIMK0Jet->FillFF(-1, -1, jetPt, incrementJetPt);
      //  }
      
      //____fetch reconstructed K0s in cone around jet axis:_______________________________________________________________________________
      
      jetConeK0list->Clear();

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
	Double_t fEta=0;
	fEta = v0->Eta();
	
	CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class


	if(fAnalysisMC){
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,jetPtSmear);

	  fh2MC2K0Cone->Fill(jetPtSmear, trackPt); //fill MC truth for jet pt smearing

	  if(incrementJetPt == kTRUE){fh1IMK0ConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 

	  }

	if(incrementJetPt==kTRUE){
	  fh1IMK0Cone->Fill(jetPt);}//normalisation by number of selected jets

	//fFFHistosIMK0Cone->FillFF(trackPt, invMK0s, jetPt, incrementJetPt);
	
	Double_t vK0sCone[4] = {jetPt, invMK0s,trackPt,fEta};
	fhnK0sCone->Fill(vK0sCone);
      }
      
      
      if(jetConeK0list->GetSize() == 0){ // no K0: increment jet pt spectrum 
	
	
	Bool_t incrementJetPt = kTRUE;//jets without K0s will be only filled in TH1F only once, so no increment needed 
	//fFFHistosIMK0Cone->FillFF(-1, -1, jetPt, incrementJetPt);
	Double_t vK0sCone[4] = {jetPt, -1, -1, -1};
	fhnK0sCone->Fill(vK0sCone);

	if(incrementJetPt==kTRUE){
	  fh1IMK0Cone->Fill(jetPt);}//normalisation by number of selected jets

	if(fAnalysisMC){
	  Double_t jetPtSmear = -1;  
	  SmearJetPt(jetPt,jetPtSmear);  
	  if(incrementJetPt == kTRUE){fh1IMK0ConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	  }
      }    
      
      //Random cones________________________________________________________________________
     

      if(IsRCEvt == kFALSE){//fetch random cone V0s only once per event
	

	IsRCEvt = kTRUE;//set boolean to kTRUE once a random cone is placed per event
	
	AliAODJet* jetRC = 0;
	jetRC = GetRandomCone(fJetsRecCuts, fCutjetEta, 2*GetFFRadius());//fetch one random cone for each event 
	
	fListK0sRC->Clear();//list for K0s in random cone (RC), one RC per event
	fListLaRC->Clear();
	fListALaRC->Clear();
	fTracksRecCutsRC->Clear();


	Double_t sumPtK0sRC = 0;
	Double_t sumPtLaRC = 0;
	Double_t sumPtALaRC = 0;
	Double_t sumPtRecCutsRC = 0;
	Bool_t isBadJetK0sRC = kFALSE;
	Bool_t isBadJetLaRC = kFALSE;
	Bool_t isBadJetALaRC = kFALSE;
	Bool_t isBadJetRecCutsRC = kFALSE;
	
	if(jetRC != 0) {//if random cone was selected properly and fullfilling all the requirements

	//fetch V0s in RC:
	  fh1RC->Fill(1.);//for normalisation purposes

	  GetTracksInCone(fListK0s, fListK0sRC, jetRC, GetFFRadius(), sumPtK0sRC, 0, 0, isBadJetK0sRC);
	  
	  //________________fill RC with all V0s__________________
	  for(Int_t it=0; it<fListK0sRC->GetSize(); ++it){ // loop for K0s in random cone
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0sRC->At(it));
	    if(!v0) continue;
	    
	    Double_t invMK0s =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
	    fEta = v0->Eta();
		
	    CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vK0sRC[3] = {invMK0s,trackPt,fEta};
	    fhnK0sRC->Fill(vK0sRC);
	  }
	  
	
	  
	  GetTracksInCone(fListLa, fListLaRC, jetRC, GetFFRadius(), sumPtLaRC, 0, 0, isBadJetLaRC);
	  
	  for(Int_t it=0; it<fListLaRC->GetSize(); ++it){ // loop for Lambdas in random cone
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLaRC->At(it));
	    if(!v0) continue;
	    
	    Double_t invMLa =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
	    fEta = v0->Eta();
	    
	    CalculateInvMass(v0, kLambda, invMLa, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vLaRC[3] = {invMLa,trackPt,fEta};
	    fhnLaRC->Fill(vLaRC);
	  }
	
	
	  
	  GetTracksInCone(fListALa, fListALaRC, jetRC, GetFFRadius(), sumPtALaRC, 0, 0, isBadJetALaRC);
	  
	  for(Int_t it=0; it<fListALaRC->GetSize(); ++it){ // loop for Lambdas in random cone
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALaRC->At(it));
	    if(!v0) continue;
	    
	    Double_t invMALa =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
	    fEta = v0->Eta();
		
	    CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vALaRC[3] = {invMALa,trackPt,fEta};
	    fhnALaRC->Fill(vALaRC);
	  }


	  if(isBadJetK0sRC == kFALSE){ //in case RC contains at least one K0s with minimum pT 
	    fh1RCBiasK0->Fill(1.);//for normalisation purposes

	    //________________fill RC (with trigger particle bias)_____________
	    for(Int_t it=0; it<fListK0sRC->GetSize(); ++it){ // loop for K0s in random cone
	      
	      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0sRC->At(it));
	      if(!v0) continue;
	      
	      Double_t invMK0s =0;
	      Double_t trackPt=0;
	      Double_t fEta=0;
	      fEta = v0->Eta();
	      
	      CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	      
	      //Double_t vK0sRC[3] = {invMK0s,trackPt,fEta};
	      //fhnK0sRCBias->Fill(vK0sRC);
	    }
	  }
	
	  
	  if(isBadJetLaRC == kFALSE){ //in case RC contains at least one Lambda with minimum pT 
	    fh1RCBiasLa->Fill(1.);//for normalisation purposes
	    for(Int_t it=0; it<fListLaRC->GetSize(); ++it){ // loop for Lambdas in random cone
	      
	      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLaRC->At(it));
	      if(!v0) continue;
	      
	      Double_t invMLa =0;
	      Double_t trackPt=0;
	      Double_t fEta=0;
	      fEta = v0->Eta();
	    
	      CalculateInvMass(v0, kLambda, invMLa, trackPt);  //function to calculate invMass with TLorentzVector class
	      
	      //Double_t vLaRC[3] = {invMLa,trackPt,fEta};
	      //fhnLaRCBias->Fill(vLaRC);
	    }
	  }
	
	  
	 
	  if(isBadJetALaRC == kFALSE){ //in case RC contains at least one Antilambda with minimum pT 
	    fh1RCBiasALa->Fill(1.);//for normalisation purposes
	    for(Int_t it=0; it<fListALaRC->GetSize(); ++it){ // loop for Lambdas in random cone
	      
	      AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALaRC->At(it));
	      if(!v0) continue;
	    
	      Double_t invMALa =0;
	      Double_t trackPt=0;
	      Double_t fEta=0;
	      fEta = v0->Eta();
	      
	      CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
	      
	      //Double_t vALaRC[3] = {invMALa,trackPt,fEta};
	      //fhnALaRCBias->Fill(vALaRC);
	    }
	    
	  }

	  
	  GetTracksInCone(fTracksRecCuts, fTracksRecCutsRC, jetRC, GetFFRadius(), sumPtRecCutsRC, 0, 0, isBadJetRecCutsRC);
	  
	  /*for(Int_t in=0; in<fTracksRecCutsRC->GetEntries(); ++in){
	    
	    AliVParticle* track = dynamic_cast<AliVParticle*>(fTracksRecCutsRC->At(in)); //inputlist is fListK0s, all reconstructed K0s in event
	    if(!track){std::cout<<"AliAnalysisTaskJetChem(): In random cones, charged track not found!!!"<<std::endl; continue;}

	    Double_t trackPt = track->Pt();
	    Double_t trackEta = track->Eta();
	    if(TMath::Abs(trackEta) <= fCutEta){
	      fh2ChTracksRC->Fill(trackPt,trackEta);
	    }   
	  }*/

	}
	
	fListK0sRC->Clear();
	fListLaRC->Clear();
	fListALaRC->Clear();
	fTracksRecCutsRC->Clear();
      }


      //fetch particles in perpendicular cone to estimate UE event contribution to particle spectrum
      //these perpendicular cone particle spectra serve to subtract the particles in jet cones, that are stemming from the Underlying event, on a statistical basis
      //for normalization the common jet pT spectrum is used: fh1IMK0Cone, fh1IMLaCone and fh1IMALaCone
      
      //____fetch reconstructed K0s in cone perpendicular to jet axis:_______________________________________________________________________________
      
  
      jetPerpConeK0list->Clear();
      jetPerpRecCutslist->Clear();

      Double_t sumPerpPtK0     = 0.;
      Double_t sumPerpPtRecCuts  = 0.;

      GetTracksInPerpCone(fListK0s, jetPerpConeK0list, jet, GetFFRadius(), sumPerpPtK0); //reconstructed K0s in cone around jet axis
      
      GetTracksInPerpCone(fTracksRecCuts, jetPerpRecCutslist, jet, GetFFRadius(), sumPerpPtRecCuts); //charged tracks in perpendicular cone, used for sytematic uncertainty calculation of UE subtraction

      /*  for(Int_t in=0; in<jetPerpRecCutslist->GetEntries(); ++in){
	    
	AliVParticle* track = dynamic_cast<AliVParticle*>(jetPerpRecCutslist->At(in)); //inputlist is fListK0s, all reconstructed K0s in event
	if(!track){std::cout<<"AliAnalysisTaskJetChem(): In perpendicular cone, charged track not found!!!"<<std::endl; continue;}
	
	Double_t trackPt = track->Pt();
	Double_t trackEta = track->Eta();
	fh2ChTracksPC->Fill(trackPt,trackEta);   	
      }*/

      
      if(fDebug>2)Printf("%s:%d nK0s total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetPerpConeK0list->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetPerpConeK0list->GetSize(); ++it){ // loop for K0s in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeK0list->At(it));
	if(!v0) continue;
	
	Double_t invMPerpK0s =0;
	Double_t trackPt=0;
	Double_t fEta=0;

	fEta = v0->Eta();	
	CalculateInvMass(v0, kK0, invMPerpK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	Double_t vK0sPC[4] = {jetPt, invMPerpK0s,trackPt,fEta};
	
	fhnK0sPC->Fill(vK0sPC);  //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
	
      }
      
      
      if(jetPerpConeK0list->GetSize() == 0){ // no K0s in jet cone 
	
	Double_t vK0sPC[4] = {jetPt, -1, -1 , -999};//default values for case: no K0s is found in PC
	fhnK0sPC->Fill(vK0sPC);
	
      }
      

      if(IsMCCEvt == kFALSE){//median cluster only once for event

	IsMCCEvt = kTRUE;

      // if(ij==0){

	AliAODJet* medianCluster = GetMedianCluster();

      	if(medianCluster){
	// ____ rec K0s in median cluster___________________________________________________________________________________________________________ 
	

	  jetMedianConeK0list->Clear();
	  jetMedianConeLalist->Clear();
	  jetMedianConeALalist->Clear();
	  jetMedianRecCutslist->Clear();//charged jet tracks

	  Double_t medianEta = medianCluster->Eta();
	
	if(TMath::Abs(medianEta)<=fCutjetEta){
	  
	  fh1MedianEta->Fill(medianEta);
	  fh1JetPtMedian->Fill(jetPt);
	  fh1MCC->Fill(1.);//for normalisation by total number of median cluster jets
	  Double_t sumMedianPtK0     = 0.;
	  Double_t sumPtRecCuts      = 0.;
	  
	  Bool_t isBadJetK0Median    = kFALSE; // dummy, do not use
	  Bool_t isBadJetRecCutsMedian = kFALSE; // dummy, do not use

	  GetTracksInCone(fListK0s, jetMedianConeK0list, medianCluster, GetFFRadius(), sumMedianPtK0, 0., 0., isBadJetK0Median); //reconstructed K0s in median cone around jet axis

	  GetTracksInCone(fTracksRecCuts, jetMedianRecCutslist, medianCluster, GetFFRadius(), sumPtRecCuts, 0., 0., isBadJetRecCutsMedian); //charged tracks in median cluster cone, used for sytematic uncertainty calculation of UE subtraction
	  
	  //cut parameters from Fragmentation Function task:
	  //Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value, use GetFFMinLTrackPt()
	  //Float_t fFFMaxTrackPt;    // reject jetscontaining any track with pt larger than this value, use GetFFMaxTrackPt()
	  

	  /*for(Int_t in=0; in<jetMedianRecCutslist->GetEntries(); ++in){
	    
	    AliVParticle* track = dynamic_cast<AliVParticle*>(jetMedianRecCutslist->At(in)); //inputlist is fListK0s, all reconstructed K0s in event
	    if(!track){std::cout<<"AliAnalysisTaskJetChem(): In median cluster  cones, charged track not found!!!"<<std::endl; continue;}

	    Double_t trackPt = track->Pt();
	    Double_t trackEta = track->Eta();
	    if(TMath::Abs(trackEta)<=fCutEta){
	      fh2ChTracksMCC->Fill(trackPt,trackEta);
	    }   	    
	  }*/



	  for(Int_t it=0; it<jetMedianConeK0list->GetSize(); ++it){ // loop for K0s in median cone
	    
	    AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetMedianConeK0list->At(it));
	    if(!v0) continue;
	    
	    Double_t invMMedianK0s =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
	    
	    fEta = v0->Eta();
	    CalculateInvMass(v0, kK0, invMMedianK0s, trackPt);  //function to calculate invMass with TLorentzVector class	
	    Double_t vK0sMCC[3] = {invMMedianK0s,trackPt,fEta};
	    fhnK0sMCC->Fill(vK0sMCC);
	    
	  }
	  
	  if(jetMedianConeK0list->GetSize() == 0){ // no K0s in median cluster cone 
	   
	    Double_t vK0sMCC[3] = {-1, -1, -999};
	    fhnK0sMCC->Fill(vK0sMCC);
	   
	  }
	  
	  //__________________________________________________________________________________________________________________________________________
	  // ____ rec Lambdas in median cluster___________________________________________________________________________________________________________ 
	  
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
	    Double_t fEta=0;
	    fEta = v0->Eta();

	    CalculateInvMass(v0, kLambda, invMMedianLa, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vLaMCC[3] = {invMMedianLa,trackPt,fEta};
	    fhnLaMCC->Fill(vLaMCC);
	  }
	  
	  if(jetMedianConeLalist->GetSize() == 0){ // no Lambdas in median cluster cone 
	   
	    Double_t vLaMCC[4] = {jetPt, -1, -1, -999};
	    fhnLaMCC->Fill(vLaMCC); 
	    
	  }
	  
	
	  // ____ rec Antilambdas in median cluster___________________________________________________________________________________________________________ 
	
	  
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
	    Double_t fEta=0;
	    
	    fEta = v0->Eta();

	    CalculateInvMass(v0, kAntiLambda, invMMedianALa, trackPt);  //function to calculate invMass with TLorentzVector class
	    Double_t vALaMCC[3] = {invMMedianALa,trackPt,fEta};
	    fhnALaMCC->Fill(vALaMCC); 
	    
	  }
	  
	  if(jetMedianConeALalist->GetSize() == 0){ // no Antilambdas in median cluster cone 

	    Double_t vALaMCC[4] = {jetPt, -1, -1, -999};
	    fhnALaMCC->Fill(vALaMCC); 
	    
	  }
	}//median cluster eta cut 
	
        jetMedianConeK0list->Clear();
	jetMedianConeLalist->Clear();
	jetMedianConeALalist->Clear();
	jetMedianRecCutslist->Clear();
	}//if mediancluster is existing
      }//end (IsMCCEvt == kFALSE)
      //_________________________________________________________________________________________________________________________________________
      
      //____fetch reconstructed Lambdas in cone perpendicular to jet axis:__________________________________________________________________________
      
      jetPerpConeLalist->Clear();
      Double_t sumPerpPtLa     = 0.;
      
      GetTracksInPerpCone(fListLa, jetPerpConeLalist, jet, GetFFRadius(), sumPerpPtLa); //reconstructed Lambdas in cone around jet axis //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
      
      if(fDebug>2)Printf("%s:%d nLa total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetPerpConeLalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetPerpConeLalist->GetSize(); ++it){ // loop for Lambdas in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeLalist->At(it));
	if(!v0) continue;
	
	Double_t invMPerpLa =0;
	Double_t trackPt=0;
	Double_t fEta=0;
	fEta = v0->Eta();
	
	CalculateInvMass(v0, kLambda, invMPerpLa, trackPt);  //function to calculate invMass with TLorentzVector class
	Double_t vLaPC[4] = {jetPt, invMPerpLa,trackPt,fEta};
	fhnLaPC->Fill(vLaPC);  //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!

      }
      
      
      if(jetPerpConeLalist->GetSize() == 0){ // no Lambdas in jet
	
	Double_t vLaPC[4] = {jetPt, -1, -1 , -999};//default values for case: no K0s is found in PC
	fhnLaPC->Fill(vLaPC);
	
	
      }
      
      
      //____fetch reconstructed Antilambdas in cone perpendicular to jet axis:___________________________________________________________________
 
      jetPerpConeALalist->Clear();
      Double_t sumPerpPtALa     = 0.;
      
      GetTracksInPerpCone(fListALa, jetPerpConeALalist, jet, GetFFRadius(), sumPerpPtALa); //reconstructed Antilambdas in cone around jet axis //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
      
      if(fDebug>2)Printf("%s:%d nALa total: %d, in perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetPerpConeALalist->GetEntries(),GetFFRadius());
            
      for(Int_t it=0; it<jetPerpConeALalist->GetSize(); ++it){ // loop for ALa in perpendicular cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeALalist->At(it));
	if(!v0) continue;
	
	Double_t invMPerpALa =0;
	Double_t trackPt=0;
	Double_t fEta=0;
	fEta = v0->Eta();

	CalculateInvMass(v0, kAntiLambda, invMPerpALa, trackPt);  //function to calculate invMass with TLorentzVector class
	Double_t vALaPC[4] = {jetPt, invMPerpALa,trackPt,fEta};
	fhnALaPC->Fill(vALaPC);
	
      }
      
      
      if(jetPerpConeALalist->GetSize() == 0){ // no Antilambda 

	Double_t vALaPC[4] = {jetPt, -1, -1, -999};
	fhnALaPC->Fill(vALaPC);

      }
   


      //###########################################################################################################
      //MC Analysis 
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

      //_________________________________________________________________
	for(Int_t it=0; it<fListFeeddownLaCand->GetSize(); ++it){ 
	  
	  AliAODv0* mcfd = dynamic_cast<AliAODv0*>(fListFeeddownLaCand->At(it));
	  if(!mcfd) continue;

	  Double_t invMLaFDcand = 0;
	  Double_t trackPt = 0;//pt of ass. particle, not used for the histos
	  
	  CalculateInvMass(mcfd, kLambda, invMLaFDcand, trackPt);
	  
	  //Get MC gen. Lambda transverse momentum
	  TClonesArray *st = 0x0;

	  if(!fAOD)continue;
	  TList *lt = fAOD->GetList();
	  if(!lt)continue;
	 
	  st = (TClonesArray*)lt->FindObject(AliAODMCParticle::StdBranchName()); //get MCAOD branch in data
	  if (!st)continue;
	  
	  AliAODTrack *daughtertrack = (AliAODTrack *) (mcfd->GetSecondaryVtx()->GetDaughter(0));//fetch the first of the two daughter tracks
	  Int_t AssLabel = TMath::Abs(daughtertrack->GetLabel());

	  AliAODMCParticle *mcDaughterPart =(AliAODMCParticle*)st->UncheckedAt(AssLabel);

	  Int_t v0lab = mcDaughterPart->GetMother();//get v0 particle label  

	  //  Int_t v0lab= TMath::Abs(mcfd->GetLabel());//GetLabel doesn't work for AliAODv0 class!!! Only for AliAODtrack

	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	 
	  Int_t motherlab = mcp->GetMother();  //get mother particle label of v0 particle
  
	  if(motherlab >= 0 && v0lab < st->GetEntriesFast()){                 //do safety check for mother label	    
	    
	  
 
	    Double_t genLaPt = mcp->Pt();
	    
	    Int_t iMother = -1;
	    
	    iMother = mcp->GetMother(); //Motherparticle of V0 candidate (e.g. phi particle,..)
	   
	    if((!iMother) || (iMother<0) || (iMother > st->GetEntriesFast()))continue;//validity checks
	   
	   
	    if(iMother >= 0){
	      
	    
	      AliAODMCParticle *partMother = (AliAODMCParticle*)st->UncheckedAt(iMother);
	      Int_t codeMother = -1;
	      if(!partMother) continue;
	      
	 
	      if(partMother) codeMother = TMath::Abs(partMother->GetPdgCode());
	      
	      //  3312    Xi-     -3312    Xibar+          
	      //  3322    Xi0     -3322    Xibar0 
	      
	      if((codeMother == 3312)||(codeMother == 3322)){// feeddown for Lambda coming from Xi- and Xi0
		
		Double_t XiPt = partMother->Pt();//MC gen. pt
		
		Double_t vFeedDownLa[3] = {5., invMLaFDcand, genLaPt};
		fhnFeedDownLa->Fill(vFeedDownLa);
		fh2FeedDownXiLa->Fill(XiPt,genLaPt);
	      }
	    }
	  }
	}//end loop over feeddown candidates for Lambda particles in jet cone
	//fetch MC truth in jet cones, denominator of rec. efficiency in jet cones
	//_________________________________________________________________
	for(Int_t it=0; it<jetConeFDLalist->GetSize(); ++it){ 
	  
	  AliAODv0* mcfd = dynamic_cast<AliAODv0*>(jetConeFDLalist->At(it));
	  if(!mcfd) continue;
	  
	  //std::cout<<"Cone, recLaPt:"<<mcfd->Pt()<<std::endl;
	  
	  Double_t invMLaFDcand = 0;
	  Double_t trackPt = mcfd->Pt();//pt of ass. particle, not used for the histos
	  
	  CalculateInvMass(mcfd, kLambda, invMLaFDcand, trackPt);
	  
	  //Get MC gen. Lambda transverse momentum
	  TClonesArray *st = 0x0;
	 	 	  
	  TList *lt = fAOD->GetList();
	  if(!lt)continue;
	  
	  st = (TClonesArray*)lt->FindObject(AliAODMCParticle::StdBranchName());
	  
	  AliAODTrack *daughtertrack = (AliAODTrack *) (mcfd->GetSecondaryVtx()->GetDaughter(0));//fetch the first of the two daughter tracks
	  Int_t AssLabel = TMath::Abs(daughtertrack->GetLabel());
	  
	  AliAODMCParticle *mcDaughterPart =(AliAODMCParticle*)st->UncheckedAt(AssLabel);
	  
	  Int_t v0lab = mcDaughterPart->GetMother(); 
	  
	  //std::cout<<"v0lab: "<<v0lab<<std::endl;
	  
	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	    
	  Double_t genLaPt = mcp->Pt();
	  

	  //std::cout<<"Cone FD, genLaPt:"<<genLaPt<<std::endl;

	  Double_t vFeedDownLaCone[3] = {jetPt, invMLaFDcand, genLaPt};
	  fhnFeedDownLaCone->Fill(vFeedDownLaCone);

	  
	}//end loop over feeddown candidates for Lambda particles in jet cone
	
	//_________________________________________________________________
	for(Int_t it=0; it<fListFeeddownALaCand->GetSize(); ++it){ 
	  
	  AliAODv0* mcfd = dynamic_cast<AliAODv0*>(fListFeeddownALaCand->At(it));
	  if(!mcfd) continue;

	  Double_t invMALaFDcand = 0;
	  Double_t trackPt = 0;//pt of ass. particle, not used for the histos
	  
	  CalculateInvMass(mcfd, kAntiLambda, invMALaFDcand, trackPt);
	  
	  //Get MC gen. Antilambda transverse momentum
	  TClonesArray *st = 0x0;
	 	 	  
	  TList *lt = fAOD->GetList();
	  if(!lt)continue;
	 
	  st = (TClonesArray*)lt->FindObject(AliAODMCParticle::StdBranchName());
	
	  AliAODTrack *daughtertrack = (AliAODTrack *) (mcfd->GetSecondaryVtx()->GetDaughter(0));//fetch the first of the two daughter tracks
	  Int_t AssLabel = TMath::Abs(daughtertrack->GetLabel());
	  
	  AliAODMCParticle *mcDaughterPart =(AliAODMCParticle*)st->UncheckedAt(AssLabel);
	  
	  Int_t v0lab = mcDaughterPart->GetMother(); 
	  
	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	   
	  Int_t motherlab = mcp->GetMother();  //get mother particle label of v0 particle

	  if(motherlab >= 0 && v0lab < st->GetEntriesFast()){                 //do safety check for mother label
	    
	   
	    Double_t genALaPt = mcp->Pt();
	    
	    Int_t iMother = -1;
	    
	    iMother = mcp->GetMother(); //Motherparticle of V0 candidate (e.g. phi particle,..)
	    

	    if((!iMother) || (iMother<0) || (iMother > st->GetEntriesFast()))continue;//validity checks

	    if(iMother >= 0){
	      
	      AliAODMCParticle *partMother = (AliAODMCParticle*)st->UncheckedAt(iMother);
	      Int_t codeMother = -1;
	      if(!partMother) continue;
	      
	      if(partMother) codeMother = TMath::Abs(partMother->GetPdgCode());
	      
	      //  3312    Xi-     -3312    Xibar+          
	      //  3322    Xi0     -3322    Xibar0 
	      
	      if((codeMother == -3312)||(codeMother == -3322)){// feeddown for Antilambda coming from Xibar+ and Xibar0
		
		Double_t XiPt = partMother->Pt();//MC gen. pt
		
		Double_t vFeedDownALa[3] = {5., invMALaFDcand, genALaPt};
		fhnFeedDownALa->Fill(vFeedDownALa);
		fh2FeedDownXiALa->Fill(XiPt,genALaPt);
	      }
	    }
	  }
	  
	}//end loop over feeddown candidates for Antilambda particles
	
	
	//_________________________________________________________________
	//feeddown for Antilambdas from Xi(bar)+ and Xi(bar)0 in jet cone:

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
	  
	  AliAODTrack *daughtertrack = (AliAODTrack *) (mcfd->GetSecondaryVtx()->GetDaughter(0));//fetch the first of the two daughter tracks
	  Int_t AssLabel = TMath::Abs(daughtertrack->GetLabel());
	  
	  AliAODMCParticle *mcDaughterPart =(AliAODMCParticle*)st->UncheckedAt(AssLabel);
	  
	  Int_t v0lab = mcDaughterPart->GetMother(); 
	  
	  if((!v0lab) || (v0lab<0) || (v0lab > st->GetEntriesFast()))continue;//validity checks

	  AliAODMCParticle *mcp=(AliAODMCParticle*)st->UncheckedAt(v0lab);
	    
	  Double_t genALaPt = mcp->Pt();
	
	  Double_t vFeedDownALaCone[3] = {jetPt, invMALaFDcand, genALaPt};
	  fhnFeedDownALaCone->Fill(vFeedDownALaCone);

	  
	}//end loop over feeddown candidates for Antilambda particles in jet cone
	
      	

	//____fetch MC generated K0s in cone around jet axis__(note: particles can stem from fragmentation but also from underlying event)________
	
	Double_t sumPtMCgenK0s   = 0.;
	Bool_t isBadJetMCgenK0s  = kFALSE; // dummy, do not use
	
	
	fListMCgenK0sCone->Clear(); //MC generated K0s in (only geometrical) jet cone (these are MC gen K0s falling geometrically into jet cone (R = 0.4) around jet axis, that was found by anti-kt jet finder, particles can stem from fragmentation but also from underlying event!!)
	
	//first: sampling MC gen K0s       
	
	GetTracksInCone(fListMCgenK0s, fListMCgenK0sCone, jet, GetFFRadius(), sumPtMCgenK0s, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMCgenK0s); //MC generated K0s in cone around jet axis 
	
	if(fDebug>2)Printf("%s:%d nMCgenK0s in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,fListMCgenK0sCone->GetEntries(),GetFFRadius());
	
	
        for(Int_t it=0; it<fListMCgenK0sCone->GetSize(); ++it){ // loop MC generated K0s in cone around jet axis
	  
	  AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0sCone->At(it));
	  if(!mcp0) continue;
	  
	  //Double_t fRapMCgenK0s   = MyRapidity(mcp0->E(),mcp0->Pz());//get rec. particle in cone information
	  Double_t fEtaMCgenK0s   = mcp0->Eta();
	  Double_t fPtMCgenK0s    = mcp0->Pt();
	  
	  Double_t vMCgenK0Cone[3] = {jetPt,fPtMCgenK0s,fEtaMCgenK0s};

	  fhnMCgenK0Cone->Fill(vMCgenK0Cone); 
	  //fhnMCEtagenK0Cone->Fill(jetPt,fEtaMCgenK0s);
	  
	  }
	
	//check whether the reconstructed K0s in jet cone are stemming from MC gen K0s (on MCgenK0s list):__________________________________________________
	/*
	for(Int_t ic=0; ic<jetConeK0list->GetSize(); ++ic){    //loop over all reconstructed K0s in jet cone
 	   
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

	  if(!v0c) continue;
	  
	  Bool_t daughtercheck = DaughterTrackCheck(v0c, nnum, pnum);//check daughter tracks have proper sign
	  if(daughtercheck == kFALSE)continue;
	  
	  const AliAODTrack *trackMCNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));
	  const AliAODTrack *trackMCPos=(AliAODTrack *)(v0c->GetDaughter(pnum));
	  Bool_t isinjected;	  
	  TString generatorName;
	  TList *listmc = fAOD->GetList();
	 	  
	  Bool_t mclabelcheck = MCLabelCheck(v0c, kK0, trackMCNeg, trackMCPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode, generatorName, isinjected);
	  		  
	  if(mclabelcheck == kFALSE)continue; 
	  if(fPhysicalPrimary == kFALSE)continue;  //requirements for rec. V0 associated to MC true primary particle

	  for(Int_t it=0; it<fListMCgenK0s->GetSize(); ++it){                                    // loop over MC generated K0s in event, check whether associated MC particle is part of it
	  
	    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    //AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0sCone->At(it));
	    AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenK0s->At(it));
	    if(!mcp0) continue;
	    
	    Bool_t particleMatching = IsParticleMatching(mcp0, v0Label);
	    
	    if(particleMatching == kFALSE)continue;    //if reconstructed V0 particle doesn't match to the associated MC particle go to next stack entry	    
	    CalculateInvMass(v0c, kK0, invMK0Match, fPtMCrecK0Match);
	    Double_t fEta = v0c->Eta();
	    Double_t fPtMCgenK0s    = mcp0->Pt();//pt has to be always MC truth value!
	    
	    Double_t vMCrecK0Cone[4] = {jetPt, invMK0Match,fPtMCgenK0s,fEta};
	    fhnMCrecK0Cone->Fill(vMCrecK0Cone);             //fill matching rec. K0s in 3D histogram

	    //SmearJetPt(jetPt,cl,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);       //jetPt, cent, jetRadius, ptmintrack, &jetPtSmear 
  	 
	    //Double_t vMCrecK0ConeSmear[4] = {jetPtSmear, invMK0Match,fPtMCgenK0s,fEta};
	    //fhnMCrecK0ConeSmear->Fill(vMCrecK0ConeSmear); 

	    //fill matching rec. K0s in 3D histogram, jet pT smeared according to deltaptjet distribution width  
  

	  } // end MCgenK0s / MCgenK0sCone loop
 
	  //___________
	  //check the K0s daughters contamination of the jet tracks:
	  	
	  //TClonesArray *stackMC = 0x0;
	   
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

		Double_t vK0sSecContinCone[3] = {jetPt, trackPosPt, trackPosEta};
		fhnK0sSecContinCone->Fill(vK0sSecContinCone);}           //if it's the case, fill jet pt, daughter track pt and track eta in histo 
	      
	      if(particleLabel == negAssLabel){
		AliAODMCParticle* mctrackNeg = dynamic_cast<AliAODMCParticle*>(stackMC->At(negAssLabel));
		if(!mctrackNeg) continue;
		Double_t trackNegPt = mctrackNeg->Pt();
		Double_t trackNegEta = mctrackNeg->Eta();
	
		Double_t vK0sSecContinCone[3] = {jetPt, trackNegPt, trackNegEta};
		fhnK0sSecContinCone->Fill(vK0sSecContinCone);}              //if it's the case, fill jet pt, daughter track pt and track eta in histo
	    }
	  }
	  	  
	    
	  //_______________
	  
	  
      } //end rec-K0-in-cone loop*/
	
	//________________________________________________________________________________________________________________________________________________________
	  
	fListMCgenK0sCone->Clear();
	
	
      }//end fAnalysisMC
      
      jetConeK0list->Clear();
      
      jetPerpConeK0list->Clear();
      jetPerpConeLalist->Clear();
      jetPerpConeALalist->Clear();
 

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
	//	Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	//if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	//fFFHistosIMLaJet->FillFF(trackPt, invMLa, jetPt, incrementJetPt);
 	
	if(dPhiJetLa<fh1dPhiJetLa->GetXaxis()->GetXmin()) dPhiJetLa += 2*TMath::Pi();
	fh1dPhiJetLa->Fill(dPhiJetLa);
      }

      /*   if(fListLa->GetSize() == 0){ // no La: increment jet pt spectrum 
	
	   Bool_t incrementJetPt = kTRUE;
	   fFFHistosIMLaJet->FillFF(-1, -1, jetPt, incrementJetPt);
	}*/
	
  
      // ____fetch rec. Lambdas in cone around jet axis_______________________________________________________________________________________
      
      jetConeLalist->Clear();
      Double_t sumPtLa     = 0.;
      Bool_t isBadJetLa    = kFALSE; // dummy, do not use

      GetTracksInCone(fListLa, jetConeLalist, jet, GetFFRadius(), sumPtLa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetLa);//method inherited from FF

      if(fDebug>2)Printf("%s:%d nLa total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetConeLalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetConeLalist->GetSize(); ++it){ // loop La in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeLalist->At(it));
	if(!v0) continue;                     

	Int_t nnum;
	Int_t pnum;
	
	Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
	if(daughtercheck == kFALSE)continue;                                                 
	Double_t invMLa =0;
	Double_t trackPt=0;
	Double_t fEta = 0;

	fEta = v0->Eta();

        CalculateInvMass(v0, kLambda, invMLa, trackPt); //function to calculate invMass with TLorentzVector class
	
	Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;//needed for all histos, which serve for normalisation
	
	if(fAnalysisMC){
	  
	  //Int_t negDaughterpdg;
	  //Int_t posDaughterpdg;
	  //Int_t motherType;
	  //Int_t v0Label;
	  Double_t jetPtSmear = -1;  
	  //Double_t MCPt;
	  //Bool_t fPhysicalPrimary = -1;
	  //Int_t MCv0PDGCode =0;
	  //TString generatorName;

	  SmearJetPt(jetPt,jetPtSmear); 
	  
	  fh2MC2LaCone->Fill(jetPtSmear,trackPt);

	  if(incrementJetPt == kTRUE){fh1IMLaConeSmear->Fill(jetPtSmear);


	    /*
	    const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	    const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));	
	    
	    TList *listmc = fAOD->GetList();
	    Bool_t isinjected;	    
	    Bool_t mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode, generatorName, isinjected);
	    if(mclabelcheck == kFALSE)continue;

	    //std::cout<<"generatorName: "<<generatorName<<std::endl;
	    
	    if(generatorName == "Hijing"){
	    Double_t vrecMCHijingLaCone[4] = {jetPt, invMLa,trackPt,fEta};	    
	    fhnrecMCHijingLaCone->Fill(vrecMCHijingLaCone);
	    }

	    if(isinjected == kTRUE){
	    Double_t vrecMCInjectLaCone[4] = {jetPt, invMLa,trackPt,fEta};
	    fhnrecMCInjectLaCone->Fill(vrecMCInjectLaCone);
	    }
	    */

	    }//fill TH1F for normalization purposes 
	}//end MC analysis part
	

	if(incrementJetPt==kTRUE){
	  fh1IMLaCone->Fill(jetPt);}//normalisation by number of selected jets
	
	//fFFHistosIMLaCone->FillFF(trackPt, invMLa, jetPt, incrementJetPt);   
  	Double_t vLaCone[4] = {jetPt, invMLa,trackPt,fEta};
	fhnLaCone->Fill(vLaCone);     
      }

      if(jetConeLalist->GetSize() == 0){ // no La: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	//	fFFHistosIMLaCone->FillFF(-1, -1, jetPt, incrementJetPt);
	Double_t vLaCone[4] = {jetPt, -1, -1, -1};
	fhnLaCone->Fill(vLaCone);
	
	if(incrementJetPt==kTRUE){
	  fh1IMLaCone->Fill(jetPt);}//normalisation by number of selected jets

	if(fAnalysisMC){ 
	  Double_t jetPtSmear;  
	  SmearJetPt(jetPt,jetPtSmear);  
	  if(incrementJetPt == kTRUE){
	    fh1IMLaConeSmear->Fill(jetPtSmear);
	  
	  }
	  }

      }
      
      if(fAnalysisMC){
	
	//____fetch MC generated Lambdas in cone around jet axis__(particles can stem from fragmentation but also from underlying event)_____________
	
	Double_t sumPtMCgenLa      = 0.;
	Bool_t isBadJetMCgenLa  = kFALSE; // dummy, do not use 
	
	//sampling MC gen. Lambdas in cone around reconstructed jet axis      

	fListMCgenLaCone->Clear();
	GetTracksInCone(fListMCgenLa, fListMCgenLaCone, jet, GetFFRadius(), sumPtMCgenLa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMCgenLa);//fetch MC generated Lambdas in cone of resolution parameter R around jet axis 
	
	if(fDebug>2)Printf("%s:%d nMCgenLa in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,fListMCgenLaCone->GetEntries(),GetFFRadius());
	
	for(Int_t it=0; it<fListMCgenLaCone->GetSize(); ++it){ // loop MC generated La in cone around jet axis
	  
	  AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLaCone->At(it));
	  if(!mcp0) continue;
	  
	  //Double_t fRapMCgenLa   = MyRapidity(mcp0->E(),mcp0->Pz());
	  Double_t fEtaMCgenLa   = mcp0->Eta();
	  Double_t fPtMCgenLa    = mcp0->Pt();
    
	  Double_t vMCgenLaCone[3] = {jetPt,fPtMCgenLa,fEtaMCgenLa};

	  fhnMCgenLaCone->Fill(vMCgenLaCone); 

	  //fh2MCgenLaCone->Fill(jetPt,fPtMCgenLa);
	  //fh2MCEtagenLaCone->Fill(jetPt,fEtaMCgenLa);
	  }
	
	/*
	//check whether the reconstructed La are stemming from MC gen La on fListMCgenLa List:__________________________________________________
	
	for(Int_t ic=0; ic<jetConeLalist->GetSize(); ++ic){//loop over all reconstructed La within jet cone, new definition
	  
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
	  TString generatorName;

	  AliAODv0* v0c = dynamic_cast<AliAODv0*>(jetConeLalist->At(ic));//new definition

	 
	  if(!v0c) continue;
	  
	  Bool_t daughtercheck = DaughterTrackCheck(v0c, nnum, pnum);
	  if(daughtercheck == kFALSE)continue;
	  
	  const AliAODTrack *trackMCNeg=(AliAODTrack *)(v0c->GetDaughter(nnum));
	  const AliAODTrack *trackMCPos=(AliAODTrack *)(v0c->GetDaughter(pnum));	

	  TList *listmc = fAOD->GetList();
	  Bool_t isinjected;	  
	  Bool_t mclabelcheck = MCLabelCheck(v0c, kLambda, trackMCNeg, trackMCPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode, generatorName, isinjected);

	  if(mclabelcheck == kFALSE)continue;
	  if(fPhysicalPrimary == kFALSE)continue;
	  
	  for(Int_t it=0; it<fListMCgenLa->GetSize(); ++it){//new definition                                  // loop over MC generated K0s in cone around jet axis


	    //Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    
	    AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLa->At(it));//new definition
	    //AliAODMCParticle* mcp0 = dynamic_cast<AliAODMCParticle*>(fListMCgenLaCone->At(it));//old definition
	    
	    if(!mcp0) continue;
	    
	    Bool_t particleMatching = IsParticleMatching(mcp0, v0Label);
	    	    
	    
	    if(particleMatching == kFALSE)continue; //particle doesn't match on any associated MC gen particle in cone around rec jet axis
	    
	    CalculateInvMass(v0c, kLambda, invMLaMatch, fPtMCrecLaMatch);
	  
	    Double_t fPtMCgenLa    = mcp0->Pt();
	    Double_t fEta          = v0c->Eta();//rec. MC particle
	    Double_t vMCrecLaCone[4] = {jetPt, invMLaMatch,fPtMCgenLa,fEta};
	    fhnMCrecLaCone->Fill(vMCrecLaCone); 

	    SmearJetPt(jetPt,cl,GetFFRadius(),GetFFMinLTrackPt(),jetPtSmear);

	    Double_t vMCrecLaConeSmear[4] = {jetPtSmear, invMLaMatch,fPtMCgenLa,fEta};
	    fhnMCrecLaConeSmear->Fill(vMCrecLaConeSmear);  //fill matching rec. Lambdas in 3D histogram, jet pT smeared according to deltaptjet distribution width     
	        

	  } // end MCgenLa loop
	  
	    //check the Lambda daughters contamination of the jet tracks://///////////////////////////////////////////////////////////////////////////////////////////
	  	
	  // TClonesArray *stackMC = 0x0;
	  	  
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
		Double_t vLaSecContinCone[3] = {jetPt, trackPosPt, trackPosEta};
		fhnLaSecContinCone->Fill(vLaSecContinCone);
		
	      }       //if it's the case, fill jet pt, daughter track pt and track eta in histo 
	            
	       
	      if(particleLabel == negAssLabel){

		AliAODMCParticle* mctrackNeg = dynamic_cast<AliAODMCParticle*>(stackMC->At(negAssLabel));
		if(!mctrackNeg) continue;

		Double_t trackNegPt = trackNeg->Pt();
		Double_t trackNegEta = trackNeg->Eta();
		
		Double_t vLaSecContinCone[3] = {jetPt, trackNegPt, trackNegEta};
		fhnLaSecContinCone->Fill(vLaSecContinCone);

		
	      }              //if it's the case, fill jet pt, daughter track pt and track eta in histo
	    }
	    }

	  		    
      } //end rec-La-in-cone loop
	  */
	//________________________________________________________________________________________________________________________________________________________
	
	fListMCgenLaCone->Clear();
	
      }//end fAnalysisMC
      
      jetConeLalist->Clear();
         
      
 
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
	//Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;

	//if(incrementJetPt){fh1V0JetPt->Fill(jetPt);}

	//fFFHistosIMALaJet->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
 	
	if(dPhiJetALa<fh1dPhiJetALa->GetXaxis()->GetXmin()) dPhiJetALa += 2*TMath::Pi();
	fh1dPhiJetALa->Fill(dPhiJetALa);
      }

      // if(fListALa->GetSize() == 0){ // no ALa: increment jet pt spectrum 
	
      //	Bool_t incrementJetPt = kTRUE;
	//fFFHistosIMALaJet->FillFF(-1, -1, jetPt, incrementJetPt);
      //}
	
  
      // ____fetch rec. Antilambdas in cone around jet axis_______________________________________________________________________________________
      
      jetConeALalist->Clear();
      Double_t sumPtALa     = 0.;
      Bool_t isBadJetALa    = kFALSE; // dummy, do not use

      GetTracksInCone(fListALa, jetConeALalist, jet, GetFFRadius(), sumPtALa, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetALa);//method inherited from FF
      
      if(fDebug>2)Printf("%s:%d nALa total: %d, in jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetConeALalist->GetEntries(),GetFFRadius());
      
      for(Int_t it=0; it<jetConeALalist->GetSize(); ++it){ // loop ALa in jet cone
	
	AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeALalist->At(it));
	if(!v0) continue;                    
	

	Int_t nnum;
	Int_t pnum; 

	Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
	if(daughtercheck == kFALSE)continue; 
	
                                          
	Double_t invMALa =0;
	Double_t trackPt=0;
	Double_t fEta = 0;

	fEta = v0->Eta();

        CalculateInvMass(v0, kAntiLambda, invMALa, trackPt); //function to calculate invMass with TLorentzVector class
	
	  Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;

          if(fAnalysisMC){    //jet pt smearing study for Antilambdas
	    /*
	  Int_t negDaughterpdg;
	  Int_t posDaughterpdg;
	  Int_t motherType;
	  Int_t v0Label;
	  Double_t MCPt;
	  Bool_t fPhysicalPrimary = -1;
	  Int_t MCv0PDGCode =0;
	  TString generatorName;
	    */
	  Double_t jetPtSmear = -1;

	  SmearJetPt(jetPt,jetPtSmear);
	  fh2MC2ALaCone->Fill(jetPtSmear, trackPt);//fill MC true particles for jet pt smearing reference

	  /*  
	  const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	  const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));	
	    
	    TList *listmc = fAOD->GetList();
	    Bool_t isinjected;	    
	    Bool_t mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PDGCode, generatorName, isinjected);
	    if(mclabelcheck == kFALSE)continue;

	    //std::cout<<"generatorName: "<<generatorName<<std::endl;
	    
	    if(generatorName == "Hijing"){
	    Double_t vrecMCHijingALaCone[4] = {jetPt, invMALa,trackPt,fEta};	    
	    fhnrecMCHijingALaCone->Fill(vrecMCHijingALaCone);
	    }

	    if(isinjected == kTRUE){
	    Double_t vrecMCInjectALaCone[4] = {jetPt, invMALa,trackPt,fEta};
	    fhnrecMCInjectALaCone->Fill(vrecMCInjectALaCone);
	    }*/

	  if(incrementJetPt == kTRUE){fh1IMALaConeSmear->Fill(jetPtSmear);}                          //fill TH1F for normalization purposes 
	}//end fAnalysisMC
       

	if(incrementJetPt==kTRUE){
	  fh1IMALaCone->Fill(jetPt);}//normalisation by number of selected jets

	//fFFHistosIMALaCone->FillFF(trackPt, invMALa, jetPt, incrementJetPt);
	Double_t vALaCone[4] = {jetPt, invMALa,trackPt,fEta};
	fhnALaCone->Fill(vALaCone);
      }

      if(jetConeALalist->GetSize() == 0){ // no ALa: increment jet pt spectrum 
	
	Bool_t incrementJetPt = kTRUE;
	
	if(incrementJetPt==kTRUE){
	  fh1IMALaCone->Fill(jetPt);}//normalisation by number of selected jets

	//fFFHistosIMALaCone->FillFF(-1, -1, jetPt, incrementJetPt);
	Double_t vALaCone[4] = {jetPt, -1, -1, -1};
	fhnALaCone->Fill(vALaCone);

	if(fAnalysisMC){ 
	  Double_t jetPtSmear;  
	  SmearJetPt(jetPt,jetPtSmear);  
	  if(incrementJetPt == kTRUE)fh1IMALaConeSmear->Fill(jetPtSmear);}

      }
      
        
      jetConeALalist->Clear();
      jettracklist->Clear();
    }//end of if 'leading' or 'all jet' requirement
  }//end of detector level BIG jet loop
  
    //##################################################################################################################
    
    //########generated jets for embedding##############################################################################
    
  if((fBranchGenJets.Length())&&(fUseExtraTracks == 1)&&(fMatchMode == 2)&&(fIsNJEventEmb == kTRUE)){//match mode needed for V0 histograms, to be running as a seperate wagon for match mode 1 and match mode 2 and only for extra jet branch + switch for Embedding into NJ events (here: events in which all jets were rejected - events with no rec. jets at all are not used here for technical reasons (but should be small amount in PbPb anyhow))

    //match mode 1 is for detector level - detector level PYTHIA matching
    //match mode 2 is for detector level - particle level PYTHIA matching, particle jets can be plotted with true jet pT or smeared jet pT (fUseExtraJetPt)
  
    // generated jets
    for(Int_t ij=0; ij<nGenJets; ++ij){ // gen jets loop, particle level
	
	AliAODJet* jet = dynamic_cast<AliAODJet*>(fJetsGen->At(ij));
	if(!jet)continue;
	
	TList* mclist = fAOD->GetList();
	if (!mclist){std::cout<<"mclist does not exist for Embedding study: "<<std::endl; continue;}


	Double_t ptFractionEmbeddedMC = 0; 
	Double_t deltaREmbeddedMC     = 0;
	Double_t ptFractionEmbedded   = 0; 
	Double_t deltaREmbedded       = 0;
	AliAODJet* embeddedJet        = 0; // jet from detector level PYTHIA tracks 
	AliAODJet* matchedJet         = 0; // jet from UE + detector level PYTHIA tracks
	
	if(fBranchEmbeddedJets.Length()){ // find embedded jet
	  
	  Int_t indexEmbeddedMC = iGenMatchIndex[ij];
	  ptFractionEmbeddedMC  = fGenMatchPtFraction[ij]; 
	  
	  
	  fh1FractionPtEmbeddedMC->Fill(ptFractionEmbeddedMC);//yes!
	  
	  if(fDebug > 2)std::cout<<" ij: "<<ij<<" indexEmbeddedMC: "<<indexEmbeddedMC<<std::endl;
	  
	  fh1IndexEmbeddedMC->Fill(indexEmbeddedMC);//yes!

	  if(indexEmbeddedMC > -1){ 
	    
	    embeddedJet = dynamic_cast<AliAODJet*>(fJetsEmbedded->At(indexEmbeddedMC));
	    if(!embeddedJet) continue;
	    
	    deltaREmbeddedMC = jet->DeltaR((AliVParticle*) (embeddedJet)); 
	    
	    Int_t indexExtra = iEmbeddedMatchIndex[indexEmbeddedMC]; 
	    
	    if(fDebug > 2)std::cout<<" ij "<<ij<<" deltaREmbeddedMC "<<deltaREmbeddedMC<<" indexExtra "<<indexExtra<<std::endl;//yes!
	    
	    if(indexExtra > -1){
	      
	      matchedJet = dynamic_cast<AliAODJet*>(fJetsRecCuts->At(indexExtra));
	      
	      ptFractionEmbedded = fRecMatchPtFraction[indexExtra]; 
	      deltaREmbedded     = embeddedJet->DeltaR((AliVParticle*) (matchedJet)); //yes!
	      
	      if(fDebug > 2)std::cout<<"In gen. jet loop - jet matching - ij: "<<ij<<" indexExtra: "<<indexExtra<<" ptFractionEmbedded: "<<ptFractionEmbedded<<" deltaREmbedded: "<<deltaREmbedded<<std::endl;//yes! This is the last printed statement
	

	    }	
	  }
	}
	
	TList* jettrackList = new TList();//gen. jets track list
	Double_t sumPt      = 0.;
	Bool_t isBadJet     = kFALSE;

	TList* jettrackListMatch = new TList();//matched jets track list
	Double_t sumPtMatch      = 0.;
	Bool_t isBadJetMatch     = kFALSE;

	//gen. jet tracks:
	if(GetFFRadius()<=0){
	  GetJetTracksTrackrefs(jettrackList, jet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);

	}
	else GetJetTracksPointing(fTracksGen, jettrackList, jet, GetFFRadius(), sumPt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJet);
	
	if(GetFFMinNTracks()>0 && jettrackList->GetSize()<=GetFFMinNTracks()){isBadJet = kTRUE;}
	
	if(jettrackList->GetEntries() == 0){
	  
	  if(fDebug >2)std::cout<<" Generated jet loop: jettrackList is empty! "<<std::endl;
	  
	  delete jettrackList;
	  continue; 
	}
	

	jettrackListMatch->Clear();
	
	if(matchedJet){
	  //matched jets tracks:
	  if(GetFFRadius()<=0){
	    GetJetTracksTrackrefs(jettrackListMatch, matchedJet, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMatch);
	    
	  }
	  else GetJetTracksPointing(fTracksRecCuts, jettrackListMatch, matchedJet, GetFFRadius(), sumPtMatch, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetMatch);}
	
	
	if(GetFFMinNTracks()>0 && jettrackListMatch->GetSize()<=GetFFMinNTracks()){isBadJetMatch = kTRUE;}
	
	if(jettrackListMatch->GetEntries() == 0){
	  
	  if(fDebug >2)std::cout<<" Generated jet loop: jettrackListMatch is empty! "<<std::endl;
	  
	  delete jettrackListMatch;
	  continue; 
	}
	
	// embedding QA after leading pt cut
	
	if(!embeddedJet)std::cout<<"Gen. jet loop: no embedded jet (in extra branch) existing!! "<<std::endl; 
	if(!matchedJet)std::cout<<"Gen. jet loop: no matched jet (in extra branch) existing!! "<<std::endl; 

	//Double_t jetPt = jet->Pt();

        if((fDebug > 2) && embeddedJet && matchedJet) cout<<" After leading pt cut: gen jet "<<ij<<" pt "<<jet->Pt()<<" embedded jet pt "<<embeddedJet->Pt()<<" matched jet pt "<<matchedJet->Pt()
				      <<" ptFractionEmbeddedMC "<<ptFractionEmbeddedMC<<" dRMC "<<deltaREmbeddedMC
				      <<" ptFractionEmbedded "<<ptFractionEmbedded<<" dR "<<deltaREmbedded<<endl;     //no!
	
	
	if(embeddedJet){ 

	  fh2FractionPtVsEmbeddedJetPtMC->Fill(embeddedJet->Pt(),ptFractionEmbeddedMC); //yes!
	  
	  if(ptFractionEmbeddedMC>=fCutFractionPtEmbedded){    
	    
	    fh1DeltaREmbeddedMC->Fill(deltaREmbeddedMC); //yes!
	  }
	}
	
	//apply both requirements: matching from rec. extra jets to detector level PYTHIA AND particle level PYTHIA
	if(ptFractionEmbeddedMC>=fCutFractionPtEmbedded && deltaREmbeddedMC <= fCutDeltaREmbedded && 
	   ptFractionEmbedded>=fCutFractionPtEmbedded && deltaREmbedded <= fCutDeltaREmbedded) { // if no embedding: ptFraction = cutFraction = 0

	  Double_t embJetPt = embeddedJet->Pt();//jet pt detector level (matched to generator level PYTHIA jets)
	  
	  fh1JetPtEmbGenAfterMatch->Fill(embJetPt);  //no! but maybe only a matter of statistics..
	  
	  if(fDebug > 2)std::cout<<" After MatchMode 2 matching cuts - embJetPt: "<<embJetPt<<std::endl; 

	  //charged tracks as crosscheck to Olivers results	
	  for(Int_t it=0; it<jettrackList->GetSize(); ++it){
	  
	    AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(jettracklist->At(it));
	    if(!trackVP)continue;
	    TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
	    
	    Float_t jetPt = jet->Pt(); //can be set in extra branch instance of task to smear the particle level jet reference
	    if(fUseExtraJetPt){
	      if(matchedJet) jetPt = matchedJet->Pt();
	      else jetPt = 0;
	    }
	    
	    if(fDebug > 2)std::cout<<" After MatchMode 2 matching cuts - jetPt: "<<jetPt<<std::endl; 

	    Float_t trackPt = trackV->Pt();
	    
	    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
	    
	    if (ij==0) fFFHistosGen->FillFF(trackPt, jetPt, incrementJetPt);

	    //fFFHistosGenInc->FillFF(trackPt, jetPt, incrementJetPt);//can be either PYTHIA particle level jet pT, or smeared with matched jet pT
	   
	    
	    delete trackV;
	  }
	  
	  
	  //V0 analyse with 'gen. PYTHIA - rec. extra jets' - matching ###################################################
	  
	  FillEmbeddedHistos(embeddedJet, matchedJet, nK0s, nLa, nALa, jettrackListMatch);//fill all V0 embedding histos for match mode 2

	  //Fill gen. jet V0s:
	  
	  Double_t sumPtK0EmbMC     = 0.;
	  Bool_t isBadJetK0EmbMC    = kFALSE; // dummy, do not use
 

	
	  if(fListK0sMC->GetEntries() > 0){GetTracksInCone(fListK0sMC, jetConeK0EmbMClist, jet, GetFFRadius(), sumPtK0EmbMC, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0EmbMC);} //reconstructed K0s in cone around jet axis
	   	
	  
	  if(fDebug>2)Printf("%s:%d nK0s total: %d, in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetConeK0EmbMClist->GetEntries(),GetFFRadius());

	  if(fUseExtraTracks == 1){//only for extra particles used	
	    
	    
	    //MC gen PYTHIA Antilambdas in jet cone
	    
	    
	    for(int it =0; it<jetConeK0EmbMClist->GetEntries(); it++){//loop over particles
	      
	      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(jetConeK0EmbMClist->At(it));
	      
	      if(!part)continue;
	      Double_t genPt = part->Pt(); 
	      
	      Float_t jetPt = jet->Pt(); //can be set in extra branch instance of task to smear the particle level jet reference
	      if(fUseExtraJetPt){
		if(matchedJet) jetPt = matchedJet->Pt();
		else jetPt = 0;
	      }
	   
	      if(fDebug > 2)std::cout<<" gen. EmbCone K0s candidate - partPt: "<<genPt<<" jet Pt: "<<jetPt<<std::endl;
	      
	      fh2MCEmbK0sJetPt->Fill(jetPt,genPt);   
	      
	    } 
	    
	    

	    //Lambdas in particle level jet cone
	    jetConeLaEmbMClist->Clear();
	    
	    Double_t sumPtLaEmbMC     = 0.;
	    
	    Bool_t isBadJetLaEmbMC    = kFALSE; // dummy, do not use
	    
	    
	    if(fMatchMode == 2){
	      if(fListLaMC->GetEntries() > 0){GetTracksInCone(fListLaMC, jetConeLaEmbMClist, jet, GetFFRadius(), sumPtLaEmbMC, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetLaEmbMC);} //reconstructed La in cone around jet axis
	    } 
	    
	    if(fDebug>2)Printf("%s:%d nLa total: %d, in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetConeLaEmbMClist->GetEntries(),GetFFRadius());
	    
	    
	    //MC gen PYTHIA Lambdas in jet cone
	      
	    
	    for(int it =0; it<jetConeLaEmbMClist->GetEntries(); it++){//loop over particles
	      
	      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(jetConeLaEmbMClist->At(it));
	      
		if(!part)continue;
		Double_t genPt = part->Pt();
		
		Float_t jetPt = jet->Pt(); //can be set in extra branch instance of task to smear the particle level jet reference
		if(fUseExtraJetPt){
		  if(matchedJet) jetPt = matchedJet->Pt();
		  else jetPt = 0;
		}
		
		
		if(fDebug > 2)std::cout<<" gen. EmbCone Lambda candidate - partPt: "<<genPt<<" jet Pt: "<<jetPt<<std::endl;
		
		fh2MCEmbLaJetPt->Fill(jetPt,genPt);   
		
	    } 
	    
	    
	    jetConeALaEmbMClist->Clear();
	    
	    Double_t sumPtALaEmbMC     = 0.;
	    
	    Bool_t isBadJetALaEmbMC    = kFALSE; // dummy, do not use
	    
	    
	    if(fListALaMC->GetEntries() > 0){GetTracksInCone(fListALaMC, jetConeALaEmbMClist, jet, GetFFRadius(), sumPtALaEmbMC, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetALaEmbMC);} //reconstructed La in cone around jet axis
	    
	    
	    if(fDebug>2)Printf("%s:%d nALa total: %d, in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetConeALaEmbMClist->GetEntries(),GetFFRadius());
	    
	    if(fMatchMode == 2) {//MC gen PYTHIA Antilambdas in jet cone
	      
	      for(int it =0; it<jetConeALaEmbMClist->GetEntries(); it++){//loop over particles
		
		AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(jetConeALaEmbMClist->At(it));
		
		if(!part)continue;
         
		Float_t jetPt = jet->Pt(); //can be set in extra branch instance of task to smear the particle level jet reference
		if(fUseExtraJetPt){
		  if(matchedJet) jetPt = matchedJet->Pt();
		  else jetPt = 0;
		}
		
		Double_t genPt = part->Pt();   
		if(fDebug > 2)std::cout<<" gen. EmbCone ALambda candidate - partPt: "<<genPt<<" jet Pt: "<<jetPt<<std::endl;
		
		fh2MCEmbALaJetPt->Fill(jetPt,genPt);    
		
	      } 
	    }
	    
	  } //end extra
	}
	
	delete jettrackList;
	delete jettrackListMatch;
    }
  }
  
  //########end of generated jets for embedding
  //#########################################################################
  
 
  jettracklist->Clear();
  jetConeK0list->Clear();
  jetConeLalist->Clear();
  jetConeALalist->Clear();
  jetConeK0Emblist->Clear();
  jetConeLaEmblist->Clear();
  jetConeALaEmblist->Clear();

  jetConeK0EmbStlist->Clear();
  jetConeLaEmbStlist->Clear();
  jetConeALaEmbStlist->Clear();

  jetConeK0EmbMClist->Clear();
  jetConeLaEmbMClist->Clear();
  jetConeALaEmbMClist->Clear();

  jetPerpConeK0list->Clear();
  jetPerpConeK0Emblist->Clear();
  jetPerpConeLalist->Clear();
  jetPerpConeLaEmblist->Clear();
  jetPerpConeALalist->Clear();
  jetPerpConeALaEmblist->Clear();
  jetMedianConeK0list->Clear();
  jetMedianConeLalist->Clear();
  jetMedianConeALalist->Clear();
  jetMedianRecCutslist->Clear();
  fListK0sRC->Clear();
  fListLaRC->Clear();
  fListALaRC->Clear();
  fTracksRecCutsRC->Clear();
  fTracksRecCuts->Clear();
  fTracksPerpCone->Clear();
  fTracksGen->Clear();
  fJetsRecCuts->Clear();
  fJetsGen->Clear();
  fJetsEmbedded->Clear();
  fBckgJetsRec->Clear();
  fListK0s->Clear();
  fListLa->Clear();
  fListALa->Clear();

  fListK0sMC->Clear();
  fListLaMC->Clear();
  fListALaMC->Clear();

  fListK0sStandard->Clear();
  fListLaStandard->Clear();
  fListALaStandard->Clear();
  fListFeeddownLaCand->Clear();
  fListFeeddownALaCand->Clear();
  jetConeFDLalist->Clear();
  jetConeFDALalist->Clear();
  fListMCgenK0s->Clear();
  fListMCgenLa->Clear();
  fListMCgenALa->Clear();
  fListMCgenK0sCone->Clear();
  fListMCgenLaCone->Clear();
  fListMCgenALaCone->Clear();
  
  //Post output data.
  PostData(1, fCommonHistList); 
  //end of event loop
   
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
Int_t AliAnalysisTaskJetChem::GetListOfV0s(TList *list, const Int_t type, const Int_t particletype, const Int_t tracktype, AliAODVertex* primVertex, AliAODEvent* aod)
{
  // fill list of V0s selected according to type
  
  if(!list){
    if(fDebug>1) Printf("%s:%d no input list", (char*)__FILE__,__LINE__);
    return -1;
  }
    
  if(fDebug>5){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s(): type: "<<type<<" particletype: "<<particletype<<"aod: "<<aod<<std::endl;
    if(type==kV0TypeUndef){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s(): kV0TypeUndef!! "<<std::endl;}
  }
   
  if(type==kV0TypeUndef) return 0;
  Int_t iCount = 0;

  if(!primVertex) return 0;// this is real vertex of PbPb event
    
    Double_t lPrimaryVtxPosition[3];
    Double_t lV0Position[3];
    lPrimaryVtxPosition[0] = primVertex->GetX();
    lPrimaryVtxPosition[1] = primVertex->GetY();
    lPrimaryVtxPosition[2] = primVertex->GetZ();

  if(fDebug>5){ std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s(): aod->GetNumberOfV0s: "<<aod->GetNumberOfV0s()<<std::endl;}

  //First the embedded v0s:########################################################################################################

  if((tracktype==kTrackAODExtraonlyCuts)||(tracktype==kTrackAODExtraCuts)){//carried out only when extra or extraonly tracks are requested

    TClonesArray *aodExtrav0s = 0x0;
        
    if(particletype == kK0){aodExtrav0s = dynamic_cast<TClonesArray*>(fAOD->FindListObject("aodExtraK0s"));}//search for embedded V0s

    if(particletype == kLambda){aodExtrav0s = dynamic_cast<TClonesArray*>(fAOD->FindListObject("aodExtraLa"));}

    if(particletype == kAntiLambda){aodExtrav0s = dynamic_cast<TClonesArray*>(fAOD->FindListObject("aodExtraALa"));}

    if(!aodExtrav0s){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s(): aodExtraV0s list object not found!!!!"<<std::endl; return iCount;}

    //std::cout<<"aodExtrav0s->GetEntries(): "<<aodExtrav0s->GetEntries()<<std::endl;

    //be careful, handle properly both vertices (PYTHIA, real)! PYTHIA v0s are cutted already in AliAnalysisFastEmbedding class since they have their own primary vertex.

    for(int it =0; it<aodExtrav0s->GetEntries(); it++) {
      AliVParticle *extrav0 = dynamic_cast<AliVParticle*> ((*aodExtrav0s)[it]);
      if (!extrav0){continue;}
      
      AliAODv0 *v0 = dynamic_cast<AliAODv0*> (extrav0);
      if(!v0)continue;
	
       Bool_t isOnFly = v0->GetOnFlyStatus();
       
       if(!isOnFly &&  (type == kOnFly || type == kOnFlyPID || type == kOnFlydEdx || type == kOnFlyPrim)) continue; 
       if( isOnFly &&  (type == kOffl  || type == kOfflPID  || type == kOffldEdx  || type == kOfflPrim))  continue; 
       

       Double_t invM = 0;
       Double_t invMK0s=0;
       Double_t invMLa=0;
       Double_t invMALa=0;
       Double_t trackPt=0;
       Int_t nnum = -1;
       Int_t pnum = -1;
        
       Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
       
       if(daughtercheck == kFALSE)continue;
 
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
       

       if(particletype == kK0)fh1V0PtCandidate->Fill(trackPt);//only used for x-checks


       /////////////////////////////////////////////////////////////
       //V0 and track Cuts:
              
       if(fDebug>20){if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa))){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s: invM not in selected mass window "<<std::endl;}}

       if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa)))continue; 
       
       //all other PYTHIA v0 particle cuts are already applied in AliAnalysisTaskFastEmbedding.cxx ->cut values are set in AddTask_aod_FastEmbedding.C! pay attention!  
       //since also pdg code check is applied before embedding, in principle no cuts are needed for the PYTHIA V0s
      list->Add(v0);
      iCount++;
      
    }
  }

  //##############################################################################################################

  if(tracktype==kV0AODMCExtraonlyCuts){//carried out only for Match Mode 2: particle to detector level matching

    TClonesArray *aodExtraMCparticles = 0x0;
        
    aodExtraMCparticles = dynamic_cast<TClonesArray*>(fAOD->FindListObject("aodExtraMCparticles"));//search for MC gen. V0s in PYTHIA

 
    if(!aodExtraMCparticles){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s(): aodExtraMCparticles array not found!!!!"<<std::endl; return iCount;}


    if(fDebug > 2){std::cout<<"aodExtraMCparticles->GetEntries(): "<<aodExtraMCparticles->GetEntries()<<std::endl;}

    for(int it =0; it<aodExtraMCparticles->GetEntries(); it++){//loop over particles

      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(aodExtraMCparticles->At(it));

      if(!part){std::cout<<"MC stack extra particle does not exist!!"<<std::endl;continue;}
      if(!part->IsPhysicalPrimary()){std::cout<<"MC stack extra particle fails IsPrimaryParticle check!!"<<std::endl;continue;}
     
      Int_t pdgCode = part->GetPdgCode();

     if(fDebug > 5)std::cout<<"MC stack extra particles pdgCode: "<<pdgCode<<std::endl;

      //keep only true primary particles from PYTHIA particle level
      if((particletype == kK0)&&(pdgCode != 310))continue;
      if((particletype == kLambda)&&(pdgCode != 3122))continue;
      if((particletype == kAntiLambda)&&(pdgCode != -3122))continue;
 
      Double_t fEta = part->Eta();      

      if(TMath::Abs(fEta) > fCutEta) continue;
    

      list->Add(part);
      iCount++;
    }
  }


    //Last the standard V0s ######################################################################################

  if((tracktype==kTrackAODExtraCuts)||(tracktype==kTrackAODCuts)){//all v0s (standard plus embedded v0s) or: standard v0s (no embedding)
    
    // std::cout<<"standard aod->GetNumberOfV0s(): "<<aod->GetNumberOfV0s()<<std::endl;
    
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
       
       
       if(fDebug>20){if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa))){std::cout<<"AliAnalysisTaskJetChem::GetListOfV0s: invM not in selected mass window "<<std::endl;}}
       
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
      
       Double_t fEta = -999;
       Double_t fV0cosPointAngle = -999;
       Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
       Double_t fV0Rap = -999;
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
       
       if(particletype == kK0)         {//fRap = v0->RapK0Short();
	 fEta = v0->PseudoRapV0();}
       if(particletype == kLambda)     {//fRap = v0->RapLambda();
	 fEta = v0->PseudoRapV0();}
       if(particletype == kAntiLambda) {//fRap = v0->Y(-3122);
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
       //Double_t PosPt;
       //Double_t NegPt;   
       
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
	 
	 //PosPt = v0Pos.Perp(v0totMom);             //longitudinal momentum of positive charged daughter track
	 PosPl = v0Pos.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of positive charged daughter track
	 
	 //NegPt = v0Neg.Perp(v0totMom);             //longitudinal momentum of negative charged daughter track
	 NegPl = v0Neg.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of nergative charged daughter track
	 
	 ArmenterosAlpha = 1.-2./(1+(PosPl/NegPl));  
	 ArmenterosPt= v0->PtArmV0();
	 
       }      
       
       if(particletype == kK0){//only cut on K0s histos
	 
	 if(IsArmenterosSelected == 1){// Armenteros Cut to reject Lambdas contamination in K0s inv. massspectrum
	   fh2ArmenterosBeforeCuts->Fill(ArmenterosAlpha,ArmenterosPt);
	 }
	 if (fV0cosPointAngle < fCutK0cosPointAngle)	continue;                                       //cospointangle cut
	 
       }
       
       //some more cuts on v0s and daughter tracks:
       
       
       if((TMath::Abs(PosEta)>fCutPostrackEta) || (TMath::Abs(NegEta)>fCutNegtrackEta))continue;   //Daughters pseudorapidity cut
       
       if(particletype != kK0){//only cut on La and ALa histos
	 if (fV0cosPointAngle < fCutLacosPointAngle)continue;                                       //cospointangle cut
       }

       if(particletype == kK0){  
	 fV0Rap = v0->RapK0Short();
       }
       
       if((particletype == kLambda)&&(particletype == kAntiLambda)){  
	 fV0Rap = v0->RapLambda();
       }

          
       //fV0Rap   = MyRapidity(v0->AliAODRecoDecay::E(),v0->Pz());//does not work anymore in class

       if ((fCutRap > 0) && (TMath::Abs(fV0Rap) > fCutRap))continue;      //V0 Rapidity Cut

       if ((fCutEta > 0) && (TMath::Abs(fEta) > fCutEta))continue;     //V0 Eta Cut
                                                    
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
       
       
       //cut on 2D DistOverTransMom
       if(particletype == kK0){//the cut on Lambdas you can find above
	 
	 fh2ProperLifetimeK0sVsPtBeforeCut->Fill(trackPt,fMROverPtK0s); //fill these histos after all other cuts
	 if(fMROverPtK0s > (fCutV0DecayMax * avDecayLengthK0s))continue;
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
	   TString generatorName;
	   Bool_t isinjected;
	   if(particletype == kLambda){
	     
	     mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
	     
	     
	     if((motherType == 3312)||(motherType == 3322)){//mother of v0 is neutral or negative Xi
	    
	       fListFeeddownLaCand->Add(v0); //fill TList with ass. particles, stemming from feeddown from Xi(bar) decays       	    
	     }
	   }
	   
	   if(particletype == kAntiLambda){
	     
	     mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
	     
	     if((motherType == -3312)||(motherType == -3322)){
	       fListFeeddownALaCand->Add(v0); //fill TList with ass. particles, stemming from feeddown from Xi(bar) decays	      	          
	     }
	   }
	 }
	 
	 
	 //_only true primary particles survive the following checks_______________________________________________________________________________________________
	 TString generatorName;
	 Bool_t isinjected;
	 if(particletype == kK0){
	   
	   mclabelcheck = MCLabelCheck(v0, kK0, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
	   if(mclabelcheck == kFALSE)continue;	 
	 }
	 if(particletype == kLambda){
	   
	   mclabelcheck = MCLabelCheck(v0, kLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
	   if(mclabelcheck == kFALSE)continue;	
	 }
	 if(particletype == kAntiLambda){
	   
	   mclabelcheck = MCLabelCheck(v0, kAntiLambda, trackNeg, trackPos, listmc, negDaughterpdg, posDaughterpdg, motherType, v0Label, MCPt, fPhysicalPrimary, MCv0PdgCode, generatorName, isinjected);
	   if(mclabelcheck == kFALSE)continue;	 
	 }
	 
	 if(fPhysicalPrimary != 1)continue; //V0 candidate (K0s, Lambda or Antilambda) must be physical primary, this means there is no mother particle existing
	 
       }
       
       
       list->Add(v0);

     } 
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
  //Double_t fPtCurrentPart   = 0; //get transverse momentum
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
      
      //fRapCurrentPart   = MyRapidity(p0->E(),p0->Pz());

      fRapCurrentPart = p0->Y();

      fEtaCurrentPart   = p0->Eta();
      //fPtCurrentPart    = p0->Pt();
            

      if ((fCutRap > 0) && (TMath::Abs(fRapCurrentPart) >= fCutRap))continue;


      if ((fCutEta > 0) && (TMath::Abs(fEtaCurrentPart) >= fCutEta))continue;
      	
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


void AliAnalysisTaskJetChem::GetTracksInCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, 
								const Double_t radius, Double_t& sumPt, const Double_t minPt, const Double_t maxPt, Bool_t& isBadPt)
{
  // fill list of V0 tracks in cone around jet axis  

  sumPt = 0;
  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;   

  Double_t jetMom[3];
  if(!jet)return;
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);


  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){//loop over all K0s found in event

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){//fill all the V0s inside cone into outputlist, radius is return value of GetFFRadius() 

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
 
  Double_t jetMom[3];             //array for entries in TVector3
  Double_t perpjetplusMom[3];     //array for entries in TVector3
  Double_t perpjetnegMom[3];
 
  if(!jet)return;

  jet->PxPyPz(jetMom);   //get 3D jet momentum
  Double_t jetPerpPt = jet->Pt(); //original jet pt, invariant under rotations
  Double_t jetPhi = jet->Phi(); //original jet phi

  Double_t jetPerpposPhi = jetPhi + ((TMath::Pi())*0.5);//get new perp. jet axis phi clockwise
  Double_t jetPerpnegPhi = jetPhi - ((TMath::Pi())*0.5);//get new perp. jet axis phi counterclockwise

  TVector3 jet3mom(jetMom); //3-Vector for original rec. jet axis
 
  //Double_t phitest  = jet3mom.Phi();
  
  perpjetplusMom[0]=(TMath::Cos(jetPerpposPhi)*jetPerpPt); //x coordinate (sidewards - when looking in beam direction)
  perpjetplusMom[1]=(TMath::Sin(jetPerpposPhi)*jetPerpPt); //y coordinate (upwards - when looking in beam direction)
  perpjetplusMom[2]=jetMom[2];                          //z coordinate (along beam axis), invariant under azimuthal rotation
  
  perpjetnegMom[0]=(TMath::Cos(jetPerpnegPhi)*jetPerpPt); //x coordinate (sidewards - when looking in beam direction)
  perpjetnegMom[1]=(TMath::Sin(jetPerpnegPhi)*jetPerpPt); //y coordinate (upwards - when looking in beam direction)
  perpjetnegMom[2]=jetMom[2];                          //z coordinate (along beam axis), invariant under azimuthal rotation       
  
  
  TVector3 perpjetplus3mom(perpjetplusMom);  //3-Vector for new perp. jet axis, clockwise rotated
  TVector3 perpjetneg3mom(perpjetnegMom);    //3-Vector for new perp. jet axis, counterclockwise rotated
  
  //for crosscheck TVector3 rotation method

  //Double_t jetMomplusTest[3];         
  //Double_t jetMomminusTest[3]; 

  //jet3mom.RotateZ(TMath::Pi()*0.5);//rotate original jet axis around +90 degrees in phi

  //perpjetminus3momTest  = jet3mom.RotateZ((-1)*TMath::Pi()*0.5);

  // jet3mom.RotateZ(TMath::Pi()*0.5);
  // jet3mom.RotateZ((-1)*TMath::Pi()*0.5);
  
  //jetMomplusTest[0] = jet3mom.X();  //fetching perp. axis coordinates
  //jetMomplusTest[1] = jet3mom.Y();
  //jetMomplusTest[2] = jet3mom.Z();

  //TVector3 perpjetplus3momTest(jetMomplusTest);   //new TVector3 for +90deg rotated jet axis with rotation method from ROOT
  //TVector3 perpjetminus3momTest(jetMomminusTest);   //new TVector3 for -90deg rotated jet axis with rotation method from ROOT


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

Bool_t AliAnalysisTaskJetChem::MCLabelCheck(AliAODv0* v0, Int_t particletype,const AliAODTrack* trackNeg, const AliAODTrack* trackPos, TList *listmc, Int_t& negDaughterpdg, Int_t& posDaughterpdg, Int_t& motherType, Int_t& v0Label, Double_t& MCPt, Bool_t& fPhysicalPrimary, Int_t& MCv0PDGCode, TString& generatorName, Bool_t& isinjected){
                               
  if(!v0)return kFALSE;

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

      

      //mc label checks
     
      if(negAssLabel>=0 && negAssLabel < stackmc->GetEntriesFast() && posAssLabel>=0 && posAssLabel < stackmc->GetEntriesFast()){//safety check if label has valid value of stack

	AliAODMCParticle *mcNegPart =(AliAODMCParticle*)stackmc->UncheckedAt(negAssLabel);//fetch the, with one MC truth track associated (reconstructed), negative charged track 
	v0Label = mcNegPart->GetMother();// mother of negative charged particle is v0, get v0 label here
	negDaughterpdg = mcNegPart->GetPdgCode();
	AliAODMCParticle *mcPosPart =(AliAODMCParticle*)stackmc->UncheckedAt(posAssLabel);//fetch the, with one MC truth track associated (reconstructed), positive charged track 
	Int_t v0PosLabel = mcPosPart->GetMother();                                        //get mother label of positive charged track label
	posDaughterpdg = mcPosPart->GetPdgCode();

	if(v0Label >= 0 && v0Label < stackmc->GetEntriesFast() && v0Label == v0PosLabel){//first v0 mc label check, then: check if both daughters are stemming from same particle
  
	  AliAODMCParticle *mcv0 = (AliAODMCParticle *)stackmc->UncheckedAt(v0Label);  //fetch MC ass. particle to v0 (mother of the both charged daughter tracks)
	 
	  Float_t fDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)
	  
	  // Get the distance between production point of the MC mother particle and the primary vertex
	 
	  Double_t dx = mcXv-mcv0->Xv();//mc primary vertex - mc particle production vertex 
	  Double_t dy = mcYv-mcv0->Yv();
	  Double_t dz = mcZv-mcv0->Zv();
	  
	  Float_t fDistPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	  fPhysicalPrimary = kFALSE;//init

	  fPhysicalPrimary = (fDistPrimary < fDistPrimaryMax);
	  MCv0PDGCode = mcv0->GetPdgCode();

	  //if(fPhysicalPrimary == kTRUE){//look only at physical primary particles

	  isinjected = IsTrackInjected(v0Label, header, stackmc, generatorName);

	    //trackinjected is kFALSE if it is either Hijing or has no generator name
	    // std::cout<<" "<<std::endl;
	    // std::cout<<"#### next particle: ####"<<std::endl;
	    //std::cout<<"Is track injected: "<< trackinjected <<std::endl;
	    // std::cout<<"pdg code: "<<MCv0PDGCode<<std::endl; 
	    // std::cout<<"v0Label: "<<v0Label<<std::endl;
	  	  
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
	       
	       if(MCv0PDGCode != 310)  {return kFALSE;}
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

//______________________________________________________________________
TString AliAnalysisTaskJetChem::GetGenerator(Int_t label, AliAODMCHeader* header){
  Int_t nsumpart=0;//number of particles
  TList *lh=header->GetCocktailHeaders();//TList with all generator headers
  Int_t nh=lh->GetEntries();//number of entries in TList with all headers

   for(Int_t i=0;i<nh;i++){
     AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
     TString genname=gh->GetName();//name of particle generator
     Int_t npart=gh->NProduced();//number of stable or undecayed particles in MC stack block (?)
     if(label>=nsumpart && label<(nsumpart+npart)) return genname;
     nsumpart+=npart;
   }
   TString empty="";
   return empty;
 }

//_____________________________________________________________________
void AliAnalysisTaskJetChem::GetTrackPrimaryGenerator(Int_t lab, AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  // method to check if a particle is stemming from a given generator

  nameGen=GetGenerator(lab,header);
  
  //  Int_t countControl=0;
  
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);//get MC generated particle for particle MC label
    if(!mcpart){
      printf("AliAnalysisTaskJetChem::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();

    if(mother<0){
      printf("AliAnalysisTaskJetChem::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=GetGenerator(mother,header);

    // countControl++;
    // if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
    //   printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
    //   break;
    // }

  }
  
  return;
}


//---------------------------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskJetChem::IsTrackInjected(Int_t lab, AliAODMCHeader *header,TClonesArray *arrayMC, TString& nameGen){
  // method to check if a v0 particle comes from the signal event or from the underlying Hijing event
  //TString nameGen;

  GetTrackPrimaryGenerator(lab, header, arrayMC, nameGen);
  
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;//particle has either no info about generator or is Hijing particle, so it is not injected
  
  //std::cout<<"generator name: "<<nameGen<<std::endl;

  return kTRUE;
}

//_________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetChem::SmearJetPt(Double_t jetPt, Double_t& jetPtSmear){	   
  
  static TF1 fsmear("f1","[0]*exp(-1*(x-[1])*(x-[1])/(2*[2]*[2]))",-100.,100.);   //smearing according to gaussian function in between  +/- 10 GeV/c
 

  fsmear.SetParameters(1,0,5.10);//for 2010 PbPb jets, R=0.4, ptmintrack = 0.15 GeV/c, cent 00-10%, delta-pt width estimated via single track embedding

  
  Double_t r = fsmear.GetRandom();
  jetPtSmear = jetPt + r;
  
  if(fDebug > 3){
    std::cout<<"jetPt: "<<jetPt<<std::endl;
    std::cout<<"jetPtSmear: "<<jetPtSmear<<std::endl;
    std::cout<<"r: "<<r<<std::endl;
  }
  
  return jetPtSmear;
}


//______________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskJetChem::IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const
{
// decides whether a particle is inside a jet cone, or has a minimum distance to a second track/axis
  if (!part1 || !part2)
    return kFALSE;

  TVector3 vecMom2(part2->Px(),part2->Py(),part2->Pz());
  TVector3 vecMom1(part1->Px(),part1->Py(),part1->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  if(dR<dRMax) // momentum vectors of part1 and part2 are closer than dRMax
    return kTRUE;
  return kFALSE;
}
//__________________________________________________________________________________________________________________


Bool_t AliAnalysisTaskJetChem::IsRCJCOverlap(TList* recjetlist, const AliVParticle* part, Double_t dDistance) const{
  
  if(!recjetlist) return kFALSE;
  if(!part) return kFALSE;
  if(!dDistance) return kFALSE;
  Int_t nRecJetsCuts = fJetsRecCuts->GetEntries(); 
  
  for(Int_t i=0; i<nRecJetsCuts; ++i){   //loop over all reconstructed jets in events      
    AliAODJet* jet = (AliAODJet*) (recjetlist->At(i));
    if(!jet){if(fDebug>2)std::cout<<"AliAnalysisTaskJetChem::IsRCJCOverlap jet pointer invalid!"<<std::endl;continue;}
    if(IsParticleInCone(jet, part, dDistance) == kTRUE)return kTRUE;//RC and JC are overlapping
    
  }//end loop testing RC-JC overlap
  return kFALSE;//RC and JC are not overlapping -> good!
}

//_______________________________________________________________________________________________________________________
AliAODJet* AliAnalysisTaskJetChem::GetRandomCone(TList* jetlist, Double_t dEtaConeMax, Double_t dDistance) const
{ 
  TLorentzVector vecRdCone;
  AliAODJet* jetRC = 0;//random cone candidate
  Double_t dEta, dPhi; //random eta and phi value for RC
  Bool_t IsRCoutJC = kFALSE;//check whether RC is not overlapping with any selected jet cone in event
  Int_t iRCTrials = 10;//search at maximum 10 times for random cone that doesn't overlap with jet cone
  
  for(Int_t i=0; i<iRCTrials; iRCTrials++){
    
    dEta = dEtaConeMax*(2*fRandom->Rndm()-1.); //random eta value in range: [-dEtaConeMax,+dEtaConeMax]
    dPhi = TMath::TwoPi()*fRandom->Rndm(); //random phi value in range: [0,2*Pi]
    vecRdCone.SetPtEtaPhiM(1.,dEta,dPhi,0.);
    jetRC = new AliAODJet(vecRdCone);//new RC candidate

    if (!IsRCJCOverlap(jetlist,jetRC,dDistance))
        {
          IsRCoutJC = kTRUE; //std::cout<<"RC and JC are not overlapping!!!"<<std::endl;
          break;
        }
      else
        delete jetRC; //RC is overlapping with JC, delete this RC candidate
         
  }
  if(!IsRCoutJC) {jetRC = 0;}//in case no random cone was selected

  return jetRC;
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

  if(TMath::Odd(nBckgClusters)){

    Int_t medianIndex = indices[(Int_t) (0.5*(nBckgClusters+1))];

    medianCluster = (AliAODJet*)(fBckgJetsRec->At(medianIndex));
    
    //Double_t clusterPt = medianCluster->Pt();
    //Double_t area      = medianCluster->EffectiveAreaCharged();
  }
  else{

    Int_t medianIndex1 = indices[(Int_t) (0.5*nBckgClusters)];
    Int_t medianIndex2 = indices[(Int_t) (0.5*nBckgClusters+1)];

    AliAODJet* medianCluster1 = (AliAODJet*)(fBckgJetsRec->At(medianIndex1));
    AliAODJet* medianCluster2 = (AliAODJet*)(fBckgJetsRec->At(medianIndex2));
    
    // Double_t density1 = 0;
    //Double_t clusterPt1 = medianCluster1->Pt();
    //Double_t area1      = medianCluster1->EffectiveAreaCharged();
    //if(area1>0) Double_t density1 = clusterPt1/area1;
    
    // Double_t density2 = 0;
    //Double_t clusterPt2 = medianCluster2->Pt();
    //Double_t area2      = medianCluster2->EffectiveAreaCharged();
    // if(area2>0) Double_t density2 = clusterPt2/area2;
    
    medianCluster = ( (gRandom->Rndm()>0.5) ? medianCluster1 : medianCluster2 );  // select one randomly to avoid adding areas
  }
    
  delete[] bgrDensity;
  delete[] indices; 

  return medianCluster;
}    

//____________________________________________________________________________________________

Double_t AliAnalysisTaskJetChem::AreaCircSegment(Double_t dRadius, Double_t dDistance) const
{
// calculate area of a circular segment defined by the circle radius and the (oriented) distance between the secant line and the circle centre
  Double_t dEpsilon = 1e-2;
  Double_t dR = dRadius;
  Double_t dD = dDistance;
  if (TMath::Abs(dR)<dEpsilon)
    {
      if(fDebug>0) printf("AliAnalysisTaskJetChem::AreaCircSegment: Error: Too small radius: %f < %f\n",dR,dEpsilon);
      return 0.;
    }
  if (dD>=dR)
    return 0.;
  if (dD<=-dR)
    return TMath::Pi()*dR*dR;
  return dR*dR*TMath::ACos(dD/dR)-dD*TMath::Sqrt(dR*dR-dD*dD);
}


//_____________________________________________________________________
void AliAnalysisTaskJetChem::FillEmbeddedHistos(const AliAODJet* embeddedJet, const AliAODJet* jet, Int_t nK0s, Int_t nLa, Int_t nALa, TList* Jettracklist)

{//Jet matching cuts (FractionPtEmbedded(MC) = 0.5, DeltaREmb(MC) = 0.75*R) are already applied
  
  //for extraonly branch the embedded jet is empty and 'jet' is detector level PYTHIA jet
   
  if(!Jettracklist){std::cout<<"FillEmbeddedHistos(): Jettracklist does not exist!"<<std::endl;}
  if(Jettracklist->GetEntries() == 0)std::cout<<"FillEmbeddedHistos(): Jettracklist() is empty!!"<<std::endl;

  Double_t jetPt = 0;

  if(jet) jetPt = jet->Pt();//getting matched jet pt

  if((fUseExtraTracks == 1)&&(fUseEmbeddedJetPt == kTRUE)){
  if(jet) jetPt = embeddedJet->Pt();//getting matched jet pt
  if((fDebug > 2)&&(embeddedJet))std::cout<<"jetPt is embedded JetPt:  "<<jetPt<<std::endl;
  }

  if((fUseExtraTracks == 1)&&(fUseEmbeddedJetPt == kFALSE)){
  if(jet) jetPt = jet->Pt();//getting matched jet pt
  if((fDebug > 2)&&(jet))std::cout<<"jetPt is matched JetPt:  "<<jetPt<<std::endl;
  }

  fh1PtEmbAfterMatch->Fill(jetPt);
  
  for(Int_t it=0; it<Jettracklist->GetSize(); ++it){
    
    AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(Jettracklist->At(it));
    if(!trackVP)continue;
    
    AliAODTrack * aodtrack  = dynamic_cast<AliAODTrack*>(Jettracklist->At(it));
    if(!aodtrack) continue;
    
    TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
     
    Float_t trackPt = trackV->Pt();
    
    Bool_t incrementJetPt = (it==0) ? kTRUE : kFALSE;
    
    fFFHistosRecCuts->FillFF(trackPt, jetPt, incrementJetPt);//fill charged tracks into RecCuts histos
    
    delete trackV;
  }
  
  //perpendicular cone method for charged particles in UE
  
  Double_t sumPerpPt = 0.;
  
  GetTracksInPerpCone(fTracksRecCuts,fTracksPerpCone, jet, GetFFRadius(), sumPerpPt);//get charged tracks from UE in perp. cone
  
  if(fTracksPerpCone){
    fh1PerpCone->Fill(1.);}
  
  for(Int_t it=0; it<fTracksPerpCone->GetSize(); ++it){
    
    AliVParticle*   trackVP = dynamic_cast<AliVParticle*>(fTracksPerpCone->At(it));
    if(!trackVP)continue;
    TLorentzVector* trackV  = new TLorentzVector(trackVP->Px(),trackVP->Py(),trackVP->Pz(),trackVP->P());
    
    
    Float_t trackPt = trackV->Pt();
    fh2TracksPerpCone->Fill(jetPt,trackPt);
    
  }//end loop perp. cone charged tracks
  
  if(fTracksPerpCone->GetSize() == 0){ // no K0s in jet cone 
    
    
    fh2TracksPerpCone->Fill(jetPt,0);	    
  }
  
  
  fTracksPerpCone->Clear();
  
  //###############################
  //Charged tracks in outside cones
  
  //Double_t dRadiusExCone = 2*GetFFRadius();//2*radius of jet cone
  
  //Int_t nChargedTracks = fTracksRecCuts->GetEntries(); //number of all charged tracks in event
  
  //###############################
  Jettracklist->Clear();
  
  
  
  //###############################	
  //V0 analysis part for embedding study
  
  
  //____fetch reconstructed K0s in cone around jet axis:_______________________________________________________________________________
  
  jetConeK0Emblist->Clear();
  jetConeLaEmbStlist->Clear();//standard V0s in embedded jets, take smeared jet energy
  jetConeALaEmbMClist->Clear();
  
  Double_t sumPtK0Emb     = 0.;	
  Bool_t isBadJetK0Emb    = kFALSE; // dummy, do not use
  Double_t sumPtK0EmbSt     = 0.;
  Bool_t isBadJetK0EmbSt    = kFALSE; // dummy, do not use
  
  GetTracksInCone(fListK0s, jetConeK0Emblist, jet, GetFFRadius(), sumPtK0Emb, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0Emb); //reconstructed K0s in cone around jet axis
  
  
  if(fListK0sStandard->GetEntries() > 0){GetTracksInCone(fListK0sStandard, jetConeK0EmbStlist, jet, GetFFRadius(), sumPtK0EmbSt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetK0EmbSt);} //reconstructed K0s in cone around jet axis
  
 
  //######
  
  if(fDebug>2)Printf("%s:%d nK0s total: %d, in emb. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetConeK0Emblist->GetEntries(),GetFFRadius());
  if(fDebug>2)Printf("%s:%d nK0s total: %d, standard K0s in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetConeK0EmbStlist->GetEntries(),GetFFRadius());

  if(fUseExtraTracks == 1){//only for extra particles used	
      
    //standard V0s for UE V0 subtraction
    
    for(Int_t it=0; it<jetConeK0EmbStlist->GetSize(); ++it){ // loop for K0s in jet cone
      
      AliAODv0* v0st = dynamic_cast<AliAODv0*>(jetConeK0EmbStlist->At(it));
      if(!v0st) continue;
      
      
      Double_t invMK0s =0;
      Double_t trackPt=0;
      Double_t fEta=0;
      
      fEta = v0st->Eta();
      
      CalculateInvMass(v0st, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
      
      //for embedding V0 UE subtraction
      
      Double_t vK0sEmbConeStandard[4] = {jetPt, invMK0s,trackPt,fEta};
      
      if(fDebug > 7)std::cout<<" v0 standard particle pt: "<<v0st->Pt()<< "eta: "<<fEta<<std::endl;
      
      fhnK0sEmbConeStandard->Fill(vK0sEmbConeStandard);
      
    }	  
      
    if(jetConeK0EmbStlist->GetSize() == 0){ 	
      
      Double_t vK0sEmbConeStandard[4] = {jetPt, 0., 0., -100.};
      fhnK0sEmbConeStandard->Fill(vK0sEmbConeStandard);
    }
  }// only for extra tracks used
  
  //######
  //V0s in embedded jet cones, both extra tracks and extraonly tracks
  
  for(Int_t it=0; it<jetConeK0Emblist->GetSize(); ++it){ // loop for K0s in jet cone
    
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeK0Emblist->At(it));
    if(!v0) continue;
    
    TString generatorName;
    
    Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
    Double_t invMK0s =0;
    Double_t trackPt=0;
    Double_t fEta=0;
    
    fEta = v0->Eta();
    
    CalculateInvMass(v0, kK0, invMK0s, trackPt);  //function to calculate invMass with TLorentzVector class
	  
    //=======
    //for embedding signal calculation and extraonly particle reference
    if(incrementJetPt==kTRUE){
      fh1IMK0EmbCone->Fill(jetPt);
    }//normalisation by number of selected jets
	  
    Double_t vK0sEmbCone[4] = {jetPt, invMK0s,trackPt,fEta};
    
    if(fDebug > 2)std::cout<<"EmbCone K0s candidate: "<<"invMK0s: "<<invMK0s<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
    
    fhnK0sEmbCone->Fill(vK0sEmbCone);
    
    //=========
    //for particle reference via embedding flags and using the extra particles
    if(fUseExtraTracks == 1){
      
      Int_t nnum;
      Int_t pnum;

      Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
      if(daughtercheck == kFALSE)continue;
	    	    
      const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
      const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
      
      if (!trackPos || !trackNeg) {
	if(fDebug > 1) Printf("AliAnalysisTaskJetChem::PYTHIAEmbedding part embedding flag:: Error:Could not retrieve one of the daughter tracks\n");
	continue;
      }
      
      Bool_t isEmbeddedNeg = (Bool_t) (trackNeg->GetFlags() & AliESDtrack::kEmbedded);//check whether daughter particles are stemming from embedded PYTHIA events
      Bool_t isEmbeddedPos = (Bool_t) (trackPos->GetFlags() & AliESDtrack::kEmbedded);
      
      if (isEmbeddedNeg==kTRUE && isEmbeddedPos==kTRUE) {
	if(fDebug > 2) Printf("AliAnalysisTaskJetChem info: PYTHIAEmbedding charged daughters both successfully found with embedding flag!\n");	     
      }
      
      if (isEmbeddedNeg == kTRUE && isEmbeddedPos == kTRUE) {
	
	Double_t vK0sEmbConeRef[4] = {jetPt,invMK0s,trackPt,fEta};
	
	
	if(fDebug > 2)std::cout<<" extra EmbConeRef K0s candidate: "<<"invMK0s: "<<invMK0s<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	
	
	fhnK0sEmbConeRef->Fill(vK0sEmbConeRef);
      }	    
    }
    
    //==========
    
    
  }//end jetConeK0Emblist loop, if list is not empty
  
  //case: if list is empty:
  if(jetConeK0Emblist->GetSize() == 0){ // no K0: increment jet pt spectrum 
    
    Bool_t incrementJetPt = kTRUE;//jets without K0s will be only filled in TH1F only once, so no increment needed 
    //fFFHistosIMK0Cone->FillFF(-1, -1, jetPt, incrementJetPt);
    Double_t vK0sEmbCone[4] = {jetPt, -1, -1, -1};
    fhnK0sEmbCone->Fill(vK0sEmbCone);
    fhnK0sEmbConeRef->Fill(vK0sEmbCone);
    
    if(incrementJetPt==kTRUE){
      fh1IMK0EmbCone->Fill(jetPt);}//normalisation by number of selected jets
  }    
  
  //######
  
  jetPerpConeK0Emblist->Clear();
  Double_t sumPerpPtK0Emb     = 0.;
  
  GetTracksInPerpCone(fListK0s, jetPerpConeK0Emblist, jet, GetFFRadius(), sumPerpPtK0Emb); //reconstructed K0s in cone around jet axis
  
  if(fDebug>2)Printf("%s:%d nK0s total: %d, in emb. perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nK0s,jetPerpConeK0Emblist->GetEntries(),GetFFRadius());
  
  for(Int_t it=0; it<jetPerpConeK0Emblist->GetSize(); ++it){ // loop for K0s in perpendicular cone
    
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeK0Emblist->At(it));
    if(!v0) continue;
    
    Double_t invMPerpK0s =0;
    Double_t trackPt=0;
    Double_t fEta=0;
    
    fEta = v0->Eta();	
    CalculateInvMass(v0, kK0, invMPerpK0s, trackPt);  //function to calculate invMass with TLorentzVector class
    Double_t vK0sPCEmb[4] = {jetPt, invMPerpK0s,trackPt,fEta};
    
    fhnK0sEmbPC->Fill(vK0sPCEmb);  //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
    
  }
  
  
  if(jetPerpConeK0Emblist->GetSize() == 0){ // no K0s in jet cone 
	  
    Double_t vK0sPCEmb[4] = {jetPt, -1, -1 , -999};//default values for case: no K0s is found in PC
    fhnK0sEmbPC->Fill(vK0sPCEmb);
    
  }//end: embedded jets K0s background via perpendicular cones
	
  
  
  jetConeK0Emblist->Clear();
  jetConeK0EmbStlist->Clear();
  jetPerpConeK0Emblist->Clear();
  

  //###########EMBEDDED LAMBDAS AND ANTILAMBDAS#############################################################################################
  
  
  //____fetch reconstructed Lambdas in cone around jet axis:_______________________________________________________________________________
  
  jetConeLaEmblist->Clear();
  jetConeLaEmbStlist->Clear();
  

  Double_t sumPtLaEmb     = 0.;
  
  Bool_t isBadJetLaEmb    = kFALSE; // dummy, do not use
  
  Double_t sumPtLaEmbSt     = 0.;
	
  Bool_t isBadJetLaEmbSt    = kFALSE; // dummy, do not use
  
 
  GetTracksInCone(fListLa, jetConeLaEmblist, jet, GetFFRadius(), sumPtLaEmb, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetLaEmb); //reconstructed La in cone around jet axis
  if(fListLaStandard->GetEntries() > 0){GetTracksInCone(fListLaStandard, jetConeLaEmbStlist, jet, GetFFRadius(), sumPtLaEmbSt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetLaEmbSt);} //reconstructed La in cone around jet axis of extra jet
  
 
  //######
	
  if(fDebug>2)Printf("%s:%d nLa total: %d, in embedded jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetConeLaEmblist->GetEntries(),GetFFRadius());
  if(fDebug>2)Printf("%s:%d nLa total: %d, standard La in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetConeLaEmbStlist->GetEntries(),GetFFRadius());
  
  if(fUseExtraTracks == 1){//only for extra particles used, V0 UE subtraction method
    
   	  //####################

	  for(Int_t it=0; it<jetConeLaEmbStlist->GetSize(); ++it){ // loop for La in jet cone
	    
	    AliAODv0* v0st = dynamic_cast<AliAODv0*>(jetConeLaEmbStlist->At(it));
	    if(!v0st) continue;
	    
        
	    Double_t invMLa =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
	
	    
	    fEta = v0st->Eta();
	    
	    CalculateInvMass(v0st, kLambda, invMLa, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vLaEmbCone[4] = {jetPt, invMLa,trackPt,fEta};
	    	    
	    if(fDebug > 2)std::cout<<" EmbCone Lambda candidate: "<<"invMLa: "<<invMLa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	    

	    fhnLaEmbConeStandard->Fill(vLaEmbCone);
	    
	  }
	  
	  
	  if(jetConeLaEmbStlist->GetSize() == 0){ // no La: increment jet pt entry of spectrum for normalisation
	    
	   
	    Double_t vLaEmbStCone[4] = {jetPt, -1., -1., -100.};
	      
	    fhnLaEmbConeStandard->Fill(vLaEmbStCone);	  
	    	
	  }    

	}


	//########
	  //extra V0s

	for(Int_t it=0; it<jetConeLaEmblist->GetSize(); ++it){ // loop for La in jet cone
	  
	  AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeLaEmblist->At(it));
	  if(!v0) continue;
	
	  TString generatorName;
	
	  Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  Double_t invMLa =0;
	  Double_t trackPt=0;
	  Double_t fEta=0;

	  fEta = v0->Eta();
	  
	  CalculateInvMass(v0, kLambda, invMLa, trackPt);  //function to calculate invMass with TLorentzVector class
	  
      	  //=======
	  //for embedding signal calculation and extraonly particle referenz
	  if(incrementJetPt==kTRUE){
	    fh1IMLaEmbCone->Fill(jetPt);
	  }//normalisation by number of selected jets
	  
	  
	  Double_t vLaEmbCone[4] = {jetPt, invMLa,trackPt,fEta};
	  
	  if(fDebug > 2)std::cout<<"EmbCone Lambda candidate: "<<"invMLa: "<<invMLa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	  
	  
	  fhnLaEmbCone->Fill(vLaEmbCone);
	  
	  	
	  //=========
	  //for particle referenz via embedding flags and using the extra particles
	  if(fUseExtraTracks == 1){
	    	    
	    Int_t nnum;
	    Int_t pnum;

	    Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
	    if(daughtercheck == kFALSE)continue;
	    	    
	    const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	    const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));

	    if (!trackPos || !trackNeg) {
	      if(fDebug > 1) Printf("AliAnalysisTaskJetChem::PYTHIAEmbedding part embedding flag:: Error:Could not retrieve one of the daughter tracks\n");
	      continue;
	    }

	    Bool_t isEmbeddedNeg = (Bool_t) (trackNeg->GetFlags() & AliESDtrack::kEmbedded);//check whether daughter particles are stemming from embedded PYTHIA events
            Bool_t isEmbeddedPos = (Bool_t) (trackPos->GetFlags() & AliESDtrack::kEmbedded);

	    if (isEmbeddedNeg == kTRUE && isEmbeddedPos == kTRUE) {
	          
	      Double_t vLaEmbConeRef[4] = {jetPt,invMLa,trackPt,fEta};
	      
	      if(fDebug > 2)std::cout<<"extra EmbConeRef Lambda candidate: "<<"invMLa: "<<invMLa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	    
	      fhnLaEmbConeRef->Fill(vLaEmbConeRef);
	    }	    
	  }

	  //==========
	  
	}
	
	
	if(jetConeLaEmblist->GetSize() == 0){ // no La: increment jet pt spectrum 
	  
	  Bool_t incrementJetPt = kTRUE;//jets without La will be only filled in TH1F only once, so no increment needed 
	  //fFFHistosIMLaCone->FillFF(-1, -1, jetPt, incrementJetPt);
	  Double_t vLaEmbCone[4] = {jetPt, -1, -1, -1};
	  fhnLaEmbCone->Fill(vLaEmbCone);
	  fhnLaEmbConeRef->Fill(vLaEmbCone);

	  if(incrementJetPt==kTRUE){
	    fh1IMLaEmbCone->Fill(jetPt);}//normalisation by number of selected jets
	}    
	
	
	//######
	
	jetPerpConeLaEmblist->Clear();
	Double_t sumPerpPtLaEmb     = 0.;
	
	GetTracksInPerpCone(fListLa, jetPerpConeLaEmblist, jet, GetFFRadius(), sumPerpPtLaEmb); //reconstructed La in cone around jet axis
	
	if(fDebug>2)Printf("%s:%d nLa total: %d, in emb. perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nLa,jetPerpConeLaEmblist->GetEntries(),GetFFRadius());
	
	for(Int_t it=0; it<jetPerpConeLaEmblist->GetSize(); ++it){ // loop for La in perpendicular cone
	  
	  AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeLaEmblist->At(it));
	  if(!v0) continue;
	  
	  Double_t invMPerpLa =0;
	  Double_t trackPt=0;
	  Double_t fEta=0;
	  
	  fEta = v0->Eta();	
	  CalculateInvMass(v0, kLambda, invMPerpLa, trackPt);  //function to calculate invMass with TLorentzVector class
	  Double_t vLaPCEmb[4] = {jetPt, invMPerpLa,trackPt,fEta};
	  
	  fhnLaEmbPC->Fill(vLaPCEmb);  //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
	  
	}
	
	
	if(jetPerpConeLaEmblist->GetSize() == 0){ // no La in jet cone 
	  
	  Double_t vLaPCEmb[4] = {jetPt, -1, -1 , -999};//default values for case: no La is found in PC
	  fhnLaEmbPC->Fill(vLaPCEmb);
	
	}//end: embedded jets La background via perpendicular cones
	
	

	jetConeLaEmblist->Clear();
	jetPerpConeLaEmblist->Clear();
	

	//____fetch reconstructed Antilambdas in cone around jet axis:_______________________________________________________________________________
	
	jetConeALaEmblist->Clear();
	jetConeALaEmbStlist->Clear();


	Double_t sumPtALaEmb     = 0.;
	
	Bool_t isBadJetALaEmb    = kFALSE; // dummy, do not use

	Double_t sumPtALaEmbSt     = 0.;
	
	Bool_t isBadJetALaEmbSt    = kFALSE; // dummy, do not use


	GetTracksInCone(fListALa, jetConeALaEmblist, jet, GetFFRadius(), sumPtALaEmb, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetALaEmb); //reconstructed ALa in cone around jet axis
	if(fListALaStandard->GetEntries() > 0){GetTracksInCone(fListALaStandard, jetConeALaEmbStlist, jet, GetFFRadius(), sumPtALaEmbSt, GetFFMinLTrackPt(), GetFFMaxTrackPt(), isBadJetALaEmbSt);}//reconstructed ALa in cone around jet axis
	


	//######
	
	if(fDebug>2)Printf("%s:%d nALa total: %d, in embedded jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetConeALaEmblist->GetEntries(),GetFFRadius());
	if(fDebug>2)Printf("%s:%d nALa total: %d, standard ALa in gen. jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetConeALaEmbStlist->GetEntries(),GetFFRadius());


	//standard Antilambdas in jets for UE V0 subtraction

	if(fUseExtraTracks == 1){//only for extra particles used

       	 
	  //####################

	  for(Int_t it=0; it<jetConeALaEmbStlist->GetSize(); ++it){ // loop for ALa in jet cone
	    
	    AliAODv0* v0st = dynamic_cast<AliAODv0*>(jetConeALaEmbStlist->At(it));
	    if(!v0st) continue;
	    
	    Double_t invMALa =0;
	    Double_t trackPt=0;
	    Double_t fEta=0;
        
    	    
	    fEta = v0st->Eta();
	    
	    CalculateInvMass(v0st, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
	    
	    Double_t vALaEmbStCone[4] = {jetPt, invMALa,trackPt,fEta};
	    	    
	    if(fDebug > 2)std::cout<<" standard EmbCone ALambda candidate: "<<"invMALa: "<<invMALa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	    
	    fhnALaEmbConeStandard->Fill(vALaEmbStCone);
	    
	  }
	  

	  if(jetConeALaEmbStlist->GetSize() == 0){ // no ALa: increment jet pt entry of spectrum for normalisation
	    
	    
	    Double_t vALaEmbStCone[4] = {jetPt, -1., -1., -100.};

	    fhnALaEmbConeStandard->Fill(vALaEmbStCone);
	   

	  }    
	  	  
	}


	//################


	for(Int_t it=0; it<jetConeALaEmblist->GetSize(); ++it){ // loop for ALa in jet cone
	  
	  AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetConeALaEmblist->At(it));
	  if(!v0) continue;
	
	  TString generatorName;

	  Bool_t   incrementJetPt = (it==0) ? kTRUE : kFALSE;
	  Double_t invMALa =0;
	  Double_t trackPt=0;
	  Double_t fEta=0;


	  fEta = v0->Eta();
	  
	  CalculateInvMass(v0, kAntiLambda, invMALa, trackPt);  //function to calculate invMass with TLorentzVector class
	  	  
	  //=======
	  //for embedding signal calculation and extraonly particle referenz

	  if(incrementJetPt==kTRUE){
	    fh1IMALaEmbCone->Fill(jetPt);
	  }//normalisation by number of selected jets
	  

	    Double_t vALaEmbCone[4] = {jetPt, invMALa,trackPt,fEta};

	    if(fDebug > 2)std::cout<<"EmbCone ALambda candidate: "<<"invMALa: "<<invMALa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	    
	    fhnALaEmbCone->Fill(vALaEmbCone);
	 
	  
	  //=========
	  //for particle referenz via embedding flags and using the extra particles
	  if(fUseExtraTracks == 1){
	    	    
	    Int_t nnum;
	    Int_t pnum;

	    Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
	    if(daughtercheck == kFALSE)continue;
	    	    
	    const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	    const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));

	    if (!trackPos || !trackNeg) {
	      if(fDebug > 1) Printf("AliAnalysisTaskJetChem::PYTHIAEmbedding part embedding flag:: Error:Could not retrieve one of the daughter tracks\n");
	      continue;
	    }

	    Bool_t isEmbeddedNeg = (Bool_t) (trackNeg->GetFlags() & AliESDtrack::kEmbedded);//check whether daughter particles are stemming from embedded PYTHIA events
            Bool_t isEmbeddedPos = (Bool_t) (trackPos->GetFlags() & AliESDtrack::kEmbedded);

	    if (isEmbeddedNeg == kTRUE && isEmbeddedPos == kTRUE) {
	      
	    
	      Double_t vALaEmbConeRef[4] = {jetPt,invMALa,trackPt,fEta};
	      
	      if(fDebug > 2)std::cout<<" EmbConeRef ALambda candidate: "<<"invMALa: "<<invMALa<<" trackPt: "<<trackPt<<" fEta: "<<fEta<<std::endl;
	    
	      fhnALaEmbConeRef->Fill(vALaEmbConeRef);
	    }	    
	  }

	  //==========
	 
	}
	
	
	if(jetConeALaEmblist->GetSize() == 0){ // no ALa: increment jet pt spectrum 
	  
	  Bool_t incrementJetPt = kTRUE;//jets without ALa will be only filled in TH1F only once, so no increment needed 
	  //fFFHistosIMALaCone->FillFF(-1, -1, jetPt, incrementJetPt);
	  Double_t vALaEmbCone[4] = {jetPt, -1, -1, -1};
	  fhnALaEmbCone->Fill(vALaEmbCone);
	  fhnALaEmbConeRef->Fill(vALaEmbCone);	  

	  if(incrementJetPt==kTRUE){
	    fh1IMALaEmbCone->Fill(jetPt);}//normalisation by number of selected jets
	}    
	
	
	//######
	
	jetPerpConeALaEmblist->Clear();
	Double_t sumPerpPtALaEmb     = 0.;
	
	GetTracksInPerpCone(fListALa, jetPerpConeALaEmblist, jet, GetFFRadius(), sumPerpPtALaEmb); //reconstructed ALa in cone around jet axis
	
	if(fDebug>2)Printf("%s:%d nALa total: %d, in emb. perp jet cone: %d,FFRadius %f ",(char*)__FILE__,__LINE__,nALa,jetPerpConeALaEmblist->GetEntries(),GetFFRadius());
	
	for(Int_t it=0; it<jetPerpConeALaEmblist->GetSize(); ++it){ // loop for ALa in perpendicular cone
	  
	  AliAODv0* v0 = dynamic_cast<AliAODv0*>(jetPerpConeALaEmblist->At(it));
	  if(!v0) continue;
	  
	  Double_t invMPerpALa =0;
	  Double_t trackPt=0;
	  Double_t fEta=0;
	  
	  fEta = v0->Eta();	
	  CalculateInvMass(v0, kAntiLambda, invMPerpALa, trackPt);  //function to calculate invMass with TLorentzVector class
	  Double_t vALaPCEmb[4] = {jetPt, invMPerpALa,trackPt,fEta};
	  
	  fhnALaEmbPC->Fill(vALaPCEmb);  //(x,y,z) //pay attention, this histogram contains the V0 content of both (+/- 90 degrees) perp. cones!!
	  
	}
	
	
	if(jetPerpConeALaEmblist->GetSize() == 0){ // no ALa in jet cone 
	  
	  Double_t vALaPCEmb[4] = {jetPt, -1, -1 , -999};//default values for case: no ALa is found in PC
	  fhnALaEmbPC->Fill(vALaPCEmb);
	  
	}//end: embedded jets ALa background via perpendicular cones
	
	
	
	jetConeALaEmblist->Clear();
	jetPerpConeALaEmblist->Clear();
	
	
	//END EMBEDDED V0 particles
	
}

//_____________________________________________________________________
