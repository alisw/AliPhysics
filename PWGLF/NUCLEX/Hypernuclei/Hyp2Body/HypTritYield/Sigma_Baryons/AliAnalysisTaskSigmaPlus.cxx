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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Sigma Plus Class                                                      //
//                                                                        //
//  This class is used for the reconstruction of Sigma+ baryons in the    //
//  Sigma+ -> p+ + pi0 decay channel.                                     //
//                                                                        //
//  Author: B.Heybeck (b.heybeck@cern.ch)                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include <TRandom3.h>
#include "TPDGCode.h"
#include <TDatabasePDG.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliEventPoolManager.h"
#include "AliV0vertexer.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliVertexerTracks.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "TChain.h"
#include "AliStack.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisTaskSigmaPlus.h"
#include "AliAODVZERO.h"
#include "AliKFParticleBase.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliAnalysisUtils.h"

#include <thread>         // std::this_thread::sleep_for()
#include <chrono>         // std::chrono::seconds()

class TTree;
class TParticle;
class TVector3;
class AliAODVertex;
class AliAODv0;

class AliAnalysisTaskSigmaPlus;    // analysis class

using std::cout;            

// ClassImp: necessary for root
ClassImp(AliAnalysisTaskSigmaPlus)
ClassImp(AliAODTrackreduced)
ClassImp(AliAODTrackcorrelation)
ClassImp(AliAODClusterreduced)

AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus() : AliAnalysisTaskSE(), 
fOutputList(0), aodEvent(0x0), mcEvent(0x0), AODMCTrackArray(0x0), 
fPIDResponse(0), fEvPoolMgr(0x0), fEvPoolMgr2(0x0), fEvPoolMgr3(0x0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaPairTreeSE(0x0), fSigmaPairTreeME(0x0), fSigmaMEBackgroundTree(0x0),
fSigmaPHOSCandTree(0x0), fSigmaPairTreePHOSSE(0x0), fSigmaPairTreePHOSME(0x0), fSigmaPHOSMEBkgTree(0x0), 
fSigmaCalcCandTree(0x0), fProtonTree(0x0), fThetaFunc(0x0), fPhiFunc(0x0),
cElectronMass(0), cProtonMass(0), cSigmaMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(-999), primaryVtxPosY(-999), primaryVtxPosZ(-999),
primaryVtxPosXMC(-999), primaryVtxPosYMC(-999), primaryVtxPosZMC(-999), 
nTracks(0), Centrality(0), EventTriggers(0),
fRefMultComb05(0), fRefMultComb08(0), fRefMultComb10(0), fGlobalEventID(0),
fEventhasSigma(kFALSE), fEventhasProton(kFALSE), fEventhasSigmaCand(kFALSE), fDebug(kFALSE), fRemoveGenPileup(kTRUE),

fProcessMCParticles(kTRUE), 

fProcessProtons(kTRUE),
fSaveProtons(kTRUE), 
fProcessV0s(kTRUE), 
fProcessReco(kFALSE), 
fSavePartCand(kFALSE),
fFillPairTreeSE(kFALSE),
fFillPairTreeME(kFALSE),
fSaveMixedBackground(kFALSE),

fProcessClusters(kTRUE),
fMapClusterstoTracks(kTRUE),
fProcessRecoPHOS(kTRUE),
fProcessRecoCalc(kFALSE),
fSavePartCandPHOS(kTRUE),
fSavePartCandCalc(kFALSE),
fFillPHOSPairTreeSE(kFALSE),
fFillPHOSPairTreeME(kFALSE),
fSavePHOSMixedBackground(kTRUE),
fSaveAdditionalBranches(kFALSE),
fSaveMCBranches(kTRUE),
fSaveAddMCBranches(kTRUE),

fMaxVertexZ(10),

fEvPoolSize(10),
fEvTrackSize(10000),
fCentralityBins(29),
fMinCentBin(-5),
fMaxCentBin(140),
fZvtxBins(10),
fMinZBin(0),
fMaxZBin(10),

fEvPoolSize2(20),
fEvTrackSize2(10000),
fCentralityBins2(29),
fMinCentBin2(-5),
fMaxCentBin2(140),
fZvtxBins2(10),
fMinZBin2(0),
fMaxZBin2(10),

fRequireSigma(kFALSE),
fRequireProton(kFALSE),
fRequireSigmaCand(kFALSE),
fUseAbsZ(kFALSE),
fUseAbsZCorr(kTRUE),

fRejectNegIDs(kTRUE),
fRejectZeroFilterBit(kFALSE),
fMaxProtEta(1),    
fMinTPCClustProt(40),
fMaxNsigProtTPC(4),
fRequireProtonTPC(kTRUE),
fRequireProtonTOF(kFALSE),
fRequireProtonTOFforPairs(kTRUE),
fMaxNsigProtTOF(6),  
fMaxpOnlyTPCPID(0.9),
fMinProtpt(0), 
fMaxProtpt(15),   

fStrictMaxProtEta(0.9),
fStrictMinTPCClustProt(60),
fStrictMaxNsigProtTPC(3),
fStrictMaxNsigProtTOF(5),
fStrictMaxpOnlyTPCPID(0.8),
fStrictMinProtpt(0),
fStrictMaxProtpt(5),

fMaxMCEta(0.9),

fMaxDaughtEta(1.2),
fMinTPCClustDaught(20),
fMaxNsigDaughtTPC(6),
fMaxalpha(1.1),
fMaxqt(0.04),
fMaxopenangle(0.4),
fMaxdeltatheta(0.15),
fMinV0CPA(0.8),
fMinV0Radius(1), 
fMaxV0Radius(250),
fMaxphotonmass(0.08),

fRequirePHOS(kFALSE),
fMinClusterBeta(0),
fMinClusterDy(8),    
fMaxClusterM02(4),

fCleanAutoCorr(kTRUE),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.05),
fMaxSigmaY(0.9),  
fMaxSigmaMass(1.4),
fMinProtonDCAxy(0.005),
fMinProtonDCAz(-1),
fMaxProtonDCAxy(5),
fMaxProtonDCAz(9999),
fRequireDCACut(kFALSE),

fMinPi0MassPHOS(0.09), 
fMaxPi0MassPHOS(0.16),  
fMaxSigmaPAPHOS(0.02),
fMaxSigmaPAPHOSHM(0.02),
fMinSigmaAntiPAPHOS(0.005),
fMaxProtPhotDCA(2),
fMinSigmaDCAtoPVPHOS(0.5),
fMaxSigmaDCAtoPVPHOS(60),
fMaxSigmaYPHOS(0.9),  
fMaxSigmaMassPHOS(1.35),
fMinProtonDCAxyPHOS(0.005),
fMinProtonDCAzPHOS(-1),
fMaxProtonDCAxyPHOS(5),
fMaxProtonDCAzPHOS(9999),
fRequireDCACutPHOS(kFALSE),

fMinCorrPi0Mass(0.1),    
fMaxCorrPi0Mass(0.16),    
fMaxCorrSigmaPA(0.06),     
fMinCorrSigmaMass(1.13),   
fMaxCorrSigmaMass(1.25),   
fMinCorrProtonDCAxy(0.005),
fMaxCorrPairProtonDCAxy(0.2),
fMaxCorrPairProtonDCAz(0.4),
fMaxCorrkstar(600),

fMinCorrPi0MassPHOS(0.09),
fMaxCorrPi0MassPHOS(0.16),
fMaxCorrSigmaPAPHOS(0.02),

fCheckProtonV0IDs(kTRUE),
fMaxSigmaMassCalc(1.35),
fMinSigmaAntiPACalc(0.005),
fMaxSigmaAntiPACalc(0.6),
fMaxIter(10),
fRejectNegPar1(kTRUE),
fMaxSigmaYCalc(0.9),
fThetaPar0(0.000133299),
fThetaPar1(0.0016761),
fThetaRange(0.02),
fPhiPar0(0.000129845),
fPhiPar1(0.00199688),
fPhiRange(0.02),

fMaxNsigPionTPC(2),
fMaxNsigKaonTPC(2),
fMaxNsigElecTPC(2),
fMaxNsigPionTOF(2),
fMaxNsigKaonTOF(2),
fMaxNsigElecTOF(2),
fRejOtherTOF(kTRUE),

fIsMCSigma(kFALSE),
fIsMCDalitz(kFALSE),
fIsMCPrimary(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fHas4DiffIDs(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(0),
fSigMCLabel(0),
fSigProtonID(-999),
fSigProtonStatus(0),
fSigProtonFilterMap(0),
fSigEventID(0),
fSigCentrality(-999),
fSigRefMultComb05(-999),
fSigRefMultComb08(-999),
fSigRefMultComb10(-999),
fSigBField(-999),
fInvSigMass(-999),
fInvSigpropMass(-999),
fInvSigMassUncorr(-999),
fSigY(-999),               
fSigYprop(-999),               
fSigPA(-999),
fSigPAprop(-999),
fSigAntiPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fSigPt(-999),        
fSigPxprop(-999),
fSigPyprop(-999),
fSigPzprop(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
fSigFlightDist(-999),        
fSigDecayVertXMC(-999),
fSigDecayVertYMC(-999),
fSigDecayVertZMC(-999),
fSigPxMC(-999),        
fSigPyMC(-999),        
fSigPzMC(-999),
fPrimVertXMC(-999),
fPrimVertYMC(-999),
fPrimVertZMC(-999),
fPhoton1Px(-999),
fPhoton1Py(-999),
fPhoton1Pz(-999),
fPhoton2Px(-999),
fPhoton2Py(-999),
fPhoton2Pz(-999),
fPhotonDaughtMaxEta(-999),
fPhotonsMaxDeltaTheta(-999),
fPhoton1CPA(-999),
fPhoton2CPA(-999),        
fPhoton1Radius(-999),
fPhoton2Radius(-999),
fPhoton1DCAPV(-999),
fPhoton2DCAPV(-999),
fPhoton1DCASV(-999),
fPhoton2DCASV(-999),
fTrackDCASV(-999),
fTrackDCASVKF(-999),
fKFChi2(-999),
fPhotonsMinCluster(-999),
fPhotonsMinITSCluster(-999),
fPhotonsMaxalpha(-999),
fPhotonsMaxqt(-999),
fPhotonsMaxOpenAngle(-999),
fPhotonsMaxinvmass(-999),
fPhotonsMaxNSigTPC(-999),
fPhotonsMaxChi2(-999),
fInvPi0Mass(-999),
fPi0Px(-999),
fPi0Py(-999),
fPi0Pz(-999),
fPi0DecayVertX(-999),
fPi0DecayVertY(-999),
fPi0DecayVertZ(-999),
fPi0PhotPhotDCA(-999),
fProtonPx(-999),
fProtonPy(-999),
fProtonPz(-999),
fProtonX(-999),
fProtonY(-999),
fProtonZ(-999),
fProtonEta(-999),
fProtonpropPx(-999),
fProtonpropPy(-999),
fProtonpropPz(-999),
fProtonDCAtoPVxy(-999),
fProtonDCAtoPVz(-999),
fProtonPi0DCA(-999),
fProtonNSigTPC(-999),
fProtonNSigTOF(-999),
fProtonNCluster(-999),
fProtonNITSCluster(-999),
fProtonChi2(-999),  
fProtonNSigTPCPion(-999),
fProtonNSigTPCKaon(-999),
fProtonNSigTPCElec(-999),
fProtonNSigTOFPion(-999),
fProtonNSigTOFKaon(-999),
fProtonNSigTOFElec(-999),

fTrackLabel(-999),
fTrackPDGCode(-999),
fTrackMotherID(-999),
fTrackMotherPDGCode(-999),
fTrackMotherMCPx(-999),
fTrackMotherMCPy(-999),
fTrackMotherMCPz(-999),
fTrackMCPx(-999),
fTrackMCPy(-999),
fTrackMCPz(-999),
fPhoton1Label(-999),
fPhoton1PDGCode(-999),
fPhoton1MotherID(-999),
fPhoton1MotherPDGCode(-999),
fPhoton1GMotherID(-999),
fPhoton1GMotherPDGCode(-999),
fPhoton1MCPx(-999),
fPhoton1MCPy(-999),
fPhoton1MCPz(-999),
fPhoton2Label(-999),
fPhoton2PDGCode(-999),
fPhoton2MotherID(-999),
fPhoton2MotherPDGCode(-999),
fPhoton2GMotherID(-999),
fPhoton2GMotherPDGCode(-999),
fPhoton2MCPx(-999),
fPhoton2MCPy(-999),
fPhoton2MCPz(-999),

fIsClusterEMCAL(kFALSE),
fIsClusterPHOS(kFALSE), 
fCaloPhotonX(-999),
fCaloPhotonY(-999),
fCaloPhotonZ(-999),
fConvPhotonX(-999),
fConvPhotonY(-999),
fConvPhotonZ(-999),
fConvPhotonSecPA(-999),
fCaloPhotonPxMC(-999),
fCaloPhotonPyMC(-999),
fCaloPhotonPzMC(-999),
fCaloPhotonE(-999),
fCaloPhotonELead(-999),
fCaloPhotonEcorr(-999),
fCaloPhotonEMC(-999),
fClustNLabels(-999),
fClustPDGCode(-999),
fClustDispersion(-999),
fClustM20(-999),
fClustM02(-999),
fClustNTracksMatched(-999),
fClustTrackDx(-999),
fClustTrackDz(-999),
fClustTrackD(-999),
fClustTOF(-999),
fClustBeta(-999),
fClustNCells(-999),
fClustDisttoBC(-999),

fPairProtonIsMC(kFALSE),
fPairProtonPDG(-999),
fPairProtonIsPrimary(kFALSE),
fPairProtonPx(-999),
fPairProtonPy(-999),
fPairProtonPz(-999),
fPairProtonPxMC(-999),
fPairProtonPyMC(-999),
fPairProtonPzMC(-999),
fPairProtonP(-999),
fPairProtonEta(-999),
fPairProtonCharge(-999),
fPairProtonDCAtoPVxy(-999),
fPairProtonDCAtoPVz(-999),       
fPairProtonNSigTPC(-999),      
fPairProtonNSigTOF(-999),
fPairProtNSigTPCPion(-999),        
fPairProtNSigTPCKaon(-999),
fPairProtNSigTPCElec(-999),
fPairProtNSigTOFPion(-999),        
fPairProtNSigTOFKaon(-999),
fPairProtNSigTOFElec(-999),
fPairProtonChi2(-999),
fPairProtonCluster(-999),
fPairProtonITSCluster(-999),
fPairProtonID(-999),
fPairProtonStatus(0),
fPairProtonFilterMap(0),
fNIter(-999),
fDeltaPhi(-999),
fDeltaTheta(-999),
fSigmaProtonkstar(-999),
fSigmaProtonpropkstar(-999)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus(const char* name) : AliAnalysisTaskSE(name),
fOutputList(0), aodEvent(0x0), mcEvent(0x0), AODMCTrackArray(0x0), 
fPIDResponse(0), fEvPoolMgr(0x0), fEvPoolMgr2(0x0), fEvPoolMgr3(0x0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaPairTreeSE(0x0), fSigmaPairTreeME(0x0), fSigmaMEBackgroundTree(0x0),
fSigmaPHOSCandTree(0x0), fSigmaPairTreePHOSSE(0x0), fSigmaPairTreePHOSME(0x0), fSigmaPHOSMEBkgTree(0x0), 
fSigmaCalcCandTree(0x0), fProtonTree(0x0), fThetaFunc(0x0), fPhiFunc(0x0),
cElectronMass(0), cProtonMass(0), cSigmaMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(-999), primaryVtxPosY(-999), primaryVtxPosZ(-999),
primaryVtxPosXMC(-999), primaryVtxPosYMC(-999), primaryVtxPosZMC(-999), 
nTracks(0), Centrality(0), EventTriggers(0),
fRefMultComb05(0), fRefMultComb08(0), fRefMultComb10(0), fGlobalEventID(0),
fEventhasSigma(kFALSE), fEventhasProton(kFALSE), fEventhasSigmaCand(kFALSE), fDebug(kFALSE), fRemoveGenPileup(kTRUE),

fProcessMCParticles(kTRUE), 

fProcessProtons(kTRUE),
fSaveProtons(kTRUE), 
fProcessV0s(kTRUE), 
fProcessReco(kFALSE), 
fSavePartCand(kFALSE),
fFillPairTreeSE(kFALSE),
fFillPairTreeME(kFALSE),
fSaveMixedBackground(kFALSE),

fProcessClusters(kTRUE),
fMapClusterstoTracks(kTRUE),
fProcessRecoPHOS(kTRUE),
fProcessRecoCalc(kFALSE),
fSavePartCandPHOS(kTRUE),
fSavePartCandCalc(kFALSE),
fFillPHOSPairTreeSE(kFALSE),
fFillPHOSPairTreeME(kFALSE),
fSavePHOSMixedBackground(kTRUE),
fSaveAdditionalBranches(kFALSE),
fSaveMCBranches(kTRUE),
fSaveAddMCBranches(kTRUE),

fMaxVertexZ(10),

fEvPoolSize(10),
fEvTrackSize(10000),
fCentralityBins(29),
fMinCentBin(-5),
fMaxCentBin(140),
fZvtxBins(10),
fMinZBin(0),
fMaxZBin(10),

fEvPoolSize2(20),
fEvTrackSize2(10000),
fCentralityBins2(29),
fMinCentBin2(-5),
fMaxCentBin2(140),
fZvtxBins2(10),
fMinZBin2(0),
fMaxZBin2(10),

fRequireSigma(kFALSE),
fRequireProton(kFALSE),
fRequireSigmaCand(kFALSE),
fUseAbsZ(kFALSE),
fUseAbsZCorr(kTRUE),

fRejectNegIDs(kTRUE),
fRejectZeroFilterBit(kFALSE),
fMaxProtEta(1),    
fMinTPCClustProt(40),
fMaxNsigProtTPC(4),
fRequireProtonTPC(kTRUE),
fRequireProtonTOF(kFALSE),
fRequireProtonTOFforPairs(kTRUE),
fMaxNsigProtTOF(6),  
fMaxpOnlyTPCPID(0.9),
fMinProtpt(0), 
fMaxProtpt(15),   

fStrictMaxProtEta(0.9),
fStrictMinTPCClustProt(60),
fStrictMaxNsigProtTPC(3),
fStrictMaxNsigProtTOF(5),
fStrictMaxpOnlyTPCPID(0.8),
fStrictMinProtpt(0),
fStrictMaxProtpt(5),

fMaxMCEta(0.9),

fMaxDaughtEta(1.2),
fMinTPCClustDaught(20),
fMaxNsigDaughtTPC(6),
fMaxalpha(1.1),
fMaxqt(0.04),
fMaxopenangle(0.4),
fMaxdeltatheta(0.15),
fMinV0CPA(0.8),
fMinV0Radius(1), 
fMaxV0Radius(250),
fMaxphotonmass(0.08),

fRequirePHOS(kFALSE),
fMinClusterBeta(0),
fMinClusterDy(8),    
fMaxClusterM02(4),

fCleanAutoCorr(kTRUE),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.05),
fMaxSigmaY(0.9),  
fMaxSigmaMass(1.4),
fMinProtonDCAxy(0.005),
fMinProtonDCAz(-1),
fMaxProtonDCAxy(5),
fMaxProtonDCAz(9999),
fRequireDCACut(kFALSE),

fMinPi0MassPHOS(0.09), 
fMaxPi0MassPHOS(0.16),  
fMaxSigmaPAPHOS(0.02),
fMaxSigmaPAPHOSHM(0.02),
fMinSigmaAntiPAPHOS(0.005),
fMaxProtPhotDCA(2),
fMinSigmaDCAtoPVPHOS(0.5),
fMaxSigmaDCAtoPVPHOS(60),
fMaxSigmaYPHOS(0.9),  
fMaxSigmaMassPHOS(1.35),
fMinProtonDCAxyPHOS(0.005),
fMinProtonDCAzPHOS(-1),
fMaxProtonDCAxyPHOS(5),
fMaxProtonDCAzPHOS(9999),
fRequireDCACutPHOS(kFALSE),

fMinCorrPi0Mass(0.1),    
fMaxCorrPi0Mass(0.16),    
fMaxCorrSigmaPA(0.06),     
fMinCorrSigmaMass(1.13),   
fMaxCorrSigmaMass(1.25),   
fMinCorrProtonDCAxy(0.005),
fMaxCorrPairProtonDCAxy(0.2),
fMaxCorrPairProtonDCAz(0.4),
fMaxCorrkstar(600),

fMinCorrPi0MassPHOS(0.09),
fMaxCorrPi0MassPHOS(0.16),
fMaxCorrSigmaPAPHOS(0.02),

fCheckProtonV0IDs(kTRUE),
fMaxSigmaMassCalc(1.35),
fMinSigmaAntiPACalc(0.005),
fMaxSigmaAntiPACalc(0.6),
fMaxIter(10),
fRejectNegPar1(kTRUE),
fMaxSigmaYCalc(0.9),
fThetaPar0(0.000133299),
fThetaPar1(0.0016761),
fThetaRange(0.02),
fPhiPar0(0.000129845),
fPhiPar1(0.00199688),
fPhiRange(0.02),

fMaxNsigPionTPC(2),
fMaxNsigKaonTPC(2),
fMaxNsigElecTPC(2),
fMaxNsigPionTOF(2),
fMaxNsigKaonTOF(2),
fMaxNsigElecTOF(2),
fRejOtherTOF(kTRUE),

fIsMCSigma(kFALSE),
fIsMCDalitz(kFALSE),
fIsMCPrimary(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fHas4DiffIDs(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(0),
fSigMCLabel(0),
fSigProtonID(-999),
fSigProtonStatus(0),
fSigProtonFilterMap(0),
fSigEventID(0),
fSigCentrality(-999),
fSigRefMultComb05(-999),
fSigRefMultComb08(-999),
fSigRefMultComb10(-999),
fSigBField(-999),
fInvSigMass(-999),
fInvSigpropMass(-999),
fInvSigMassUncorr(-999),
fSigY(-999),               
fSigYprop(-999),               
fSigPA(-999),
fSigPAprop(-999),
fSigAntiPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fSigPt(-999),        
fSigPxprop(-999),
fSigPyprop(-999),
fSigPzprop(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
fSigFlightDist(-999),        
fSigDecayVertXMC(-999),
fSigDecayVertYMC(-999),
fSigDecayVertZMC(-999),
fSigPxMC(-999),        
fSigPyMC(-999),        
fSigPzMC(-999),
fPrimVertXMC(-999),
fPrimVertYMC(-999),
fPrimVertZMC(-999),
fPhoton1Px(-999),
fPhoton1Py(-999),
fPhoton1Pz(-999),
fPhoton2Px(-999),
fPhoton2Py(-999),
fPhoton2Pz(-999),
fPhotonDaughtMaxEta(-999),
fPhotonsMaxDeltaTheta(-999),
fPhoton1CPA(-999),
fPhoton2CPA(-999),        
fPhoton1Radius(-999),
fPhoton2Radius(-999),
fPhoton1DCAPV(-999),
fPhoton2DCAPV(-999),
fPhoton1DCASV(-999),
fPhoton2DCASV(-999),
fTrackDCASV(-999),
fTrackDCASVKF(-999),
fKFChi2(-999),
fPhotonsMinCluster(-999),
fPhotonsMinITSCluster(-999),
fPhotonsMaxalpha(-999),
fPhotonsMaxqt(-999),
fPhotonsMaxOpenAngle(-999),
fPhotonsMaxinvmass(-999),
fPhotonsMaxNSigTPC(-999),
fPhotonsMaxChi2(-999),
fInvPi0Mass(-999),
fPi0Px(-999),
fPi0Py(-999),
fPi0Pz(-999),
fPi0DecayVertX(-999),
fPi0DecayVertY(-999),
fPi0DecayVertZ(-999),
fPi0PhotPhotDCA(-999),
fProtonPx(-999),
fProtonPy(-999),
fProtonPz(-999),
fProtonX(-999),
fProtonY(-999),
fProtonZ(-999),
fProtonEta(-999),
fProtonpropPx(-999),
fProtonpropPy(-999),
fProtonpropPz(-999),
fProtonDCAtoPVxy(-999),
fProtonDCAtoPVz(-999),
fProtonPi0DCA(-999),
fProtonNSigTPC(-999),
fProtonNSigTOF(-999),
fProtonNCluster(-999),
fProtonNITSCluster(-999),
fProtonChi2(-999),  
fProtonNSigTPCPion(-999),
fProtonNSigTPCKaon(-999),
fProtonNSigTPCElec(-999),
fProtonNSigTOFPion(-999),
fProtonNSigTOFKaon(-999),
fProtonNSigTOFElec(-999),

fTrackLabel(-999),
fTrackPDGCode(-999),
fTrackMotherID(-999),
fTrackMotherPDGCode(-999),
fTrackMotherMCPx(-999),
fTrackMotherMCPy(-999),
fTrackMotherMCPz(-999),
fTrackMCPx(-999),
fTrackMCPy(-999),
fTrackMCPz(-999),
fPhoton1Label(-999),
fPhoton1PDGCode(-999),
fPhoton1MotherID(-999),
fPhoton1MotherPDGCode(-999),
fPhoton1GMotherID(-999),
fPhoton1GMotherPDGCode(-999),
fPhoton1MCPx(-999),
fPhoton1MCPy(-999),
fPhoton1MCPz(-999),
fPhoton2Label(-999),
fPhoton2PDGCode(-999),
fPhoton2MotherID(-999),
fPhoton2MotherPDGCode(-999),
fPhoton2GMotherID(-999),
fPhoton2GMotherPDGCode(-999),
fPhoton2MCPx(-999),
fPhoton2MCPy(-999),
fPhoton2MCPz(-999),

fIsClusterEMCAL(kFALSE),
fIsClusterPHOS(kFALSE), 
fCaloPhotonX(-999),
fCaloPhotonY(-999),
fCaloPhotonZ(-999),
fConvPhotonX(-999),
fConvPhotonY(-999),
fConvPhotonZ(-999),
fConvPhotonSecPA(-999),
fCaloPhotonPxMC(-999),
fCaloPhotonPyMC(-999),
fCaloPhotonPzMC(-999),
fCaloPhotonE(-999),
fCaloPhotonELead(-999),
fCaloPhotonEcorr(-999),
fCaloPhotonEMC(-999),
fClustNLabels(-999),
fClustPDGCode(-999),
fClustDispersion(-999),
fClustM20(-999),
fClustM02(-999),
fClustNTracksMatched(-999),
fClustTrackDx(-999),
fClustTrackDz(-999),
fClustTrackD(-999),
fClustTOF(-999),
fClustBeta(-999),
fClustNCells(-999),
fClustDisttoBC(-999),

fPairProtonIsMC(kFALSE),
fPairProtonPDG(-999),
fPairProtonIsPrimary(kFALSE),
fPairProtonPx(-999),
fPairProtonPy(-999),
fPairProtonPz(-999),
fPairProtonPxMC(-999),
fPairProtonPyMC(-999),
fPairProtonPzMC(-999),
fPairProtonP(-999),
fPairProtonEta(-999),
fPairProtonCharge(-999),
fPairProtonDCAtoPVxy(-999),
fPairProtonDCAtoPVz(-999),       
fPairProtonNSigTPC(-999),      
fPairProtonNSigTOF(-999),
fPairProtNSigTPCPion(-999),        
fPairProtNSigTPCKaon(-999),
fPairProtNSigTPCElec(-999),
fPairProtNSigTOFPion(-999),        
fPairProtNSigTOFKaon(-999),
fPairProtNSigTOFElec(-999),
fPairProtonChi2(-999),
fPairProtonCluster(-999),
fPairProtonITSCluster(-999),
fPairProtonID(-999),
fPairProtonStatus(0),
fPairProtonFilterMap(0),
fNIter(-999),
fDeltaPhi(-999),
fDeltaTheta(-999),
fSigmaProtonkstar(-999),
fSigmaProtonpropkstar(-999)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!

    DefineOutput(2, TTree::Class());    //Tree with Sigma candidates
    DefineOutput(3, TTree::Class());    //Tree with Sigma Proton Pairs Same Event
    DefineOutput(4, TTree::Class());    //Tree with Sigma Proton Pairs Mixed Event
    DefineOutput(5, TTree::Class());    //Tree with Sigma candidate Mixed Event Background
    DefineOutput(6, TTree::Class());    //Tree with Sigma candidates from PHOS
    DefineOutput(7, TTree::Class());    //Tree with Sigma Proton Pairs Same Event
    DefineOutput(8, TTree::Class());    //Tree with Sigma Proton Pairs Mixed Event
    DefineOutput(9, TTree::Class());    //Tree with Sigma candidate Mixed Event Background
    DefineOutput(10,TTree::Class());    //Tree with Sigma candidates calculated
    DefineOutput(11,TTree::Class());    //Tree Protons
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::~AliAnalysisTaskSigmaPlus()
{
    // destructor
    // delete objects from memory at the end of the task
    //Histograms
    if(fOutputList) delete fOutputList;     
    //Trees
    if(fSigmaCandTree) delete fSigmaCandTree;
    if(fSigmaPairTreeSE) delete fSigmaPairTreeSE;
    if(fSigmaPairTreeME) delete fSigmaPairTreeME;
    if(fSigmaMEBackgroundTree) delete fSigmaMEBackgroundTree;
    if(fSigmaPHOSCandTree) delete fSigmaPHOSCandTree;
    if(fSigmaPairTreePHOSSE) delete fSigmaPairTreePHOSSE;
    if(fSigmaPairTreePHOSME) delete fSigmaPairTreePHOSME;
    if(fSigmaPHOSMEBkgTree) delete fSigmaPHOSMEBkgTree;
    if(fSigmaCalcCandTree) delete fSigmaCalcCandTree;
    if(fProtonTree) delete fProtonTree;
    //Functions
    if(fThetaFunc) delete fThetaFunc;
    if(fPhiFunc) delete fPhiFunc;
    //Pool managers
    if(fEvPoolMgr) delete fEvPoolMgr;
    if(fEvPoolMgr2) delete fEvPoolMgr2;
    if(fEvPoolMgr3) delete fEvPoolMgr3;
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlus::UserCreateOutputObjects()
{ 

    //Check analysis settings. Not all combinations make sense
    if(!fProcessProtons||!fProcessV0s) fProcessReco = kFALSE;
    if(!fProcessReco) {fSavePartCand = kFALSE; fFillPairTreeSE = kFALSE; fFillPairTreeME = kFALSE; fSaveMixedBackground = kFALSE;}
    if(!fSavePartCand) fSaveMixedBackground = kFALSE;
    if(!fFillPairTreeSE) fFillPairTreeME = kFALSE;
    if(!fProcessProtons||!fProcessV0s||!fProcessClusters) fProcessRecoPHOS = kFALSE;
    if(!fProcessRecoPHOS) {fSavePartCandPHOS = kFALSE; fFillPHOSPairTreeSE = kFALSE; fFillPHOSPairTreeME = kFALSE; fSavePHOSMixedBackground = kFALSE;}
    if(!fSavePartCandPHOS) fSavePHOSMixedBackground = kFALSE;
    if(!fFillPHOSPairTreeSE) fFillPHOSPairTreeME = kFALSE;
    if(!fProcessProtons||!fProcessV0s) fProcessRecoCalc = kFALSE;

    // Create output objects. Called once at start of the analysis (RUNTIME). 
    
    fOutputList = new TList();          // List which contains all Histograms.                           
    fOutputList->SetOwner(kTRUE);       // The list is owner of all objects it contains and will delete them if requested.

    // PID Response
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();               //Get Analysis Manager
      if(man)
      {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());  //Get Input Handler
        if(inputHandler) fPIDResponse = inputHandler->GetPIDResponse();    //Retrieve the AliPIDResponse object from the analysis manager
        else AliWarning("No Input Handler!");
      }
      else AliWarning("No Analysis Manager!");

    const Int_t nCentralityBins = fCentralityBins;
    const Int_t nZvtxBins = fZvtxBins;

    if(fSaveMixedBackground){
      //Create Pool Manager for Event Mixing
      Int_t poolSize = fEvPoolSize;
      Int_t trackDepth = fEvTrackSize;
      Double_t centBinz[nCentralityBins+1];
      for(Int_t b=0;b<=nCentralityBins;b++){centBinz[b]=fMinCentBin+b*(fMaxCentBin-fMinCentBin)/nCentralityBins;}
      Double_t vertexBinz[nZvtxBins+1];
      for(Int_t b=0;b<=nZvtxBins;b++){vertexBinz[b]=fMinZBin+b*(fMaxZBin-fMinZBin)/nZvtxBins;}
      Double_t* centBins = centBinz; 
      Double_t* vertexBins = vertexBinz;
      fEvPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);
    }

    if(fSavePHOSMixedBackground){
      //Create Pool Manager for Event Mixing
      Int_t poolSize = fEvPoolSize;
      Int_t trackDepth = fEvTrackSize;
      Double_t centBinz[nCentralityBins+1];
      for(Int_t b=0;b<=nCentralityBins;b++){centBinz[b]=fMinCentBin+b*(fMaxCentBin-fMinCentBin)/nCentralityBins;}
      Double_t vertexBinz[nZvtxBins+1];
      for(Int_t b=0;b<=nZvtxBins;b++){vertexBinz[b]=fMinZBin+b*(fMaxZBin-fMinZBin)/nZvtxBins;}
      Double_t* centBins = centBinz; 
      Double_t* vertexBins = vertexBinz;
      fEvPoolMgr3 = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);
    }

    const Int_t nCentralityBins2 = fCentralityBins2;
    const Int_t nZvtxBins2 = fZvtxBins2;

    if(fFillPairTreeME||fFillPHOSPairTreeME){
      //Create Pool Manager for Event Mixing
      Int_t poolSize2 = fEvPoolSize2;
      Int_t trackDepth2 = fEvTrackSize2;
      Double_t centBinz2[nCentralityBins2+1];
      for(Int_t b=0;b<=nCentralityBins2;b++){centBinz2[b]=fMinCentBin2+b*(fMaxCentBin2-fMinCentBin2)/nCentralityBins2;}
      Double_t vertexBinz2[nZvtxBins2+1];
      for(Int_t b=0;b<=nZvtxBins2;b++){vertexBinz2[b]=fMinZBin2+b*(fMaxZBin2-fMinZBin2)/nZvtxBins2;}
      Double_t* centBins2 = centBinz2; 
      Double_t* vertexBins2 = vertexBinz2;
      fEvPoolMgr2 = new AliEventPoolManager(poolSize2, trackDepth2, nCentralityBins2, centBins2, nZvtxBins2, vertexBins2);
    }

    if(fProcessRecoCalc){
      fThetaFunc = new TF1("fThetaFunc",Form("%.9f/(abs(x)+%.9f)",fThetaPar0,fThetaPar1),-fThetaRange,fThetaRange);
      fPhiFunc = new TF1("fPhiFunc",Form("%.9f/(abs(x)+%.9f)",fPhiPar0,fPhiPar1),-fPhiRange,fPhiRange);
    }

    // Create TTree of Sigma Candidates
    fSigmaCandTree = new TTree("fSigmaCandTree","Tree of Sigma Candidates");
    fSigmaCandTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaCandTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaCandTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaCandTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaCandTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaCandTree->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaCandTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaCandTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCandTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaCandTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCandTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaCandTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaCandTree->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaCandTree->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaCandTree->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaCandTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaCandTree->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaCandTree->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaCandTree->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaCandTree->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaCandTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCandTree->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaCandTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaCandTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    //fSigmaCandTree->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
    fSigmaCandTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaCandTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaCandTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaCandTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaCandTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaCandTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaCandTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaCandTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaCandTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaCandTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaCandTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaCandTree->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaCandTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaCandTree->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaCandTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCandTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCandTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaCandTree->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaCandTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaCandTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaCandTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaCandTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaCandTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaCandTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    if(fSaveMCBranches){
      fSigmaCandTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaCandTree->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaCandTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaCandTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaCandTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaCandTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaCandTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaCandTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaCandTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaCandTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaCandTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaCandTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaCandTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    }
    if(fSaveAddMCBranches){
      fSigmaCandTree->Branch("fTrackLabel",&fTrackLabel,"fTrackLabel/I");
      fSigmaCandTree->Branch("fTrackPDGCode",&fTrackPDGCode,"fTrackPDGCode/I");
      fSigmaCandTree->Branch("fTrackMotherID",&fTrackMotherID,"fTrackMotherID/I");
      fSigmaCandTree->Branch("fTrackMotherPDGCode",&fTrackMotherPDGCode,"fTrackMotherPDGCode/I");
      fSigmaCandTree->Branch("fTrackMotherMCPx",&fTrackMotherMCPx,"fTrackMotherMCPx/F");
      fSigmaCandTree->Branch("fTrackMotherMCPy",&fTrackMotherMCPy,"fTrackMotherMCPy/F");
      fSigmaCandTree->Branch("fTrackMotherMCPz",&fTrackMotherMCPz,"fTrackMotherMCPz/F");
      fSigmaCandTree->Branch("fTrackMCPx",&fTrackMCPx,"fTrackMCPx/F");
      fSigmaCandTree->Branch("fTrackMCPy",&fTrackMCPy,"fTrackMCPy/F");
      fSigmaCandTree->Branch("fTrackMCPz",&fTrackMCPz,"fTrackMCPz/F");
      fSigmaCandTree->Branch("fPhoton1Label",&fPhoton1Label,"fPhoton1Label/I");
      fSigmaCandTree->Branch("fPhoton1PDGCode",&fPhoton1PDGCode,"fPhoton1PDGCode/I");
      fSigmaCandTree->Branch("fPhoton1MotherID",&fPhoton1MotherID,"fPhoton1MotherID/I");
      fSigmaCandTree->Branch("fPhoton1MotherPDGCode",&fPhoton1MotherPDGCode,"fPhoton1MotherPDGCode/I");
      fSigmaCandTree->Branch("fPhoton1GMotherID",&fPhoton1GMotherID,"fPhoton1GMotherID/I");
      fSigmaCandTree->Branch("fPhoton1GMotherPDGCode",&fPhoton1GMotherPDGCode,"fPhoton1GMotherPDGCode/I");
      fSigmaCandTree->Branch("fPhoton1MCPx",&fPhoton1MCPx,"fPhoton1MCPx/F");
      fSigmaCandTree->Branch("fPhoton1MCPy",&fPhoton1MCPy,"fPhoton1MCPy/F");
      fSigmaCandTree->Branch("fPhoton1MCPz",&fPhoton1MCPz,"fPhoton1MCPz/F");
      fSigmaCandTree->Branch("fPhoton2Label",&fPhoton2Label,"fPhoton2Label/I");
      fSigmaCandTree->Branch("fPhoton2PDGCode",&fPhoton2PDGCode,"fPhoton2PDGCode/I");
      fSigmaCandTree->Branch("fPhoton2MotherID",&fPhoton2MotherID,"fPhoton2MotherID/I");
      fSigmaCandTree->Branch("fPhoton2MotherPDGCode",&fPhoton2MotherPDGCode,"fPhoton2MotherPDGCode/I");
      fSigmaCandTree->Branch("fPhoton2GMotherID",&fPhoton2GMotherID,"fPhoton2GMotherID/I");
      fSigmaCandTree->Branch("fPhoton2GMotherPDGCode",&fPhoton2GMotherPDGCode,"fPhoton2GMotherPDGCode/I");
      fSigmaCandTree->Branch("fPhoton2MCPx",&fPhoton2MCPx,"fPhoton2MCPx/F");
      fSigmaCandTree->Branch("fPhoton2MCPy",&fPhoton2MCPy,"fPhoton2MCPy/F");
      fSigmaCandTree->Branch("fPhoton2MCPz",&fPhoton2MCPz,"fPhoton2MCPz/F");
    }
    if(fSaveAdditionalBranches){
      fSigmaCandTree->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
      fSigmaCandTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaCandTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaCandTree->Branch("fSigY",&fSigY,"fSigY/F");
      fSigmaCandTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
      fSigmaCandTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
      fSigmaCandTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
      fSigmaCandTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
      fSigmaCandTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fSigmaCandTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaCandTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
      fSigmaCandTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
      fSigmaCandTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
      fSigmaCandTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
      fSigmaCandTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
      fSigmaCandTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
      fSigmaCandTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
      fSigmaCandTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
      fSigmaCandTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
      fSigmaCandTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaCandTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
      fSigmaCandTree->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaCandTree->Branch("fPhoton2DCASV",&fPhoton2DCASV,"fPhoton2DCASV/F");
      fSigmaCandTree->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      fSigmaCandTree->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaCandTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
      fSigmaCandTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
      fSigmaCandTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
      fSigmaCandTree->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
      fSigmaCandTree->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
      fSigmaCandTree->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
      fSigmaCandTree->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaCandTree->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaCandTree->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
      fSigmaCandTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
      fSigmaCandTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
      fSigmaCandTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
      fSigmaCandTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
      fSigmaCandTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
      fSigmaCandTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
      fSigmaCandTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
      fSigmaCandTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
      fSigmaCandTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }
    
    // Create TTree of Sigma Candidate Mixed Event Background
    fSigmaMEBackgroundTree = new TTree("fSigmaMEBackgroundTree","Tree of Sigma Mixed Event Background");
    fSigmaMEBackgroundTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaMEBackgroundTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaMEBackgroundTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaMEBackgroundTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaMEBackgroundTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaMEBackgroundTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaMEBackgroundTree->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaMEBackgroundTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaMEBackgroundTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaMEBackgroundTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaMEBackgroundTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaMEBackgroundTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaMEBackgroundTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaMEBackgroundTree->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaMEBackgroundTree->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaMEBackgroundTree->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaMEBackgroundTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaMEBackgroundTree->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaMEBackgroundTree->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaMEBackgroundTree->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaMEBackgroundTree->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaMEBackgroundTree->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    //fSigmaMEBackgroundTree->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaMEBackgroundTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaMEBackgroundTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaMEBackgroundTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaMEBackgroundTree->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaMEBackgroundTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaMEBackgroundTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaMEBackgroundTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaMEBackgroundTree->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaMEBackgroundTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaMEBackgroundTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaMEBackgroundTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaMEBackgroundTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaMEBackgroundTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaMEBackgroundTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    if(fSaveMCBranches){
        fSigmaMEBackgroundTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
        fSigmaMEBackgroundTree->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
        fSigmaMEBackgroundTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
        fSigmaMEBackgroundTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
        fSigmaMEBackgroundTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
        fSigmaMEBackgroundTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
        fSigmaMEBackgroundTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
        fSigmaMEBackgroundTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
        fSigmaMEBackgroundTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
        fSigmaMEBackgroundTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    }
    if(fSaveAdditionalBranches){
        fSigmaMEBackgroundTree->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
        fSigmaMEBackgroundTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
        fSigmaMEBackgroundTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
        fSigmaMEBackgroundTree->Branch("fSigY",&fSigY,"fSigY/F");
        fSigmaMEBackgroundTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
        fSigmaMEBackgroundTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
        fSigmaMEBackgroundTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
        fSigmaMEBackgroundTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
        fSigmaMEBackgroundTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
        fSigmaMEBackgroundTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
        fSigmaMEBackgroundTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
        fSigmaMEBackgroundTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
        fSigmaMEBackgroundTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
        fSigmaMEBackgroundTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
        fSigmaMEBackgroundTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
        fSigmaMEBackgroundTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
        fSigmaMEBackgroundTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
        fSigmaMEBackgroundTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
        fSigmaMEBackgroundTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
        fSigmaMEBackgroundTree->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
        fSigmaMEBackgroundTree->Branch("fPhoton2DCASV",&fPhoton2DCASV,"fPhoton2DCASV/F");
        fSigmaMEBackgroundTree->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
        fSigmaMEBackgroundTree->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
        fSigmaMEBackgroundTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
        fSigmaMEBackgroundTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
        fSigmaMEBackgroundTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
        fSigmaMEBackgroundTree->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
        fSigmaMEBackgroundTree->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
        fSigmaMEBackgroundTree->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
        fSigmaMEBackgroundTree->Branch("fProtonX",&fProtonX,"fProtonX/F");
        fSigmaMEBackgroundTree->Branch("fProtonY",&fProtonY,"fProtonY/F");
        fSigmaMEBackgroundTree->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
        fSigmaMEBackgroundTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
        fSigmaMEBackgroundTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
        fSigmaMEBackgroundTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
        fSigmaMEBackgroundTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }

    // Create TTree of Sigma Proton Pairs in Same Event
    fSigmaPairTreeSE = new TTree("fSigmaPairTreeSE","Tree of Sigma Proton Pairs in Same Event");
    fSigmaPairTreeSE->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreeSE->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaPairTreeSE->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreeSE->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreeSE->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaPairTreeSE->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreeSE->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreeSE->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPairTreeSE->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreeSE->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreeSE->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPairTreeSE->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreeSE->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreeSE->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPairTreeSE->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPairTreeSE->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPairTreeSE->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPairTreeSE->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPairTreeSE->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaPairTreeSE->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaPairTreeSE->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaPairTreeSE->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPairTreeSE->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPairTreeSE->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaPairTreeSE->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaPairTreeSE->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPairTreeSE->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPairTreeSE->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPairTreeSE->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPairTreeSE->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPairTreeSE->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaPairTreeSE->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaPairTreeSE->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPairTreeSE->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaPairTreeSE->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreeSE->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreeSE->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPairTreeSE->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPairTreeSE->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreeSE->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreeSE->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreeSE->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPairTreeSE->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPairTreeSE->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaPairTreeSE->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreeSE->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreeSE->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreeSE->Branch("fPairProtonFilterMap",&fPairProtonFilterMap,"fPairProtonFilterMap/i");
    fSigmaPairTreeSE->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/S");
    fSigmaPairTreeSE->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreeSE->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreeSE->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fSigmaPairTreeSE->Branch("fPairProtonP",&fPairProtonP,"fPairProtonP/F");
    fSigmaPairTreeSE->Branch("fPairProtonEta",&fPairProtonEta,"fPairProtonEta/F");
    fSigmaPairTreeSE->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreeSE->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreeSE->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreeSE->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreeSE->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/S");
    fSigmaPairTreeSE->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/S");
    fSigmaPairTreeSE->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    if(fSaveMCBranches){
        fSigmaPairTreeSE->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
        fSigmaPairTreeSE->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
        fSigmaPairTreeSE->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
        fSigmaPairTreeSE->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
        fSigmaPairTreeSE->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
        fSigmaPairTreeSE->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
        fSigmaPairTreeSE->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
        fSigmaPairTreeSE->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
        fSigmaPairTreeSE->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
        fSigmaPairTreeSE->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
        fSigmaPairTreeSE->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    }
    if(fSaveAdditionalBranches){
        fSigmaPairTreeSE->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
        fSigmaPairTreeSE->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
        fSigmaPairTreeSE->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
        fSigmaPairTreeSE->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
        fSigmaPairTreeSE->Branch("fSigY",&fSigY,"fSigY/F");
        fSigmaPairTreeSE->Branch("fSigPA",&fSigPA,"fSigPA/F");
        fSigmaPairTreeSE->Branch("fSigPx",&fSigPx,"fSigPx/F");
        fSigmaPairTreeSE->Branch("fSigPy",&fSigPy,"fSigPy/F");
        fSigmaPairTreeSE->Branch("fSigPz",&fSigPz,"fSigPz/F");
        fSigmaPairTreeSE->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
        fSigmaPairTreeSE->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
        fSigmaPairTreeSE->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
        fSigmaPairTreeSE->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
        fSigmaPairTreeSE->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
        fSigmaPairTreeSE->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
        fSigmaPairTreeSE->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
        fSigmaPairTreeSE->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
        fSigmaPairTreeSE->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
        fSigmaPairTreeSE->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
        fSigmaPairTreeSE->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
        fSigmaPairTreeSE->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
        fSigmaPairTreeSE->Branch("fPhoton2DCASV",&fPhoton2DCASV,"fPhoton2DCASV/F");
        fSigmaPairTreeSE->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
        fSigmaPairTreeSE->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
        //fSigmaPairTreeSE->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
        fSigmaPairTreeSE->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
        fSigmaPairTreeSE->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
        fSigmaPairTreeSE->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
        fSigmaPairTreeSE->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
        fSigmaPairTreeSE->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
        fSigmaPairTreeSE->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
        fSigmaPairTreeSE->Branch("fProtonX",&fProtonX,"fProtonX/F");
        fSigmaPairTreeSE->Branch("fProtonY",&fProtonY,"fProtonY/F");
        fSigmaPairTreeSE->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
        fSigmaPairTreeSE->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
        fSigmaPairTreeSE->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
        fSigmaPairTreeSE->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
        fSigmaPairTreeSE->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    }

    // Create TTree of Sigma Proton Pairs in Mixed Event
    fSigmaPairTreeME = new TTree("fSigmaPairTreeME","Tree of Sigma Proton Pairs in Mixed Event");
    fSigmaPairTreeME->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreeME->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaPairTreeME->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreeME->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreeME->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaPairTreeME->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreeME->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreeME->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPairTreeME->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreeME->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreeME->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPairTreeME->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreeME->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreeME->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPairTreeME->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPairTreeME->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPairTreeME->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPairTreeME->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPairTreeME->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaPairTreeME->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaPairTreeME->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaPairTreeME->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPairTreeME->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPairTreeME->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaPairTreeME->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaPairTreeME->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPairTreeME->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPairTreeME->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPairTreeME->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPairTreeME->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPairTreeME->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaPairTreeME->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaPairTreeME->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPairTreeME->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaPairTreeME->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreeME->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreeME->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPairTreeME->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPairTreeME->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreeME->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreeME->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreeME->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPairTreeME->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPairTreeME->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaPairTreeME->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreeME->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreeME->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreeME->Branch("fPairProtonFilterMap",&fPairProtonFilterMap,"fPairProtonFilterMap/i");
    fSigmaPairTreeME->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/S");
    fSigmaPairTreeME->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreeME->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreeME->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fSigmaPairTreeME->Branch("fPairProtonP",&fPairProtonP,"fPairProtonP/F");
    fSigmaPairTreeME->Branch("fPairProtonEta",&fPairProtonEta,"fPairProtonEta/F");
    fSigmaPairTreeME->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreeME->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreeME->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreeME->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreeME->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/S");
    fSigmaPairTreeME->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/S");
    fSigmaPairTreeME->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    if(fSaveMCBranches){
        fSigmaPairTreeME->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
        fSigmaPairTreeME->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
        fSigmaPairTreeME->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
        fSigmaPairTreeME->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
        fSigmaPairTreeME->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
        fSigmaPairTreeME->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
        fSigmaPairTreeME->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
        fSigmaPairTreeME->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
        fSigmaPairTreeME->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
        fSigmaPairTreeME->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
        fSigmaPairTreeME->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
        fSigmaPairTreeME->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
        fSigmaPairTreeME->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
        fSigmaPairTreeME->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    }
    if(fSaveAdditionalBranches){
        fSigmaPairTreeME->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
        fSigmaPairTreeME->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
        fSigmaPairTreeME->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
        fSigmaPairTreeME->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
        fSigmaPairTreeME->Branch("fSigY",&fSigY,"fSigY/F");
        fSigmaPairTreeME->Branch("fSigPA",&fSigPA,"fSigPA/F");
        fSigmaPairTreeME->Branch("fSigPx",&fSigPx,"fSigPx/F");
        fSigmaPairTreeME->Branch("fSigPy",&fSigPy,"fSigPy/F");
        fSigmaPairTreeME->Branch("fSigPz",&fSigPz,"fSigPz/F");
        fSigmaPairTreeME->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
        fSigmaPairTreeME->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
        fSigmaPairTreeME->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
        fSigmaPairTreeME->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
        fSigmaPairTreeME->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
        fSigmaPairTreeME->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
        fSigmaPairTreeME->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
        fSigmaPairTreeME->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
        fSigmaPairTreeME->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
        fSigmaPairTreeME->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
        fSigmaPairTreeME->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
        fSigmaPairTreeME->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
        fSigmaPairTreeME->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
        fSigmaPairTreeME->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
        fSigmaPairTreeME->Branch("fPhoton2DCASV",&fPhoton2DCASV,"fPhoton2DCASV/F");
        fSigmaPairTreeME->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
        fSigmaPairTreeME->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
        //fSigmaPairTreeME->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
        fSigmaPairTreeME->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
        fSigmaPairTreeME->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
        fSigmaPairTreeME->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
        fSigmaPairTreeME->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
        fSigmaPairTreeME->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
        fSigmaPairTreeME->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
        fSigmaPairTreeME->Branch("fProtonX",&fProtonX,"fProtonX/F");
        fSigmaPairTreeME->Branch("fProtonY",&fProtonY,"fProtonY/F");
        fSigmaPairTreeME->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
        fSigmaPairTreeME->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
        fSigmaPairTreeME->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
        fSigmaPairTreeME->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
        fSigmaPairTreeME->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    }

    // Create TTree of PHOS Sigma Candidates
    fSigmaPHOSCandTree = new TTree("fSigmaPHOSCandTree","Tree of PHOS Sigma Candidates");
    fSigmaPHOSCandTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPHOSCandTree->Branch("fIsClusterEMCAL",&fIsClusterEMCAL,"fIsClusterEMCAL/O");
    fSigmaPHOSCandTree->Branch("fIsClusterPHOS",&fIsClusterPHOS,"fIsClusterPHOS/O");
    fSigmaPHOSCandTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPHOSCandTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPHOSCandTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPHOSCandTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPHOSCandTree->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPHOSCandTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPHOSCandTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPHOSCandTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPHOSCandTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPHOSCandTree->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaPHOSCandTree->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaPHOSCandTree->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaPHOSCandTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPHOSCandTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPHOSCandTree->Branch("fInvSigMassUncorr",&fInvSigMassUncorr,"fInvSigMassUncorr/F");
    fSigmaPHOSCandTree->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPHOSCandTree->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPHOSCandTree->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPHOSCandTree->Branch("fConvPhotonSecPA",&fConvPhotonSecPA,"fConvPhotonSecPA/F");
    fSigmaPHOSCandTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPHOSCandTree->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPHOSCandTree->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPHOSCandTree->Branch("fCaloPhotonE",&fCaloPhotonE,"fCaloPhotonE/F");
    fSigmaPHOSCandTree->Branch("fCaloPhotonELead",&fCaloPhotonELead,"fCaloPhotonELead/F");
    fSigmaPHOSCandTree->Branch("fClustTrackDx",&fClustTrackDx,"fClustTrackDx/F");
    fSigmaPHOSCandTree->Branch("fClustTrackDz",&fClustTrackDz,"fClustTrackDz/F");
    fSigmaPHOSCandTree->Branch("fClustTrackD",&fClustTrackD,"fClustTrackD/F");
    fSigmaPHOSCandTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");    
    fSigmaPHOSCandTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPHOSCandTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPHOSCandTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPHOSCandTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPHOSCandTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPHOSCandTree->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPHOSCandTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPHOSCandTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPHOSCandTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPHOSCandTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPHOSCandTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPHOSCandTree->Branch("fClustM20",&fClustM20,"fClustM20/F");
    fSigmaPHOSCandTree->Branch("fClustM02",&fClustM02,"fClustM02/F");
    fSigmaPHOSCandTree->Branch("fClustNCells",&fClustNCells,"fClustNCells/I");
    fSigmaPHOSCandTree->Branch("fClustDisttoBC",&fClustDisttoBC,"fClustDisttoBC/F");
    fSigmaPHOSCandTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPHOSCandTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPHOSCandTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPHOSCandTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPHOSCandTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPHOSCandTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPHOSCandTree->Branch("fCaloPhotonX",&fCaloPhotonX,"fCaloPhotonX/F");
    fSigmaPHOSCandTree->Branch("fCaloPhotonY",&fCaloPhotonY,"fCaloPhotonY/F");
    fSigmaPHOSCandTree->Branch("fCaloPhotonZ",&fCaloPhotonZ,"fCaloPhotonZ/F");
    fSigmaPHOSCandTree->Branch("fCellsAbsId",&fCellsAbsId,"fCellsAbsId[20]/S");
    fSigmaPHOSCandTree->Branch("fCaloPhotonEcorr",&fCaloPhotonEcorr,"fCaloPhotonEcorr/F");
    fSigmaPHOSCandTree->Branch("fClustDispersion",&fClustDispersion,"fClustDispersion/F");
    fSigmaPHOSCandTree->Branch("fClustNTracksMatched",&fClustNTracksMatched,"fClustNTracksMatched/I");
    fSigmaPHOSCandTree->Branch("fClustTOF",&fClustTOF,"fClustTOF/F");
    fSigmaPHOSCandTree->Branch("fClustBeta",&fClustBeta,"fClustBeta/F");
    if(fSaveMCBranches){
      fSigmaPHOSCandTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaPHOSCandTree->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaPHOSCandTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaPHOSCandTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaPHOSCandTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaPHOSCandTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaPHOSCandTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaPHOSCandTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaPHOSCandTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaPHOSCandTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
      fSigmaPHOSCandTree->Branch("fCaloPhotonPxMC",&fCaloPhotonPxMC,"fCaloPhotonPxMC/F");
      fSigmaPHOSCandTree->Branch("fCaloPhotonPyMC",&fCaloPhotonPyMC,"fCaloPhotonPyMC/F");
      fSigmaPHOSCandTree->Branch("fCaloPhotonPzMC",&fCaloPhotonPzMC,"fCaloPhotonPzMC/F");
      fSigmaPHOSCandTree->Branch("fCaloPhotonEMC",&fCaloPhotonEMC,"fCaloPhotonEMC/F");
      fSigmaPHOSCandTree->Branch("fClustNLabels",&fClustNLabels,"fClustNLabels/I");
      fSigmaPHOSCandTree->Branch("fClustPDGCode",&fClustPDGCode,"fClustPDGCode/I");
    }
    if(fSaveAddMCBranches){
      fSigmaPHOSCandTree->Branch("fTrackLabel",&fTrackLabel,"fTrackLabel/I");
      fSigmaPHOSCandTree->Branch("fTrackPDGCode",&fTrackPDGCode,"fTrackPDGCode/I");
      fSigmaPHOSCandTree->Branch("fTrackMotherID",&fTrackMotherID,"fTrackMotherID/I");
      fSigmaPHOSCandTree->Branch("fTrackMotherPDGCode",&fTrackMotherPDGCode,"fTrackMotherPDGCode/I");
      fSigmaPHOSCandTree->Branch("fTrackMotherMCPx",&fTrackMotherMCPx,"fTrackMotherMCPx/F");
      fSigmaPHOSCandTree->Branch("fTrackMotherMCPy",&fTrackMotherMCPy,"fTrackMotherMCPy/F");
      fSigmaPHOSCandTree->Branch("fTrackMotherMCPz",&fTrackMotherMCPz,"fTrackMotherMCPz/F");
      fSigmaPHOSCandTree->Branch("fTrackMCPx",&fTrackMCPx,"fTrackMCPx/F");
      fSigmaPHOSCandTree->Branch("fTrackMCPy",&fTrackMCPy,"fTrackMCPy/F");
      fSigmaPHOSCandTree->Branch("fTrackMCPz",&fTrackMCPz,"fTrackMCPz/F");
      fSigmaPHOSCandTree->Branch("fPhoton1Label",&fPhoton1Label,"fPhoton1Label/I");
      fSigmaPHOSCandTree->Branch("fPhoton1PDGCode",&fPhoton1PDGCode,"fPhoton1PDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton1MotherID",&fPhoton1MotherID,"fPhoton1MotherID/I");
      fSigmaPHOSCandTree->Branch("fPhoton1MotherPDGCode",&fPhoton1MotherPDGCode,"fPhoton1MotherPDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton1GMotherID",&fPhoton1GMotherID,"fPhoton1GMotherID/I");
      fSigmaPHOSCandTree->Branch("fPhoton1GMotherPDGCode",&fPhoton1GMotherPDGCode,"fPhoton1GMotherPDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton1MCPx",&fPhoton1MCPx,"fPhoton1MCPx/F");
      fSigmaPHOSCandTree->Branch("fPhoton1MCPy",&fPhoton1MCPy,"fPhoton1MCPy/F");
      fSigmaPHOSCandTree->Branch("fPhoton1MCPz",&fPhoton1MCPz,"fPhoton1MCPz/F");
      fSigmaPHOSCandTree->Branch("fPhoton2Label",&fPhoton2Label,"fPhoton2Label/I");
      fSigmaPHOSCandTree->Branch("fPhoton2PDGCode",&fPhoton2PDGCode,"fPhoton2PDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton2MotherID",&fPhoton2MotherID,"fPhoton2MotherID/I");
      fSigmaPHOSCandTree->Branch("fPhoton2MotherPDGCode",&fPhoton2MotherPDGCode,"fPhoton2MotherPDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton2GMotherID",&fPhoton2GMotherID,"fPhoton2GMotherID/I");
      fSigmaPHOSCandTree->Branch("fPhoton2GMotherPDGCode",&fPhoton2GMotherPDGCode,"fPhoton2GMotherPDGCode/I");
      fSigmaPHOSCandTree->Branch("fPhoton2MCPx",&fPhoton2MCPx,"fPhoton2MCPx/F");
      fSigmaPHOSCandTree->Branch("fPhoton2MCPy",&fPhoton2MCPy,"fPhoton2MCPy/F");
      fSigmaPHOSCandTree->Branch("fPhoton2MCPz",&fPhoton2MCPz,"fPhoton2MCPz/F");
    }
    if(fSaveAdditionalBranches){
      fSigmaPHOSCandTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaPHOSCandTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaPHOSCandTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
      fSigmaPHOSCandTree->Branch("fSigY",&fSigY,"fSigY/F");
      fSigmaPHOSCandTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
      fSigmaPHOSCandTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
      fSigmaPHOSCandTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
      fSigmaPHOSCandTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
      fSigmaPHOSCandTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fSigmaPHOSCandTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
      fSigmaPHOSCandTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
      fSigmaPHOSCandTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
      fSigmaPHOSCandTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
      fSigmaPHOSCandTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
      fSigmaPHOSCandTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
      fSigmaPHOSCandTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
      fSigmaPHOSCandTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
      fSigmaPHOSCandTree->Branch("fConvPhotonX",&fConvPhotonX,"fConvPhotonX/F");
      fSigmaPHOSCandTree->Branch("fConvPhotonY",&fConvPhotonY,"fConvPhotonY/F");
      fSigmaPHOSCandTree->Branch("fConvPhotonZ",&fConvPhotonZ,"fConvPhotonZ/F");
      fSigmaPHOSCandTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaPHOSCandTree->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaPHOSCandTree->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      //fSigmaPHOSCandTree->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
      fSigmaPHOSCandTree->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaPHOSCandTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
      fSigmaPHOSCandTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
      fSigmaPHOSCandTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
      fSigmaPHOSCandTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
      fSigmaPHOSCandTree->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaPHOSCandTree->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaPHOSCandTree->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
      fSigmaPHOSCandTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
      fSigmaPHOSCandTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
      fSigmaPHOSCandTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
      fSigmaPHOSCandTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }

    // Create TTree of PHOS Sigma Candidates
    fSigmaPHOSMEBkgTree = new TTree("fSigmaPHOSMEBkgTree","Tree of Sigma PHOS Mixed Event Background");
    fSigmaPHOSMEBkgTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPHOSMEBkgTree->Branch("fIsClusterEMCAL",&fIsClusterEMCAL,"fIsClusterEMCAL/O");
    fSigmaPHOSMEBkgTree->Branch("fIsClusterPHOS",&fIsClusterPHOS,"fIsClusterPHOS/O");
    fSigmaPHOSMEBkgTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPHOSMEBkgTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPHOSMEBkgTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPHOSMEBkgTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPHOSMEBkgTree->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPHOSMEBkgTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPHOSMEBkgTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPHOSMEBkgTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPHOSMEBkgTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPHOSMEBkgTree->Branch("fInvSigMassUncorr",&fInvSigMassUncorr,"fInvSigMassUncorr/F");
    fSigmaPHOSMEBkgTree->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPHOSMEBkgTree->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPHOSMEBkgTree->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPHOSMEBkgTree->Branch("fConvPhotonSecPA",&fConvPhotonSecPA,"fConvPhotonSecPA/F");
    fSigmaPHOSMEBkgTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPHOSMEBkgTree->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPHOSMEBkgTree->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonE",&fCaloPhotonE,"fCaloPhotonE/F");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonELead",&fCaloPhotonELead,"fCaloPhotonELead/F");
    fSigmaPHOSMEBkgTree->Branch("fClustTrackDx",&fClustTrackDx,"fClustTrackDx/F");
    fSigmaPHOSMEBkgTree->Branch("fClustTrackDz",&fClustTrackDz,"fClustTrackDz/F");
    fSigmaPHOSMEBkgTree->Branch("fClustTrackD",&fClustTrackD,"fClustTrackD/F");
    fSigmaPHOSMEBkgTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");    
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPHOSMEBkgTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPHOSMEBkgTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPHOSMEBkgTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPHOSMEBkgTree->Branch("fClustM20",&fClustM20,"fClustM20/F");
    fSigmaPHOSMEBkgTree->Branch("fClustM02",&fClustM02,"fClustM02/F");
    fSigmaPHOSMEBkgTree->Branch("fClustNCells",&fClustNCells,"fClustNCells/I");
    fSigmaPHOSMEBkgTree->Branch("fClustDisttoBC",&fClustDisttoBC,"fClustDisttoBC/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPHOSMEBkgTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPHOSMEBkgTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonX",&fCaloPhotonX,"fCaloPhotonX/F");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonY",&fCaloPhotonY,"fCaloPhotonY/F");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonZ",&fCaloPhotonZ,"fCaloPhotonZ/F");
    fSigmaPHOSMEBkgTree->Branch("fCellsAbsId",&fCellsAbsId,"fCellsAbsId[20]/S");
    fSigmaPHOSMEBkgTree->Branch("fCaloPhotonEcorr",&fCaloPhotonEcorr,"fCaloPhotonEcorr/F");
    fSigmaPHOSMEBkgTree->Branch("fClustDispersion",&fClustDispersion,"fClustDispersion/F");
    fSigmaPHOSMEBkgTree->Branch("fClustNTracksMatched",&fClustNTracksMatched,"fClustNTracksMatched/I");
    fSigmaPHOSMEBkgTree->Branch("fClustTOF",&fClustTOF,"fClustTOF/F");
    fSigmaPHOSMEBkgTree->Branch("fClustBeta",&fClustBeta,"fClustBeta/F");
    if(fSaveMCBranches){
      fSigmaPHOSMEBkgTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaPHOSMEBkgTree->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaPHOSMEBkgTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaPHOSMEBkgTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
      fSigmaPHOSMEBkgTree->Branch("fCaloPhotonPxMC",&fCaloPhotonPxMC,"fCaloPhotonPxMC/F");
      fSigmaPHOSMEBkgTree->Branch("fCaloPhotonPyMC",&fCaloPhotonPyMC,"fCaloPhotonPyMC/F");
      fSigmaPHOSMEBkgTree->Branch("fCaloPhotonPzMC",&fCaloPhotonPzMC,"fCaloPhotonPzMC/F");
      fSigmaPHOSMEBkgTree->Branch("fCaloPhotonEMC",&fCaloPhotonEMC,"fCaloPhotonEMC/F");
      fSigmaPHOSMEBkgTree->Branch("fClustNLabels",&fClustNLabels,"fClustNLabels/I");
      fSigmaPHOSMEBkgTree->Branch("fClustPDGCode",&fClustPDGCode,"fClustPDGCode/I");
    }
    if(fSaveAdditionalBranches){
      fSigmaPHOSMEBkgTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaPHOSMEBkgTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
      fSigmaPHOSMEBkgTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaPHOSMEBkgTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
      fSigmaPHOSMEBkgTree->Branch("fSigY",&fSigY,"fSigY/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
      fSigmaPHOSMEBkgTree->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaPHOSMEBkgTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
      fSigmaPHOSMEBkgTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
      fSigmaPHOSMEBkgTree->Branch("fConvPhotonX",&fConvPhotonX,"fConvPhotonX/F");
      fSigmaPHOSMEBkgTree->Branch("fConvPhotonY",&fConvPhotonY,"fConvPhotonY/F");
      fSigmaPHOSMEBkgTree->Branch("fConvPhotonZ",&fConvPhotonZ,"fConvPhotonZ/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaPHOSMEBkgTree->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      //fSigmaPHOSMEBkgTree->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
      fSigmaPHOSMEBkgTree->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaPHOSMEBkgTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
      fSigmaPHOSMEBkgTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
      fSigmaPHOSMEBkgTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
      fSigmaPHOSMEBkgTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
      fSigmaPHOSMEBkgTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }

    // Create TTree of Sigma Proton Pairs in Same Event PHOS
    fSigmaPairTreePHOSSE = new TTree("fSigmaPairTreePHOSSE","Tree of Sigma Proton Pairs in Same Event PHOS");
    fSigmaPairTreePHOSSE->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreePHOSSE->Branch("fIsClusterEMCAL",&fIsClusterEMCAL,"fIsClusterEMCAL/O");
    fSigmaPairTreePHOSSE->Branch("fIsClusterPHOS",&fIsClusterPHOS,"fIsClusterPHOS/O");
    fSigmaPairTreePHOSSE->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreePHOSSE->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreePHOSSE->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreePHOSSE->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreePHOSSE->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPairTreePHOSSE->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreePHOSSE->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreePHOSSE->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreePHOSSE->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreePHOSSE->Branch("fInvSigMassUncorr",&fInvSigMassUncorr,"fInvSigMassUncorr/F");
    fSigmaPairTreePHOSSE->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPairTreePHOSSE->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPairTreePHOSSE->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPairTreePHOSSE->Branch("fConvPhotonSecPA",&fConvPhotonSecPA,"fConvPhotonSecPA/F");
    fSigmaPairTreePHOSSE->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPairTreePHOSSE->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPairTreePHOSSE->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonE",&fCaloPhotonE,"fCaloPhotonE/F");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonELead",&fCaloPhotonELead,"fCaloPhotonELead/F");
    fSigmaPairTreePHOSSE->Branch("fClustTrackDx",&fClustTrackDx,"fClustTrackDx/F");
    fSigmaPairTreePHOSSE->Branch("fClustTrackDz",&fClustTrackDz,"fClustTrackDz/F");
    fSigmaPairTreePHOSSE->Branch("fClustTrackD",&fClustTrackD,"fClustTrackD/F");
    fSigmaPairTreePHOSSE->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");    
    fSigmaPairTreePHOSSE->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPairTreePHOSSE->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPairTreePHOSSE->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPairTreePHOSSE->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPairTreePHOSSE->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPairTreePHOSSE->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreePHOSSE->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreePHOSSE->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreePHOSSE->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreePHOSSE->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPairTreePHOSSE->Branch("fClustM20",&fClustM20,"fClustM20/F");
    fSigmaPairTreePHOSSE->Branch("fClustM02",&fClustM02,"fClustM02/F");
    fSigmaPairTreePHOSSE->Branch("fClustNCells",&fClustNCells,"fClustNCells/I");
    fSigmaPairTreePHOSSE->Branch("fClustDisttoBC",&fClustDisttoBC,"fClustDisttoBC/F");
    fSigmaPairTreePHOSSE->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreePHOSSE->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPairTreePHOSSE->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    fSigmaPairTreePHOSSE->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreePHOSSE->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreePHOSSE->Branch("fPairProtonFilterMap",&fPairProtonFilterMap,"fPairProtonFilterMap/i");
    fSigmaPairTreePHOSSE->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/S");
    fSigmaPairTreePHOSSE->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonP",&fPairProtonP,"fPairProtonP/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonEta",&fPairProtonEta,"fPairProtonEta/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/S");
    fSigmaPairTreePHOSSE->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/S");
    fSigmaPairTreePHOSSE->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreePHOSSE->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreePHOSSE->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPairTreePHOSSE->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreePHOSSE->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonX",&fCaloPhotonX,"fCaloPhotonX/F");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonY",&fCaloPhotonY,"fCaloPhotonY/F");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonZ",&fCaloPhotonZ,"fCaloPhotonZ/F");
    fSigmaPairTreePHOSSE->Branch("fCellsAbsId",&fCellsAbsId,"fCellsAbsId[20]/S");
    fSigmaPairTreePHOSSE->Branch("fCaloPhotonEcorr",&fCaloPhotonEcorr,"fCaloPhotonEcorr/F");
    fSigmaPairTreePHOSSE->Branch("fClustDispersion",&fClustDispersion,"fClustDispersion/F");
    fSigmaPairTreePHOSSE->Branch("fClustNTracksMatched",&fClustNTracksMatched,"fClustNTracksMatched/I");
    fSigmaPairTreePHOSSE->Branch("fClustTOF",&fClustTOF,"fClustTOF/F");
    fSigmaPairTreePHOSSE->Branch("fClustBeta",&fClustBeta,"fClustBeta/F");
    if(fSaveMCBranches){
      fSigmaPairTreePHOSSE->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaPairTreePHOSSE->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaPairTreePHOSSE->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaPairTreePHOSSE->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaPairTreePHOSSE->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaPairTreePHOSSE->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaPairTreePHOSSE->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaPairTreePHOSSE->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
      fSigmaPairTreePHOSSE->Branch("fCaloPhotonPxMC",&fCaloPhotonPxMC,"fCaloPhotonPxMC/F");
      fSigmaPairTreePHOSSE->Branch("fCaloPhotonPyMC",&fCaloPhotonPyMC,"fCaloPhotonPyMC/F");
      fSigmaPairTreePHOSSE->Branch("fCaloPhotonPzMC",&fCaloPhotonPzMC,"fCaloPhotonPzMC/F");
      fSigmaPairTreePHOSSE->Branch("fCaloPhotonEMC",&fCaloPhotonEMC,"fCaloPhotonEMC/F");
      fSigmaPairTreePHOSSE->Branch("fClustNLabels",&fClustNLabels,"fClustNLabels/I");
      fSigmaPairTreePHOSSE->Branch("fClustPDGCode",&fClustPDGCode,"fClustPDGCode/I");
      fSigmaPairTreePHOSSE->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
      fSigmaPairTreePHOSSE->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    }
    if(fSaveAdditionalBranches){
      fSigmaPairTreePHOSSE->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaPairTreePHOSSE->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
      fSigmaPairTreePHOSSE->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaPairTreePHOSSE->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
      fSigmaPairTreePHOSSE->Branch("fSigY",&fSigY,"fSigY/F");
      fSigmaPairTreePHOSSE->Branch("fSigPA",&fSigPA,"fSigPA/F");
      fSigmaPairTreePHOSSE->Branch("fSigPx",&fSigPx,"fSigPx/F");
      fSigmaPairTreePHOSSE->Branch("fSigPy",&fSigPy,"fSigPy/F");
      fSigmaPairTreePHOSSE->Branch("fSigPz",&fSigPz,"fSigPz/F");
      fSigmaPairTreePHOSSE->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
      fSigmaPairTreePHOSSE->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
      fSigmaPairTreePHOSSE->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
      fSigmaPairTreePHOSSE->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fSigmaPairTreePHOSSE->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaPairTreePHOSSE->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
      fSigmaPairTreePHOSSE->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
      fSigmaPairTreePHOSSE->Branch("fConvPhotonX",&fConvPhotonX,"fConvPhotonX/F");
      fSigmaPairTreePHOSSE->Branch("fConvPhotonY",&fConvPhotonY,"fConvPhotonY/F");
      fSigmaPairTreePHOSSE->Branch("fConvPhotonZ",&fConvPhotonZ,"fConvPhotonZ/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaPairTreePHOSSE->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      //fSigmaPairTreePHOSSE->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
      fSigmaPairTreePHOSSE->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaPairTreePHOSSE->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
      fSigmaPairTreePHOSSE->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
      fSigmaPairTreePHOSSE->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
      fSigmaPairTreePHOSSE->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
      fSigmaPairTreePHOSSE->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
      fSigmaPairTreePHOSSE->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
      fSigmaPairTreePHOSSE->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
      fSigmaPairTreePHOSSE->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaPairTreePHOSSE->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaPairTreePHOSSE->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
      fSigmaPairTreePHOSSE->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }

    // Create TTree of Sigma Proton Pairs in Mixed Event PHOS
    fSigmaPairTreePHOSME = new TTree("fSigmaPairTreePHOSME","Tree of Sigma Proton Pairs in Mixed Event PHOS");
    fSigmaPairTreePHOSME->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreePHOSME->Branch("fIsClusterEMCAL",&fIsClusterEMCAL,"fIsClusterEMCAL/O");
    fSigmaPairTreePHOSME->Branch("fIsClusterPHOS",&fIsClusterPHOS,"fIsClusterPHOS/O");
    fSigmaPairTreePHOSME->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreePHOSME->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreePHOSME->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreePHOSME->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreePHOSME->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaPairTreePHOSME->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreePHOSME->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreePHOSME->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreePHOSME->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreePHOSME->Branch("fInvSigMassUncorr",&fInvSigMassUncorr,"fInvSigMassUncorr/F");
    fSigmaPairTreePHOSME->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaPairTreePHOSME->Branch("fSigPAprop",&fSigPAprop,"fSigPAprop/F");
    fSigmaPairTreePHOSME->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaPairTreePHOSME->Branch("fConvPhotonSecPA",&fConvPhotonSecPA,"fConvPhotonSecPA/F");
    fSigmaPairTreePHOSME->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaPairTreePHOSME->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaPairTreePHOSME->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonE",&fCaloPhotonE,"fCaloPhotonE/F");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonELead",&fCaloPhotonELead,"fCaloPhotonELead/F");
    fSigmaPairTreePHOSME->Branch("fClustTrackDx",&fClustTrackDx,"fClustTrackDx/F");
    fSigmaPairTreePHOSME->Branch("fClustTrackDz",&fClustTrackDz,"fClustTrackDz/F");
    fSigmaPairTreePHOSME->Branch("fClustTrackD",&fClustTrackD,"fClustTrackD/F");
    fSigmaPairTreePHOSME->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");    
    fSigmaPairTreePHOSME->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaPairTreePHOSME->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaPairTreePHOSME->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaPairTreePHOSME->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaPairTreePHOSME->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPairTreePHOSME->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaPairTreePHOSME->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreePHOSME->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreePHOSME->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreePHOSME->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreePHOSME->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaPairTreePHOSME->Branch("fClustM20",&fClustM20,"fClustM20/F");
    fSigmaPairTreePHOSME->Branch("fClustM02",&fClustM02,"fClustM02/F");
    fSigmaPairTreePHOSME->Branch("fClustNCells",&fClustNCells,"fClustNCells/I");
    fSigmaPairTreePHOSME->Branch("fClustDisttoBC",&fClustDisttoBC,"fClustDisttoBC/F");
    fSigmaPairTreePHOSME->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreePHOSME->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaPairTreePHOSME->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    fSigmaPairTreePHOSME->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreePHOSME->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreePHOSME->Branch("fPairProtonFilterMap",&fPairProtonFilterMap,"fPairProtonFilterMap/i");
    fSigmaPairTreePHOSME->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/S");
    fSigmaPairTreePHOSME->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonP",&fPairProtonP,"fPairProtonP/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonEta",&fPairProtonEta,"fPairProtonEta/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreePHOSME->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreePHOSME->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreePHOSME->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreePHOSME->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/S");
    fSigmaPairTreePHOSME->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/S");
    fSigmaPairTreePHOSME->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreePHOSME->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreePHOSME->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaPairTreePHOSME->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreePHOSME->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonX",&fCaloPhotonX,"fCaloPhotonX/F");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonY",&fCaloPhotonY,"fCaloPhotonY/F");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonZ",&fCaloPhotonZ,"fCaloPhotonZ/F");
    fSigmaPairTreePHOSME->Branch("fCellsAbsId",&fCellsAbsId,"fCellsAbsId[20]/S");
    fSigmaPairTreePHOSME->Branch("fCaloPhotonEcorr",&fCaloPhotonEcorr,"fCaloPhotonEcorr/F");
    fSigmaPairTreePHOSME->Branch("fClustDispersion",&fClustDispersion,"fClustDispersion/F");
    fSigmaPairTreePHOSME->Branch("fClustNTracksMatched",&fClustNTracksMatched,"fClustNTracksMatched/I");
    fSigmaPairTreePHOSME->Branch("fClustTOF",&fClustTOF,"fClustTOF/F");
    fSigmaPairTreePHOSME->Branch("fClustBeta",&fClustBeta,"fClustBeta/F");
    if(fSaveMCBranches){
      fSigmaPairTreePHOSME->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaPairTreePHOSME->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaPairTreePHOSME->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaPairTreePHOSME->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaPairTreePHOSME->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaPairTreePHOSME->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaPairTreePHOSME->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaPairTreePHOSME->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaPairTreePHOSME->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaPairTreePHOSME->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
      fSigmaPairTreePHOSME->Branch("fCaloPhotonPxMC",&fCaloPhotonPxMC,"fCaloPhotonPxMC/F");
      fSigmaPairTreePHOSME->Branch("fCaloPhotonPyMC",&fCaloPhotonPyMC,"fCaloPhotonPyMC/F");
      fSigmaPairTreePHOSME->Branch("fCaloPhotonPzMC",&fCaloPhotonPzMC,"fCaloPhotonPzMC/F");
      fSigmaPairTreePHOSME->Branch("fCaloPhotonEMC",&fCaloPhotonEMC,"fCaloPhotonEMC/F");
      fSigmaPairTreePHOSME->Branch("fClustNLabels",&fClustNLabels,"fClustNLabels/I");
      fSigmaPairTreePHOSME->Branch("fClustPDGCode",&fClustPDGCode,"fClustPDGCode/I");
      fSigmaPairTreePHOSME->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
      fSigmaPairTreePHOSME->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    }
    if(fSaveAdditionalBranches){
      fSigmaPairTreePHOSME->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaPairTreePHOSME->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
      fSigmaPairTreePHOSME->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaPairTreePHOSME->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
      fSigmaPairTreePHOSME->Branch("fSigY",&fSigY,"fSigY/F");
      fSigmaPairTreePHOSME->Branch("fSigPA",&fSigPA,"fSigPA/F");
      fSigmaPairTreePHOSME->Branch("fSigPx",&fSigPx,"fSigPx/F");
      fSigmaPairTreePHOSME->Branch("fSigPy",&fSigPy,"fSigPy/F");
      fSigmaPairTreePHOSME->Branch("fSigPz",&fSigPz,"fSigPz/F");
      fSigmaPairTreePHOSME->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
      fSigmaPairTreePHOSME->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
      fSigmaPairTreePHOSME->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
      fSigmaPairTreePHOSME->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fSigmaPairTreePHOSME->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaPairTreePHOSME->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
      fSigmaPairTreePHOSME->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
      fSigmaPairTreePHOSME->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
      fSigmaPairTreePHOSME->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
      fSigmaPairTreePHOSME->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
      fSigmaPairTreePHOSME->Branch("fConvPhotonX",&fConvPhotonX,"fConvPhotonX/F");
      fSigmaPairTreePHOSME->Branch("fConvPhotonY",&fConvPhotonY,"fConvPhotonY/F");
      fSigmaPairTreePHOSME->Branch("fConvPhotonZ",&fConvPhotonZ,"fConvPhotonZ/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaPairTreePHOSME->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      //fSigmaPairTreePHOSME->Branch("fTrackDCASVKF",&fTrackDCASVKF,"fTrackDCASVKF/F");
      fSigmaPairTreePHOSME->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaPairTreePHOSME->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
      fSigmaPairTreePHOSME->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
      fSigmaPairTreePHOSME->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
      fSigmaPairTreePHOSME->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
      fSigmaPairTreePHOSME->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
      fSigmaPairTreePHOSME->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
      fSigmaPairTreePHOSME->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
      fSigmaPairTreePHOSME->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaPairTreePHOSME->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaPairTreePHOSME->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
      fSigmaPairTreePHOSME->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    }

    fSigmaCalcCandTree = new TTree("fSigmaCalcCandTree","Tree of calculated Sigma Candidates");
    fSigmaCalcCandTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaCalcCandTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCalcCandTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCalcCandTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaCalcCandTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCalcCandTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaCalcCandTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaCalcCandTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaCalcCandTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCalcCandTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/S");
    fSigmaCalcCandTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaCalcCandTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaCalcCandTree->Branch("fSigProtonFilterMap",&fSigProtonFilterMap,"fSigProtonFilterMap/i");
    fSigmaCalcCandTree->Branch("fProtonEta",&fProtonEta,"fProtonEta/F");
    fSigmaCalcCandTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCalcCandTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaCalcCandTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/S");
    fSigmaCalcCandTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/S");
    fSigmaCalcCandTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCalcCandTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCalcCandTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaCalcCandTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaCalcCandTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaCalcCandTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaCalcCandTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaCalcCandTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaCalcCandTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaCalcCandTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaCalcCandTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaCalcCandTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaCalcCandTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaCalcCandTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaCalcCandTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaCalcCandTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaCalcCandTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaCalcCandTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/S");
    fSigmaCalcCandTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/S");
    fSigmaCalcCandTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaCalcCandTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaCalcCandTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaCalcCandTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaCalcCandTree->Branch("fSigPxprop",&fSigPxprop,"fSigPxprop/F");
    fSigmaCalcCandTree->Branch("fSigPyprop",&fSigPyprop,"fSigPyprop/F");
    fSigmaCalcCandTree->Branch("fSigPzprop",&fSigPzprop,"fSigPzprop/F");
    fSigmaCalcCandTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaCalcCandTree->Branch("fSigYprop",&fSigYprop,"fSigYprop/F");
    fSigmaCalcCandTree->Branch("fSigAntiPA",&fSigAntiPA,"fSigAntiPA/F");
    fSigmaCalcCandTree->Branch("fSigPt",&fSigPt,"fSigPt/F");
    fSigmaCalcCandTree->Branch("fSigFlightDist",&fSigFlightDist,"fSigFlightDist/F");
    fSigmaCalcCandTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaCalcCandTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaCalcCandTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");      
    fSigmaCalcCandTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaCalcCandTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaCalcCandTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaCalcCandTree->Branch("fNIter",&fNIter,"fNIter/S");
    fSigmaCalcCandTree->Branch("fDeltaPhi",&fDeltaPhi,"fDeltaPhi/F");
    fSigmaCalcCandTree->Branch("fDeltaTheta",&fDeltaTheta,"fDeltaTheta/F");
    if(fSaveMCBranches){
      fSigmaCalcCandTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
      fSigmaCalcCandTree->Branch("fIsMCDalitz",&fIsMCDalitz,"fIsMCDalitz/O");
      fSigmaCalcCandTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
      fSigmaCalcCandTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
      fSigmaCalcCandTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fSigmaCalcCandTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fSigmaCalcCandTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fSigmaCalcCandTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
      fSigmaCalcCandTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
      fSigmaCalcCandTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
      fSigmaCalcCandTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
      fSigmaCalcCandTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
      fSigmaCalcCandTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
      fSigmaCalcCandTree->Branch("fCaloPhotonPxMC",&fCaloPhotonPxMC,"fCaloPhotonPxMC/F");
      fSigmaCalcCandTree->Branch("fCaloPhotonPyMC",&fCaloPhotonPyMC,"fCaloPhotonPyMC/F");
      fSigmaCalcCandTree->Branch("fCaloPhotonPzMC",&fCaloPhotonPzMC,"fCaloPhotonPzMC/F");
    }
    if(fSaveAddMCBranches){
      fSigmaCalcCandTree->Branch("fTrackLabel",&fTrackLabel,"fTrackLabel/I");
      fSigmaCalcCandTree->Branch("fTrackPDGCode",&fTrackPDGCode,"fTrackPDGCode/I");
      fSigmaCalcCandTree->Branch("fTrackMotherID",&fTrackMotherID,"fTrackMotherID/I");
      fSigmaCalcCandTree->Branch("fTrackMotherPDGCode",&fTrackMotherPDGCode,"fTrackMotherPDGCode/I");
      fSigmaCalcCandTree->Branch("fTrackMotherMCPx",&fTrackMotherMCPx,"fTrackMotherMCPx/F");
      fSigmaCalcCandTree->Branch("fTrackMotherMCPy",&fTrackMotherMCPy,"fTrackMotherMCPy/F");
      fSigmaCalcCandTree->Branch("fTrackMotherMCPz",&fTrackMotherMCPz,"fTrackMotherMCPz/F");
      fSigmaCalcCandTree->Branch("fTrackMCPx",&fTrackMCPx,"fTrackMCPx/F");
      fSigmaCalcCandTree->Branch("fTrackMCPy",&fTrackMCPy,"fTrackMCPy/F");
      fSigmaCalcCandTree->Branch("fTrackMCPz",&fTrackMCPz,"fTrackMCPz/F");
      fSigmaCalcCandTree->Branch("fPhoton1Label",&fPhoton1Label,"fPhoton1Label/I");
      fSigmaCalcCandTree->Branch("fPhoton1PDGCode",&fPhoton1PDGCode,"fPhoton1PDGCode/I");
      fSigmaCalcCandTree->Branch("fPhoton1MotherID",&fPhoton1MotherID,"fPhoton1MotherID/I");
      fSigmaCalcCandTree->Branch("fPhoton1MotherPDGCode",&fPhoton1MotherPDGCode,"fPhoton1MotherPDGCode/I");
      fSigmaCalcCandTree->Branch("fPhoton1GMotherID",&fPhoton1GMotherID,"fPhoton1GMotherID/I");
      fSigmaCalcCandTree->Branch("fPhoton1GMotherPDGCode",&fPhoton1GMotherPDGCode,"fPhoton1GMotherPDGCode/I");
      fSigmaCalcCandTree->Branch("fPhoton1MCPx",&fPhoton1MCPx,"fPhoton1MCPx/F");
      fSigmaCalcCandTree->Branch("fPhoton1MCPy",&fPhoton1MCPy,"fPhoton1MCPy/F");
      fSigmaCalcCandTree->Branch("fPhoton1MCPz",&fPhoton1MCPz,"fPhoton1MCPz/F");
    }
    if(fSaveAdditionalBranches){
      fSigmaCalcCandTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fSigmaCalcCandTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fSigmaCalcCandTree->Branch("fConvPhotonX",&fConvPhotonX,"fConvPhotonX/F");
      fSigmaCalcCandTree->Branch("fConvPhotonY",&fConvPhotonY,"fConvPhotonY/F");
      fSigmaCalcCandTree->Branch("fConvPhotonZ",&fConvPhotonZ,"fConvPhotonZ/F");
      fSigmaCalcCandTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
      fSigmaCalcCandTree->Branch("fPhoton1DCASV",&fPhoton1DCASV,"fPhoton1DCASV/F");
      fSigmaCalcCandTree->Branch("fTrackDCASV",&fTrackDCASV,"fTrackDCASV/F");
      fSigmaCalcCandTree->Branch("fKFChi2",&fKFChi2,"fKFChi2/F");
      fSigmaCalcCandTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
      fSigmaCalcCandTree->Branch("fProtonX",&fProtonX,"fProtonX/F");
      fSigmaCalcCandTree->Branch("fProtonY",&fProtonY,"fProtonY/F");
      fSigmaCalcCandTree->Branch("fProtonZ",&fProtonZ,"fProtonZ/F");
    }

    // Create TTree of Protons
    fProtonTree = new TTree("fProtonTree","Tree of Protons");
    fProtonTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fProtonTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fProtonTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fProtonTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fProtonTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fProtonTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fProtonTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fProtonTree->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fProtonTree->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fProtonTree->Branch("fPairProtonFilterMap",&fPairProtonFilterMap,"fPairProtonFilterMap/i");
    fProtonTree->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/S");
    fProtonTree->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fProtonTree->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fProtonTree->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fProtonTree->Branch("fPairProtonP",&fPairProtonP,"fPairProtonP/F");
    fProtonTree->Branch("fPairProtonEta",&fPairProtonEta,"fPairProtonEta/F");
    fProtonTree->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fProtonTree->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fProtonTree->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fProtonTree->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fProtonTree->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fProtonTree->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fProtonTree->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fProtonTree->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fProtonTree->Branch("fPairProtNSigTPCElec",&fPairProtNSigTPCElec,"fPairProtNSigTPCElec/F");
    fProtonTree->Branch("fPairProtNSigTOFElec",&fPairProtNSigTOFElec,"fPairProtNSigTOFElec/F");
    fProtonTree->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/S");
    fProtonTree->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/S");
    fProtonTree->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    if(fSaveMCBranches){
      fProtonTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
      fProtonTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
      fProtonTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
      fProtonTree->Branch("fPairProtonPxMC",&fPairProtonPxMC,"fPairProtonPxMC/F");
      fProtonTree->Branch("fPairProtonPyMC",&fPairProtonPyMC,"fPairProtonPyMC/F");
      fProtonTree->Branch("fPairProtonPzMC",&fPairProtonPzMC,"fPairProtonPzMC/F");
      fProtonTree->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
      fProtonTree->Branch("fPairProtonPDG",&fPairProtonPDG,"fPairProtonPDG/I");
      fProtonTree->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    }
    if(fSaveAdditionalBranches){
      fProtonTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
      fProtonTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
      fProtonTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
      fProtonTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
      fSigmaCalcCandTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    }

    //Save Particle Masses and other constants for later use
    cElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();      
    cProtonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    
    cSigmaMass    = TDatabasePDG::Instance()->GetParticle(3222)->Mass();    
    cPi0Mass      = TDatabasePDG::Instance()->GetParticle(111)->Mass();     
    c = 2.99792457999999984e-02; // [cm/ps]

/**************************Histograms********************************/

    //Book Keeper for used Cuts 
    TH1D* fHistCutBookKeeper           = new TH1D("fHistCutBookKeeper", "Book Keeper for used Cuts", 117, 0.5, 117.5);

    //Event related                    
    TH1F* fHistMCGenPileup             = new TH1F("fHistMCGenPileup", "Generated pile-up;IsPileUp;", 2, -0.5, 1.5);
    TH1F* fHistVertexZ                 = new TH1F("fHistVertexZ", "Z Vertex Position;z (cm);Counts/mm", 400, -20, 20);
    TH1F* fHistVertexZMC               = new TH1F("fHistVertexZMC", "MC Z Vertex Position;z (cm);Counts/mm", 400, -20, 20);
    TH1F* fHistCentrality              = new TH1F("fHistCentrality", "Centrality Percentile;Centrality (%);Counts/(0.01 %)", 10000, 0, 100);
    TH1F* fHistCentralityHMV0          = new TH1F("fHistCentralityHMV0", "Centrality Percentile - HMV0;Centrality (%);Counts/(0.01 %)", 10000, 0, 100);
    TH1F* fHistCentralityHMSPD         = new TH1F("fHistCentralityHMSPD", "Centrality Percentile - HMSPD;Centrality (%);Counts/(0.01 %)", 10000, 0, 100);
    TH1F* fHistCentralityINT7          = new TH1F("fHistCentralityINT7", "Centrality Percentile - INT7;Centrality (%);Counts/(0.01 %)", 10000, 0, 100);
    TH1F* fHistCentralityINT7LF        = new TH1F("fHistCentralityINT7LF", "Centrality Percentile - INT7 Low Field;Centrality (%);Counts/(0.01 %)", 10000, 0, 100);
    TH1F* fHistEventCounter            = new TH1F("fHistEventCounter", "Event Counter", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdouble      = new TH1D("fHistEventCounterdouble", "Event Counter", 2, 0.5, 2.5);
    TH1F* fHistEventCounterHM          = new TH1F("fHistEventCounterHM", "Event Counter - High-Mult. V0+SPD", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdoubleHM    = new TH1D("fHistEventCounterdoubleHM", "Event Counter - High-Mult. V0+SPD", 2, 0.5, 2.5);
    TH1F* fHistEventCounterHMV0        = new TH1F("fHistEventCounterHMV0", "Event Counter - HMV0", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdoubleHMV0  = new TH1D("fHistEventCounterdoubleHMV0", "Event Counter - HMV0", 2, 0.5, 2.5);
    TH1F* fHistEventCounterHMSPD       = new TH1F("fHistEventCounterHMSPD", "Event Counter - HMSPD", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdoubleHMSPD = new TH1D("fHistEventCounterdoubleHMSPD", "Event Counter - HMSPD", 2, 0.5, 2.5);
    TH1F* fHistEventCounterINT7        = new TH1F("fHistEventCounterINT7", "Event Counter - INT7", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdoubleINT7  = new TH1D("fHistEventCounterdoubleINT7", "Event Counter - INT7", 2, 0.5, 2.5);
    TH1F* fHistEventCounterINT7LF      = new TH1F("fHistEventCounterINT7LF", "Event Counter - INT7 Low Field", 2, 0.5, 2.5);
    TH1D* fHistEventCounterdoubleINT7LF= new TH1D("fHistEventCounterdoubleINT7LF", "Event Counter - INT7 Low Field", 2, 0.5, 2.5);
    TH1F* fHistTrackMultiplicity       = new TH1F("fHistTrackMultiplicity","Number of Tracks per Event",20000,0,20000);
    TH1F* fHistRefMultComb05           = new TH1F("fHistRefMultComb05","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5",20000,0,20000);
    TH1F* fHistRefMultComb08           = new TH1F("fHistRefMultComb08","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8",20000,0,20000);
    TH1F* fHistRefMultComb10           = new TH1F("fHistRefMultComb10","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0",20000,0,20000);
    TH1F* fHistRefMultComb05HM         = new TH1F("fHistRefMultComb05HM","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5 - HM",20000,0,20000);
    TH1F* fHistRefMultComb08HM         = new TH1F("fHistRefMultComb08HM","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8 - HM",20000,0,20000);
    TH1F* fHistRefMultComb10HM         = new TH1F("fHistRefMultComb10HM","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0 - HM",20000,0,20000);
    TH1F* fHistRefMultComb05HMV0       = new TH1F("fHistRefMultComb05HMV0","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5 - HMV0",20000,0,20000);
    TH1F* fHistRefMultComb08HMV0       = new TH1F("fHistRefMultComb08HMV0","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8 - HMV0",20000,0,20000);
    TH1F* fHistRefMultComb10HMV0       = new TH1F("fHistRefMultComb10HMV0","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0 - HMV0",20000,0,20000);
    TH1F* fHistRefMultComb05HMSPD      = new TH1F("fHistRefMultComb05HMSPD","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5 - HMSPD",20000,0,20000);
    TH1F* fHistRefMultComb08HMSPD      = new TH1F("fHistRefMultComb08HMSPD","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8 - HMSPD",20000,0,20000);
    TH1F* fHistRefMultComb10HMSPD      = new TH1F("fHistRefMultComb10HMSPD","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0 - HMSPD",20000,0,20000);
    TH1F* fHistRefMultComb05INT7       = new TH1F("fHistRefMultComb05INT7","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5 - INT7",20000,0,20000);
    TH1F* fHistRefMultComb08INT7       = new TH1F("fHistRefMultComb08INT7","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8 - INT7",20000,0,20000);
    TH1F* fHistRefMultComb10INT7       = new TH1F("fHistRefMultComb10INT7","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0 - INT7",20000,0,20000);
    TH1F* fHistRefMultComb05INT7LF     = new TH1F("fHistRefMultComb05INT7LF","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.5 - INT7 Low Field",20000,0,20000);
    TH1F* fHistRefMultComb08INT7LF     = new TH1F("fHistRefMultComb08INT7LF","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8 - INT7 Low Field",20000,0,20000);
    TH1F* fHistRefMultComb10INT7LF     = new TH1F("fHistRefMultComb10INT7LF","Combined reference multiplicity (tracklets + ITSTPC) in |eta|<1.0 - INT7 Low Field",20000,0,20000);

    //Counters                         
    TH1F* fHistMCCounter               = new TH1F("fHistMCCounter", "Monte Carlo Counter", 12, 0.5, 12.5);
    TH1F* fHistV0Statistics            = new TH1F("fHistV0Statistics", "V0 Counter;Stage;Counts",14,0.5,14.5);
    TH1F* fHistV0StatisticsMC          = new TH1F("fHistV0StatisticsMC", "V0 Counter MC;Stage;Counts",14,0.5,14.5);
    TH1F* fHistV0StatisticsSigmaMC     = new TH1F("fHistV0StatisticsSigmaMC", "V0 Counter Sigma MC;Stage;Counts",14,0.5,14.5);
    TH1F* fHistProtonStatistics        = new TH1F("fHistProtonStatistics", "Proton Counter;Stage;Counts",8,-0.5,7.5);
    TH1F* fHistProtonStatisticsMC      = new TH1F("fHistProtonStatisticsMC", "Proton Counter MC;Stage;Counts",8,-0.5,7.5);
    TH1F* fHistProtonStatisticsSigmaMC = new TH1F("fHistProtonStatisticsSigmaMC", "Proton Counter Sigma MC;Stage;Counts",8,-0.5,7.5);
    TH1F* fHistAddV0Statistics         = new TH1F("fHistAddV0Statistics", "Additional V0 Counter;Stage;Counts",11,0.5,11.5);
    TH1F* fHistAddV0StatisticsMC       = new TH1F("fHistAddV0StatisticsMC", "Additional V0 Counter MC;Stage;Counts",11,0.5,11.5);
    TH1F* fHistAddV0StatisticsSigmaMC  = new TH1F("fHistAddV0StatisticsSigmaMC", "Additional V0 Counter Sigma MC;Stage;Counts",11,0.5,11.5);
    TH1F* fHistGammaPairStats          = new TH1F("fHistGammaPairStats", "Gamma Pair Counter;;Counts",6,0.5,6.5);
    TH1F* fHistGammaPairStatsOneadd    = new TH1F("fHistGammaPairStatsOneadd", "Gamma Pair Counter, One Additional;;Counts",6,0.5,6.5);
    TH1F* fHistGammaPairStatsOnlyadd   = new TH1F("fHistGammaPairStatsOnlyadd", "Gamma Pair Counter, One Additional;;Counts",6,0.5,6.5);
    TH1F* fHistSigmaCounter            = new TH1F("fHistSigmaCounter", "#Sigma Counter", 7, 0.5, 7.5);
                                       
    //MC Information              
    TH2F* fHistMCPrimSigmaPtvsRap      = new TH2F("fHistMCPrimSigmaPtvsRap","Transverse momentum of primary MC #Sigma^{+} vs. y;#it{p}_{T} (GeV/#it{c});y",1200,0,12,3000,-1.5,1.5);
    TH2F* fHistMCPrimAntiSigmaPtvsRap  = new TH2F("fHistMCPrimAntiSigmaPtvsRap","Transverse momentum of primary MC #bar#Sigma^{-} vs. y;#it{p}_{T} (GeV/#it{c});y",1200,0,12,3000,-1.5,1.5);
    TH1F* fHistMCPrimSigmaPtRap05      = new TH1F("fHistMCPrimSigmaPtRap05","Transverse momentum of primary MC #Sigma^{+}, |y|<0.5;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPtRap05  = new TH1F("fHistMCPrimAntiSigmaPtRap05","Transverse momentum of primary MC #bar#Sigma^{-}, |y|<0.5;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPtRap08      = new TH1F("fHistMCPrimSigmaPtRap08","Transverse momentum of primary MC #Sigma^{+}, |y|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPtRap08  = new TH1F("fHistMCPrimAntiSigmaPtRap08","Transverse momentum of primary MC #bar#Sigma^{-}, |y|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPtRap10      = new TH1F("fHistMCPrimSigmaPtRap10","Transverse momentum of primary MC #Sigma^{+}, |y|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPtRap10  = new TH1F("fHistMCPrimAntiSigmaPtRap10","Transverse momentum of primary MC #bar#Sigma^{-}, |y|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH2F* fHistMCSigmaPtvsRap          = new TH2F("fHistMCSigmaPtvsRap","Transverse momentum of MC #Sigma^{+} vs. y;#it{p}_{T} (GeV/#it{c});y",1200,0,12,3000,-1.5,1.5);
    TH2F* fHistMCAntiSigmaPtvsRap      = new TH2F("fHistMCAntiSigmaPtvsRap","Transverse momentum of MC #bar#Sigma^{-} vs. y;#it{p}_{T} (GeV/#it{c});y",1200,0,12,3000,-1.5,1.5);
    TH1F* fHistMCSigmaPtRap05          = new TH1F("fHistMCSigmaPtRap05","Transverse momentum of MC #Sigma^{+}, |y|<0.5;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiSigmaPtRap05      = new TH1F("fHistMCAntiSigmaPtRap05","Transverse momentum of MC #bar#Sigma^{-}, |y|<0.5;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaPtRap08          = new TH1F("fHistMCSigmaPtRap08","Transverse momentum of MC #Sigma^{+}, |y|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiSigmaPtRap08      = new TH1F("fHistMCAntiSigmaPtRap08","Transverse momentum of MC #bar#Sigma^{-}, |y|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaPtRap10          = new TH1F("fHistMCSigmaPtRap10","Transverse momentum of MC #Sigma^{+}, |y|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiSigmaPtRap10      = new TH1F("fHistMCAntiSigmaPtRap10","Transverse momentum of MC #bar#Sigma^{-}, |y|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaPt               = new TH1F("fHistMCSigmaPt","Transverse momentum of MC #Sigma^{+};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiSigmaPt           = new TH1F("fHistMCAntiSigmaPt","Transverse momentum of MC #bar#Sigma^{-};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPt           = new TH1F("fHistMCPrimSigmaPt","Transverse momentum of primary MC #Sigma^{+};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPt       = new TH1F("fHistMCPrimAntiSigmaPt","Transverse momentum of primary MC #bar#Sigma^{-};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPt08         = new TH1F("fHistMCPrimSigmaPt08","Transverse momentum of primary MC #Sigma^{+}, |#eta|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPt08     = new TH1F("fHistMCPrimAntiSigmaPt08","Transverse momentum of primary MC #bar#Sigma^{-}, |#eta|<0.8;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPt10         = new TH1F("fHistMCPrimSigmaPt10","Transverse momentum of primary MC #Sigma^{+}, |#eta|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPt10     = new TH1F("fHistMCPrimAntiSigmaPt10","Transverse momentum of primary MC #bar#Sigma^{-}, |#eta|<1;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH2F* fHistMCPrimSigmaPtvsEta      = new TH2F("fHistMCPrimSigmaPtvsEta","Transverse momentum of primary MC #Sigma^{+} vs. |#eta|;#it{p}_{T} (GeV/#it{c});|#eta|",1200,0,12,30,0,1.5);
    TH2F* fHistMCPrimAntiSigmaPtvsEta  = new TH2F("fHistMCPrimAntiSigmaPtvsEta","Transverse momentum of primary MC #bar#Sigma^{-} vs. |#eta|;#it{p}_{T} (GeV/#it{c});|#eta|",1200,0,12,30,0,1.5);
    TH2F* fHistMCSigmaPtvsEta          = new TH2F("fHistMCSigmaPtvsEta","Transverse momentum of MC #Sigma^{+} vs. |#eta|;#it{p}_{T} (GeV/#it{c});|#eta|",1200,0,12,30,0,1.5);
    TH2F* fHistMCAntiSigmaPtvsEta      = new TH2F("fHistMCAntiSigmaPtvsEta","Transverse momentum of MC #bar#Sigma^{-} vs. |#eta|;#it{p}_{T} (GeV/#it{c});|#eta|",1200,0,12,30,0,1.5);
    TH1F* fHistMCSigmaOrigin           = new TH1F("fHistMCSigmaOrigin","Origin of MC #Sigma^{+};r_{xy} (cm);Counts/mm",5000,0,500);   
    TH1F* fHistMCAntiSigmaOrigin       = new TH1F("fHistMCAntiSigmaOrigin","Origin of MC #bar#Sigma^{-};r_{xy} (cm);Counts/mm",5000,0,500);   
    TH1F* fHistMCDeltaPt               = new TH1F("fHistMCDeltaPt","Transverse momentum of MC #Delta^{+};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiDeltaPt           = new TH1F("fHistMCAntiDeltaPt","Transverse momentum of MC #bar#Delta^{-};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPi0Pt                 = new TH1F("fHistMCPi0Pt","Transverse momentum of MC #pi^{0};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCProtonPt              = new TH1F("fHistMCProtonPt","Transverse momentum of MC Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiProtonPt          = new TH1F("fHistMCAntiProtonPt","Transverse momentum of MC Anti-Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimProtonPt          = new TH1F("fHistMCPrimProtonPt","Transverse momentum primary of MC Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiProtonPt      = new TH1F("fHistMCPrimAntiProtonPt","Transverse momentum of primary MC Anti-Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPhotonPt              = new TH1F("fHistMCPhotonPt","Transverse momentum of MC Photons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaPhotonPt         = new TH1F("fHistMCSigmaPhotonPt","Transverse momentum of MC Photons from #Sigma^{+}/#bar#Sigma^{-};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCConvPhotonPt          = new TH1F("fHistMCConvPhotonPt","Transverse momentum of converted MC Photons (R<180cm);#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaConvPhotonPt     = new TH1F("fHistMCSigmaConvPhotonPt","Transverse momentum of converted MC Photons from #Sigma^{+}/#bar#Sigma^{-} (R<180cm);#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCConvRadius            = new TH1F("fHistMCConvRadius","Conversion radius of MC Photons;r_{xy} (cm);Counts/mm",5000,0,500);   
    TH2F* fHistMCConvRadiusvspt        = new TH2F("fHistMCConvRadiusvspt","Conversion radius of MC Photons vs. #it{p}_{T};r_{xy} (cm);#it{p}_{T} (GeV/#it{c})",500,0,500,500,0,5);   
    TH1F* fHistSigmaMotherPart         = new TH1F("fHistSigmaMotherPart", "#Sigma^{+} Mother Particle;Mother Particle;Counts", 25, 0.5, 25.5); 
    TH1F* fHistAntiSigmaMotherPart     = new TH1F("fHistAntiSigmaMotherPart", "#bar#Sigma^{-} Mother Particle;Mother Particle;Counts", 25, 0.5, 25.5); 

    //Track Quality                        
    TH2F* fHistTrackEtaPhi             = new TH2F("fHistTrackEtaPhi","#eta vs. #phi of Tracks;#eta;#phi (rad)",300,-1.5,1.5,300,0,2*TMath::Pi());                                
    TH1F* fHistTrackChi2               = new TH1F("fHistTrackChi2","#chi^{2}/NDF of Tracks;#chi^{2}/NDF",200,0,100);
    TH1F* fHistTrackTPCCluster         = new TH1F("fHistTrackTPCCluster","Number of TPC Clusters of Tracks",161,-0.5,160.5);
    TH1F* fHistTrackITSCluster         = new TH1F("fHistTrackITSCluster","Number of ITS Clusters of Tracks",13,-0.5,12.5);
    TH2F* fHistTrackpvsdEdx            = new TH2F("fHistTrackpvsdEdx","Momentum vs. TPC dE/dx;p (GeV/#it{c});TPC #frac{dE}{dx} (a.u.)",500,0,10,500,0,500);
    TH2F* fHistTrackpvsbeta            = new TH2F("fHistTrackpvsbeta","Momentum vs. #beta;p (GeV/#it{c});#beta",500,0,10,500,0,1.2);
    TH1F* fHistTrackpt                 = new TH1F("fHistTrackpt","Transverse momentum of Tracks;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);

    //Proton QA                        
    TH2F* fHistProtonEtaPhiMC          = new TH2F("fHistProtonEtaPhiMC","#eta vs. #phi of Protons MC;#eta;#phi (rad)",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistProtonChi2MC            = new TH1F("fHistProtonChi2MC","#chi^{2}/NDF of Protons MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistProtonTPCClusterMC      = new TH1F("fHistProtonTPCClusterMC","Number of TPC Clusters of Protons MC",161,-0.5,160.5);
    TH1F* fHistProtonITSClusterMC      = new TH1F("fHistProtonITSClusterMC","Number of ITS Clusters of Protons MC",13,-0.5,12.5);
    TH2F* fHistProtonpvsNSigmaTPC      = new TH2F("fHistProtonpvsNSigmaTPC","Momentum vs N sigma TPC Proton;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTPCMC    = new TH2F("fHistProtonpvsNSigmaTPCMC","Momentum vs N sigma TPC Proton MC;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTOF      = new TH2F("fHistProtonpvsNSigmaTOF","Momentum vs N sigma TOF Proton;p (GeV/#it{c});n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTOFMC    = new TH2F("fHistProtonpvsNSigmaTOFMC","Momentum vs N sigma TOF Proton MC;p (GeV/#it{c});n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH1F* fHistProtonDCAxy             = new TH1F("fHistProtonDCAxy","DCA_{xy} of Protons;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAz              = new TH1F("fHistProtonDCAz","DCA_{z} of Protons;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAxyMC           = new TH1F("fHistProtonDCAxyMC","DCA_{xy} of Protons MC;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAzMC            = new TH1F("fHistProtonDCAzMC","DCA_{z} of Protons MC;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAxyMCSigma      = new TH1F("fHistProtonDCAxyMCSigma","DCA_{xy} of Protons from MC #Sigma^{+}/#bar#Sigma^{-};DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAzMCSigma       = new TH1F("fHistProtonDCAzMCSigma","DCA_{z} of Protons from MC #Sigma^{+}/#bar#Sigma^{-};DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAxyMCPrimSig    = new TH1F("fHistProtonDCAxyMCPrimSig","DCA_{xy} of Protons from primary MC #Sigma^{+}/#bar#Sigma^{-};DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonDCAzMCPrimSig     = new TH1F("fHistProtonDCAzMCPrimSig","DCA_{z} of Protons from primary MC #Sigma^{+}/#bar#Sigma^{-};DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistPrimProtonDCAxyMC       = new TH1F("fHistPrimProtonDCAxyMC","DCA_{xy} of primary Protons MC;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistPrimProtonDCAzMC        = new TH1F("fHistPrimProtonDCAzMC","DCA_{z} of primary Protons MC;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistMaterialProtonDCAxyMC   = new TH1F("fHistMaterialProtonDCAxyMC","DCA_{xy} of Protons from Material MC;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistMaterialProtonDCAzMC    = new TH1F("fHistMaterialProtonDCAzMC","DCA_{z} of Protons from Material MC;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistWeakProtonDCAxyMC       = new TH1F("fHistWeakProtonDCAxyMC","DCA_{xy} of Protons from weak decay MC;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistWeakProtonDCAzMC        = new TH1F("fHistWeakProtonDCAzMC","DCA_{z} of Protons from weak decay MC;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistLambdaProtonDCAxyMC     = new TH1F("fHistLambdaProtonDCAxyMC","DCA_{xy} of Protons from #Lambda MC;DCA_{xy} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistLambdaProtonDCAzMC      = new TH1F("fHistLambdaProtonDCAzMC","DCA_{z} of Protons from #Lambda MC;DCA_{z} (cm);Counts/(10 #mum))",10000,0,10);
    TH1F* fHistProtonptMC              = new TH1F("fHistProtonptMC","Transverse momentum of MC Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistProtonptwCutsMC         = new TH1F("fHistProtonptwCutsMC","Transverse momentum of MC Protons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistProtonptwCuts           = new TH1F("fHistProtonptwCuts","Transverse momentum of Protons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);   

    //V0 QA                            
    TH2F* fHistV0OnflyvsOffline        = new TH2F("fHistV0OnflyvsOffline", "Number of V0s found by On-the-fly Finder vs. Offline Finder", 201, -0.5, 200.5, 201, -0.5, 200.5);
    TH2F* fHistV0DaughtEtaPhi          = new TH2F("fHistV0DaughtEtaPhi","#eta vs. #phi of V0 Daughters;#eta;#phi (rad)",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH2F* fHistV0DaughtEtaPhiMC        = new TH2F("fHistV0DaughtEtaPhiMC","#eta vs. #phi of V0 Daughters MC;#eta;#phi (rad)",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistV0DaughtChi2            = new TH1F("fHistV0DaughtChi2","#chi^{2}/NDF of V0 Daughters;#chi^{2}/NDF",200,0,100);
    TH1F* fHistV0DaughtChi2MC          = new TH1F("fHistV0DaughtChi2MC","#chi^{2}/NDF of V0 Daughters MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistV0DaughtTPCClust        = new TH1F("fHistV0DaughtTPCClust","Number of TPC Clusters of V0 Daughters",161,-0.5,160.5);
    TH1F* fHistV0DaughtITSClust        = new TH1F("fHistV0DaughtITSClust","Number of ITS Clusters of V0 Daughters",13,-0.5,12.5);
    TH1F* fHistV0DaughtTPCClustMC      = new TH1F("fHistV0DaughtTPCClustMC","Number of TPC Clusters of V0 Daughters MC",161,-0.5,160.5);
    TH1F* fHistV0DaughtITSClustMC      = new TH1F("fHistV0DaughtITSClustMC","Number of ITS Clusters of V0 Daughters MC",13,-0.5,12.5);
    TH2F* fHistV0DaughtpvsNSigmaTPC    = new TH2F("fHistV0DaughtpvsNSigmaTPC","Momentum vs N sigma TPC Electron;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH1F* fHistV0DaughtDCAtoPV         = new TH1F("fHistV0DaughtDCAtoPV","DCA to PV of V0 Daughters;DCA (cm);Counts/mm",200,0,20);
    TH1F* fHistV0DaughtDCAtoPVMC       = new TH1F("fHistV0DaughtDCAtoPVMC","DCA to PV of V0 Daughters MC;DCA (cm);Counts/mm",200,0,20);
    TH1F* fHistV0DaughtDCA             = new TH1F("fHistV0DaughtDCA","DCA between V0 Daughters;DCA (cm);Counts/(0.5mm)",200,0,10);
    TH1F* fHistV0DaughtDCAMC           = new TH1F("fHistV0DaughtDCAMC","DCA between V0 Daughters MC;DCA (cm);Counts/(0.5mm)",200,0,10);
    TH1F* fHistV0CPA                   = new TH1F("fHistV0CPA","Cosine of Pointing Angle of V0s;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0CPAMC                 = new TH1F("fHistV0CPAMC","Cosine of Pointing Angle of V0s MC;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0CPAMCSigma            = new TH1F("fHistV0CPAMCSigma","Cosine of Pointing Angle of V0s MC (#gamma from #Sigma);cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0Radius                = new TH1F("fHistV0Radius","Radius of V0s;r_{xy} (cm);Counts/mm",3000,0,300);   
    TH1F* fHistV0RadiusMC              = new TH1F("fHistV0RadiusMC","Radius of V0s MC;r_{xy} (cm);Counts/mm",3000,0,300);   
    TH2F* fHistV0Position2D            = new TH2F("fHistV0Position2D","Position of V0s;x (cm);y (cm)",1000,-180,180,1000,-180,180);   
    TH2F* fHistV0Position2DMC          = new TH2F("fHistV0Position2DMC","Position of V0s MC;x (cm);y (cm)",1000,-180,180,1000,-180,180);   
    TH1F* fHistV0PhotonDCAPV           = new TH1F("fHistV0PhotonDCAPV","DCA to PV of Photons;DCA (cm);Counts/mm)",100,0,10);
    TH1F* fHistV0PhotonDCAPVMC         = new TH1F("fHistV0PhotonDCAPVMC","DCA to PV of MC Photons;DCA (cm);Counts/mm)",100,0,10);
    TH1F* fHistV0PhotonDCAPVMCSigma    = new TH1F("fHistV0PhotonDCAPVMCSigma","DCA to PV of MC Photons from #Sigma^{+}/#bar#Sigma^{-};DCA (cm);Counts/mm)",100,0,10);
    TH2F* fHistV0ArmPod                = new TH2F("fHistV0ArmPod", "Armenteros-Podolanski Plot of V0s;#alpha;q_{t} (GeV/#it{c})",300,-1.5,1.5,200,-0.05,0.4);       
    TH2F* fHistV0ArmPodMC              = new TH2F("fHistV0ArmPodMC", "Armenteros-Podolanski Plot of V0s MC;#alpha;q_{t} (GeV/#it{c})",300,-1.5,1.5,200,-0.05,0.4);
    TH1F* fHistV0OpenAngle             = new TH1F("fHistV0OpenAngle","Total Opening Angle of V0s;#xi (rad);Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistV0OpenAngleMC           = new TH1F("fHistV0OpenAngleMC","Total Opening Angle of V0s MC;#xi (rad);Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistV0DeltaTheta            = new TH1F("fHistV0DeltaTheta","#Delta#Theta of V0 Daughters;#Delta#Theta (rad);Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistV0DeltaThetaMC          = new TH1F("fHistV0DeltaThetaMC","#Delta#Theta of V0 Daughters MC;#Delta#Theta (rad);Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistV0InvMass               = new TH1F("fHistV0InvMass", "Invariant mass of V0s;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistV0InvMassMC             = new TH1F("fHistV0InvMassMC", "Invariant mass of V0s MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistV0ptMC                  = new TH1F("fHistV0ptMC","Transverse momentum of MC V0 Photons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistV0ptwCutsMC             = new TH1F("fHistV0ptwCutsMC","Transverse momentum of MC V0 Photons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistV0SigmaptMC             = new TH1F("fHistV0SigmaptMC","Transverse momentum of MC V0 Photons from #Sigma^{+}/#bar#Sigma^{-};#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistV0SigmaptwCutsMC        = new TH1F("fHistV0SigmaptwCutsMC","Transverse momentum of MC V0 Photons from #Sigma^{+}/#bar#Sigma^{-} with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistV0ptwCuts               = new TH1F("fHistV0ptwCuts","Transverse momentum of V0 Photons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);    

    //Gamma Gamma QA                      
    TH1F* fHistGammaPairInvMass        = new TH1F("fHistGammaPairInvMass", "Invariant mass of Photon Pairs;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMC      = new TH1F("fHistGammaPairInvMassMC", "Invariant mass of Photon Pairs MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnfly   = new TH1F("fHistGammaPairInvMassOnfly", "Invariant mass of On-the-fly Photon Pairs;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMCOnfly = new TH1F("fHistGammaPairInvMassMCOnfly", "Invariant mass of On-the-fly Photon Pairs MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAdd  = new TH1F("fHistGammaPairInvMassOneAdd", "Invariant mass of Photon Pairs, one additional Photon;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAddMC= new TH1F("fHistGammaPairInvMassOneAddMC", "Invariant mass of On-fly Photon Pairs, one additional Photon MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAdd = new TH1F("fHistGammaPairInvMassOnlyAdd", "Invariant mass of Photon Pairs, only additional Photons;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAddMC= new TH1F("fHistGammaPairInvMassOnlyAddMC", "Invariant mass of Photon Pairs, only additional Photons MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairDCA            = new TH1F("fHistGammaPairDCA","DCA between Photon Pairs;DCA (cm);Counts/(0.5mm)",200,0,10); 
    TH1F* fHistGammaPairDCAMC          = new TH1F("fHistGammaPairDCAMC","DCA between Photon Pairs MC;DCA (cm);Counts/(0.5mm)",200,0,10);
    TH1F* fHistGammaPairInvMassPHOS    = new TH1F("fHistGammaPairInvMassPHOS","Invariant mass of Photon Pairs;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMCPHOS  = new TH1F("fHistGammaPairInvMassMCPHOS","Invariant mass of Photon Pairs MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnflyPHOS = new TH1F("fHistGammaPairInvMassOnflyPHOS","Invariant mass of On-the-fly Photon Pairs;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMCOnflyPHOS = new TH1F("fHistGammaPairInvMassMCOnflyPHOS","Invariant mass of On-the-fly Photon Pairs MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);

    //Gamma Gamma QA - check auto correlations                      
    TH1F* fHistGammaPairInvMass2       = new TH1F("fHistGammaPairInvMass2", "Invariant mass of Photon Pairs (check auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnfly2  = new TH1F("fHistGammaPairInvMassOnfly2", "Invariant mass of On-the-fly Photon Pairs (check auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAdd2 = new TH1F("fHistGammaPairInvMassOneAdd2", "Invariant mass of Photon Pairs, one additional Photon (check auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAdd2= new TH1F("fHistGammaPairInvMassOnlyAdd2", "Invariant mass of Photon Pairs, only additional Photons (check auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairDCA2           = new TH1F("fHistGammaPairDCA2","DCA between Photon Pairs (check auto-corr);DCA (cm);Counts/(0.5mm)",200,0,10); 

    //Gamma Gamma QA - likely auto correlations                      
    TH1F* fHistGammaPairInvMass3       = new TH1F("fHistGammaPairInvMass3", "Invariant mass of Photon Pairs (likely auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnfly3  = new TH1F("fHistGammaPairInvMassOnfly3", "Invariant mass of On-the-fly Photon Pairs (likely auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAdd3 = new TH1F("fHistGammaPairInvMassOneAdd3", "Invariant mass of Photon Pairs, one additional Photon (likely auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAdd3= new TH1F("fHistGammaPairInvMassOnlyAdd3", "Invariant mass of Photon Pairs, only additional Photons (likely auto-corr);m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairDCA3           = new TH1F("fHistGammaPairDCA3","DCA between Photon Pairs (likely auto-corr);DCA (cm);Counts/(0.5mm)",200,0,10); 

    //MC Sigma Topology
    TH1F* fHistKFSigmaVertexResX       = new TH1F("fHistKFSigmaVertexResX","KF #Sigma decay vertex X - MC #Sigma decay vertex X;#DeltaX (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistKFSigmaVertexResY       = new TH1F("fHistKFSigmaVertexResY","KF #Sigma decay vertex Y - MC #Sigma decay vertex Y;#DeltaY (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistKFSigmaVertexResZ       = new TH1F("fHistKFSigmaVertexResZ","KF #Sigma decay vertex Z - MC #Sigma decay vertex Z;#DeltaZ (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistPi0VertexvsMC           = new TH1F("fHistPi0VertexvsMC","#pi^{0} KF Decay Radius - MC Decay Radius;#Deltar (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistMCSigmaPA               = new TH1F("fHistMCSigmaPA","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistMCSigmaPArot            = new TH1F("fHistMCSigmaPArot","#Sigma^{+}/#bar#Sigma^{-} PA rotated;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistMCPrimSigmaPA           = new TH1F("fHistMCPrimSigmaPA","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistMCPrimSigmaPArot        = new TH1F("fHistMCPrimSigmaPArot","#Sigma^{+}/#bar#Sigma^{-} PA rotated;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistSigmaPA                 = new TH1F("fHistSigmaPA","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaPAmix              = new TH1F("fHistSigmaPAmix","#Sigma^{+}/#bar#Sigma^{-} PA mixed;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaPArot              = new TH1F("fHistSigmaPArot","#Sigma^{+}/#bar#Sigma^{-} PA rotated;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaY                  = new TH1F("fHistSigmaY","#Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistMCSigmaY                = new TH1F("fHistMCSigmaY","#Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistSigmaYrot               = new TH1F("fHistSigmaYrot","#Sigma^{+}/#bar#Sigma^{-} Y rotated;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistMCSigmaYrot             = new TH1F("fHistMCSigmaYrot","#Sigma^{+}/#bar#Sigma^{-} Y rotated;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistSigmaYmix               = new TH1F("fHistSigmaYmix","#Sigma^{+}/#bar#Sigma^{-} Y mixed;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistInvSigmaMass            = new TH1F("fHistInvSigmaMass","Invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistMCInvSigmaMass          = new TH1F("fHistMCInvSigmaMass","Invariant mass of MC #Sigma^{+}/#bar#Sigma^{-};m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 200, 1.1, 1.3);
    TH1F* fHistMCInvSigmaMassrot       = new TH1F("fHistMCInvSigmaMassrot","Invariant mass of MC #Sigma^{+}/#bar#Sigma^{-};m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 200, 1.1, 1.3);
    TH1F* fHistInvSigmaMassmix         = new TH1F("fHistInvSigmaMassmix","Invariant mass of mixed #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistMCOneGammaSigmaPA       = new TH1F("fHistMCOneGammaSigmaPA","#Sigma^{+} PA, One #gamma MC;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistMCPrimOneGammaSigmaPA   = new TH1F("fHistMCPrimOneGammaSigmaPA","#Sigma^{+} PA, One #gamma primary MC;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistOneGammaSigmaPA         = new TH1F("fHistOneGammaSigmaPA","#Sigma^{+} PA, One #gamma;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistPi0VertexMC             = new TH1F("fHistPi0VertexMC","#pi^{0} MC Decay Radius;r (cm);Counts/(mm)", 200, 0, 20);
    TH1F* fHistKFSigmaVertexResXPHOS   = new TH1F("fHistKFSigmaVertexResXPHOS","KF #Sigma decay vertex X - MC #Sigma decay vertex X;#DeltaX (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistKFSigmaVertexResYPHOS   = new TH1F("fHistKFSigmaVertexResYPHOS","KF #Sigma decay vertex Y - MC #Sigma decay vertex Y;#DeltaY (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistKFSigmaVertexResZPHOS   = new TH1F("fHistKFSigmaVertexResZPHOS","KF #Sigma decay vertex Z - MC #Sigma decay vertex Z;#DeltaZ (cm);Counts/(mm)", 400, -20, 20);
    TH1F* fHistInvSigmaMassPHOS        = new TH1F("fHistInvSigmaMassPHOS","Invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistMCInvSigmaMassPHOS      = new TH1F("fHistMCInvSigmaMassPHOS","Invariant mass of MC #Sigma^{+}/#bar#Sigma^{-};m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 200, 1.1, 1.3);
    TH1F* fHistMCInvSigmaMassPHOSrot   = new TH1F("fHistMCInvSigmaMassPHOSrot","Invariant mass of MC #Sigma^{+}/#bar#Sigma^{-};m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 200, 1.1, 1.3);
    TH1F* fHistInvSigmaMassPHOSUncorr   = new TH1F("fHistInvSigmaMassPHOSUncorr","Uncorrected invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 400, 1, 1.4);
    TH1F* fHistMCInvSigmaMassPHOSUncorr = new TH1F("fHistMCInvSigmaMassPHOSUncorr","Uncorrected invariant mass of MC #Sigma^{+}/#bar#Sigma^{-};m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 200, 1.1, 1.3);
    TH1F* fHistSigmaYPHOS              = new TH1F("fHistSigmaYPHOS","#Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistMCSigmaYPHOS            = new TH1F("fHistMCSigmaYPHOS","#Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistMCSigmaPAPHOS           = new TH1F("fHistMCSigmaPAPHOS","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.001)",3000,0,3);
    TH1F* fHistMCSigmaAntiPAPHOS       = new TH1F("fHistMCSigmaAntiPAPHOS","#Sigma^{+}/#bar#Sigma^{-} Anti-PA;PA (rad);Counts/(0.001)",3000,0,3);
    TH1F* fHistMCPrimSigmaPAPHOS       = new TH1F("fHistMCPrimSigmaPAPHOS","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.001)",3000,0,3);
    TH1F* fHistMCPrimSigmaPAPHOSrot    = new TH1F("fHistMCPrimSigmaPAPHOSrot","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.001)",3000,0,3);
    TH1F* fHistMCPrimSigmaAntiPAPHOS   = new TH1F("fHistMCPrimSigmaAntiPAPHOS","#Sigma^{+}/#bar#Sigma^{-} Anti-PA;PA (rad);Counts/(0.001)",3000,0,3);
    TH1F* fHistSigmaPAPHOS             = new TH1F("fHistSigmaPAPHOS","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaPAPHOSrot          = new TH1F("fHistSigmaPAPHOSrot","#Sigma^{+}/#bar#Sigma^{-} PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaAntiPAPHOS         = new TH1F("fHistSigmaAntiPAPHOS","#Sigma^{+}/#bar#Sigma^{-} Anti-PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistInvSigmaMassPHOSmix     = new TH1F("fHistInvSigmaMassPHOSmix","Invariant mass of mixed #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistSigmaYPHOSmix           = new TH1F("fHistSigmaYPHOSmix","#Sigma^{+}/#bar#Sigma^{-} Y mixed;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistSigmaPAPHOSmix          = new TH1F("fHistSigmaPAPHOSmix","#Sigma^{+}/#bar#Sigma^{-} PA mixed;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaAntiPAPHOSmix      = new TH1F("fHistSigmaAntiPAPHOSmix","#Sigma^{+}/#bar#Sigma^{-} Anti-PA mixed;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistProtPhotonDCAPHOS       = new TH1F("fHistProtPhotonDCAPHOS","#Sigma^{+}/#bar#Sigma^{-} Photon Proton DCA KF;DCA (cm);Counts/(100 #mum)",500,0,5);
    TH1F* fHistProtPhotonDCAPHOSMC     = new TH1F("fHistProtPhotonDCAPHOSMC","MC #Sigma^{+}/#bar#Sigma^{-} Photon Proton DCA KF;DCA (cm);Counts/(100 #mum)",500,0,5);
    TH1F* fHistSigmaDCAtoPVPHOS        = new TH1F("fHistSigmaDCAtoPVPHOS","#Sigma^{+}/#bar#Sigma^{-} Dist to PV KF;DCA (cm);Counts/(mm)",500,0,50);
    TH1F* fHistSigmaRadiusPHOS         = new TH1F("fHistSigmaRadiusPHOS","#Sigma^{+}/#bar#Sigma^{-} Radius KF;R (cm);Counts/(mm)",500,0,50);
    TH1F* fHistSigmaDCAtoPVPHOSMC      = new TH1F("fHistSigmaDCAtoPVPHOSMC","MC #Sigma^{+}/#bar#Sigma^{-} Dist to PV KF;DCA (cm);Counts/(mm)",500,0,50);
    TH1F* fHistPhotonSecPAPHOS         = new TH1F("fHistPhotonSecPAPHOS","Photon Secondary PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistPhotonSecPAPHOSMC       = new TH1F("fHistPhotonSecPAPHOSMC","Photon Secondary PA MC;PA (rad);Counts/(0.005)",700,0,3.5);

    //Sigma Momentum Resolution
    TH1F* fHistSigmaPxResnoprop        = new TH1F("fHistSigmaPxResnoprop","#Sigma p_{x} - MC #Sigma p_{x}, no propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResnoprop        = new TH1F("fHistSigmaPyResnoprop","#Sigma p_{y} - MC #Sigma p_{y}, no propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResnoprop        = new TH1F("fHistSigmaPzResnoprop","#Sigma p_{z} - MC #Sigma p_{z}, no propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxResprop          = new TH1F("fHistSigmaPxResprop","#Sigma p_{x} - MC #Sigma p_{x}, ETP propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResprop          = new TH1F("fHistSigmaPyResprop","#Sigma p_{y} - MC #Sigma p_{y}, ETP propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResprop          = new TH1F("fHistSigmaPzResprop","#Sigma p_{z} - MC #Sigma p_{z}, ETP propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxResnopropPHOS    = new TH1F("fHistSigmaPxResnopropPHOS","#Sigma p_{x} - MC #Sigma p_{x}, no propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResnopropPHOS    = new TH1F("fHistSigmaPyResnopropPHOS","#Sigma p_{y} - MC #Sigma p_{y}, no propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResnopropPHOS    = new TH1F("fHistSigmaPzResnopropPHOS","#Sigma p_{z} - MC #Sigma p_{z}, no propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxRespropPHOS      = new TH1F("fHistSigmaPxRespropPHOS","#Sigma p_{x} - MC #Sigma p_{x}, ETP propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyRespropPHOS      = new TH1F("fHistSigmaPyRespropPHOS","#Sigma p_{y} - MC #Sigma p_{y}, ETP propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzRespropPHOS      = new TH1F("fHistSigmaPzRespropPHOS","#Sigma p_{z} - MC #Sigma p_{z}, ETP propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);

    //Event Mixing
    TH3F* fHistPairNMixedEvents        = new TH3F("fHistPairNMixedEvents","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins2, fMinCentBin2, fMaxCentBin2, nZvtxBins2, fMinZBin2, fMaxZBin2, fEvPoolSize2+1, -0.5, fEvPoolSize2+0.5);
    TH3F* fHistBkgNMixedEvents         = new TH3F("fHistBkgNMixedEvents","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins, fMinCentBin, fMaxCentBin, nZvtxBins, fMinZBin, fMaxZBin, fEvPoolSize+1, -0.5, fEvPoolSize+0.5);
    TH3F* fHistPairNMixedEventsPHOS    = new TH3F("fHistPairNMixedEventsPHOS","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins2, fMinCentBin2, fMaxCentBin2, nZvtxBins2, fMinZBin2, fMaxZBin2, fEvPoolSize2+1, -0.5, fEvPoolSize2+0.5);
    TH3F* fHistBkgNMixedEventsPHOS     = new TH3F("fHistBkgNMixedEventsPHOS","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins, fMinCentBin, fMaxCentBin, nZvtxBins, fMinZBin, fMaxZBin, fEvPoolSize+1, -0.5, fEvPoolSize+0.5);

    //KF Checks
    TH1F* fHistPhotonKFCheck           = new TH1F("fHistPhotonKFCheck","Check Photon Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistPi0KFCheck              = new TH1F("fHistPi0KFCheck","Check Pion Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistSigmaKFCheck            = new TH1F("fHistSigmaKFCheck","Check Sigma Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistPhotonKFCheckPHOS       = new TH1F("fHistPhotonKFCheckPHOS","Check Photon Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistSigmaKFCheckPHOS        = new TH1F("fHistSigmaKFCheckPHOS","Check Sigma Daughters to avoid floating point exceptions", 4, 0.5, 4.5);

    //PHOS
    TH1F* fHistNCaloPhotons            = new TH1F("fHistNCaloPhotons","Number of Clusters per Event;;Counts",1000,0,1000);
    TH1F* fHistNClusters               = new TH1F("fHistNClusters","Number of selected Cluster Photons per Event;;Counts",1000,0,1000);
    TH1F* fHistClusterStatistics       = new TH1F("fHistClusterStatistics", "Cluster Counter;Stage;Counts",5,0.5,5.5);
    TH1F* fHistClusterStatisticsMC     = new TH1F("fHistClusterStatisticsMC", "Cluster Counter MC;Stage;Counts",5,0.5,5.5);
    TH1F* fHistClusterStatisticsMCSig  = new TH1F("fHistClusterStatisticsMCSig", "Cluster Counter Sigma MC;Stage;Counts",5,0.5,5.5);
    TH1F* fHistClusterType             = new TH1F("fHistClusterType","Cluster Type;Type;Counts",3,0.5,3.5);
    TH1F* fHistClusterTypeMC           = new TH1F("fHistClusterTypeMC","Cluster Type MC;Type;Counts",3,0.5,3.5);
    TH1F* fHistClusterTypeMCSig        = new TH1F("fHistClusterTypeMCSig","Cluster Type Sigma MC;Type;Counts",3,0.5,3.5);
    TH1F* fHistPHOSDisttoBC            = new TH1F("fHistPHOSDisttoBC","Distance to Bad Channel PHOS;Distance;Counts/0.1",200,0,20);
    TH1F* fHistPHOSM02                 = new TH1F("fHistPHOSM02","M02 PHOS;M02;Counts/0.1",200,0,20);
    TH1F* fHistPHOSM20                 = new TH1F("fHistPHOSM20","M20 PHOS;M20;Counts/0.1",200,0,20);
    TH1F* fHistPHOSTOF                 = new TH1F("fHistPHOSTOF","TOF PHOS;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistPHOSBeta                = new TH1F("fHistPHOSBeta","#beta PHOS;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistPHOSE                   = new TH1F("fHistPHOSE","Energy PHOS;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistPHOSNTracks             = new TH1F("fHistPHOSNTracks","N Tracks matched to Cluster PHOS;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistPHOSNCells              = new TH1F("fHistPHOSNCells","N Cells of Cluster PHOS;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistPHOSDx                  = new TH1F("fHistPHOSDx","Dx (Closest Track in #phi) PHOS;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSDz                  = new TH1F("fHistPHOSDz","Dz (Closest Track in z (#eta)) PHOS;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSDy                  = new TH1F("fHistPHOSDy","Dy (Closest Track in #eta#phi) PHOS;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSDisp                = new TH1F("fHistPHOSDisp","Dispersion PHOS;Dispersion;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCDisttoBC          = new TH1F("fHistPHOSMCDisttoBC","Distance to Bad Channel PHOS MC;Distance;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCM02               = new TH1F("fHistPHOSMCM02","M02 PHOS MC;M02;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCM20               = new TH1F("fHistPHOSMCM20","M20 PHOS MC;M20;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCTOF               = new TH1F("fHistPHOSMCTOF","TOF PHOS MC;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistPHOSMCBeta              = new TH1F("fHistPHOSMCBeta","#beta PHOS MC;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistPHOSMCE                 = new TH1F("fHistPHOSMCE","Energy PHOS MC;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistPHOSMCNTracks           = new TH1F("fHistPHOSMCNTracks","N Tracks matched to Cluster PHOS MC;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistPHOSMCNCells            = new TH1F("fHistPHOSMCNCells","N Cells of Cluster PHOS MC;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistPHOSMCDx                = new TH1F("fHistPHOSMCDx","Dx (Closest Track in #phi) PHOS MC;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCDz                = new TH1F("fHistPHOSMCDz","Dz (Closest Track in z (#eta)) PHOS MC;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCDy                = new TH1F("fHistPHOSMCDy","Dy (Closest Track in #eta#phi) PHOS MC;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCDisp              = new TH1F("fHistPHOSMCDisp","Dispersion PHOS MC;Dispersion;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCSigDisttoBC       = new TH1F("fHistPHOSMCSigDisttoBC","Distance to Bad Channel PHOS MC Sigma;Distance;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCSigM02            = new TH1F("fHistPHOSMCSigM02","M02 PHOS MC Sigma;M02;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCSigM20            = new TH1F("fHistPHOSMCSigM20","M20 PHOS MC Sigma;M20;Counts/0.1",200,0,20);
    TH1F* fHistPHOSMCSigTOF            = new TH1F("fHistPHOSMCSigTOF","TOF PHOS MC Sigma;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistPHOSMCSigBeta           = new TH1F("fHistPHOSMCSigBeta","#beta PHOS MC Sigma;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistPHOSMCSigE              = new TH1F("fHistPHOSMCSigE","Energy PHOS MC Sigma;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistPHOSMCSigNTracks        = new TH1F("fHistPHOSMCSigNTracks","N Tracks matched to Cluster PHOS MC Sigma;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistPHOSMCSigNCells         = new TH1F("fHistPHOSMCSigNCells","N Cells of Cluster PHOS MC Sigma;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistPHOSMCSigDx             = new TH1F("fHistPHOSMCSigDx","Dx (Closest Track in #phi) PHOS MC Sigma;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCSigDz             = new TH1F("fHistPHOSMCSigDz","Dz (Closest Track in z (#eta)) PHOS MC Sigma;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCSigDy             = new TH1F("fHistPHOSMCSigDy","Dy (Closest Track in #eta#phi) PHOS MC Sigma;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistPHOSMCSigDisp           = new TH1F("fHistPHOSMCSigDisp","Dispersion PHOS MC Sigma;Dispersion;Counts/0.1",200,0,20);
    TH1F* fHistEMCALDisttoBC           = new TH1F("fHistEMCALDisttoBC","Distance to Bad Channel EMCAL;Distance;Counts/0.1",200,0,20);
    TH1F* fHistEMCALM02                = new TH1F("fHistEMCALM02","M02 EMCAL;M02;Counts/0.1",200,0,20);
    TH1F* fHistEMCALM20                = new TH1F("fHistEMCALM20","M20 EMCAL;M20;Counts/0.1",200,0,20);
    TH1F* fHistEMCALTOF                = new TH1F("fHistEMCALTOF","TOF EMCAL;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistEMCALBeta               = new TH1F("fHistEMCALBeta","#beta EMCAL;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistEMCALE                  = new TH1F("fHistEMCALE","Energy EMCAL;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistEMCALNTracks            = new TH1F("fHistEMCALNTracks","N Tracks matched to Cluster EMCAL;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistEMCALNCells             = new TH1F("fHistEMCALNCells","N Cells of Cluster EMCAL;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistEMCALDx                 = new TH1F("fHistEMCALDx","Dx (Closest Track in #phi) EMCAL;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALDz                 = new TH1F("fHistEMCALDz","Dz (Closest Track in z (#eta)) EMCAL;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALDy                 = new TH1F("fHistEMCALDy","Dy (Closest Track in #eta#phi) EMCAL;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALDisp               = new TH1F("fHistEMCALDisp","Dispersion EMCAL;Dispersion;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCDisttoBC         = new TH1F("fHistEMCALMCDisttoBC","Distance to Bad Channel EMCAL MC;Distance;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCM02              = new TH1F("fHistEMCALMCM02","M02 EMCAL MC;M02;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCM20              = new TH1F("fHistEMCALMCM20","M20 EMCAL MC;M20;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCTOF              = new TH1F("fHistEMCALMCTOF","TOF EMCAL MC;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistEMCALMCBeta             = new TH1F("fHistEMCALMCBeta","#beta EMCAL MC;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistEMCALMCE                = new TH1F("fHistEMCALMCE","Energy EMCAL MC;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistEMCALMCNTracks          = new TH1F("fHistEMCALMCNTracks","N Tracks matched to Cluster EMCAL MC;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistEMCALMCNCells           = new TH1F("fHistEMCALMCNCells","N Cells of Cluster EMCAL MC;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistEMCALMCDx               = new TH1F("fHistEMCALMCDx","Dx (Closest Track in #phi) EMCAL MC;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCDz               = new TH1F("fHistEMCALMCDz","Dz (Closest Track in z (#eta)) EMCAL MC;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCDy               = new TH1F("fHistEMCALMCDy","Dy (Closest Track in #eta#phi) EMCAL MC;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCDisp             = new TH1F("fHistEMCALMCDisp","Dispersion EMCAL MC;Dispersion;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCSigDisttoBC      = new TH1F("fHistEMCALMCSigDisttoBC","Distance to Bad Channel EMCAL MC Sigma;Distance;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCSigM02           = new TH1F("fHistEMCALMCSigM02","M02 EMCAL MC Sigma;M02;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCSigM20           = new TH1F("fHistEMCALMCSigM20","M20 EMCAL MC Sigma;M20;Counts/0.1",200,0,20);
    TH1F* fHistEMCALMCSigTOF           = new TH1F("fHistEMCALMCSigTOF","TOF EMCAL MC Sigma;TOF (ns);Counts/(100 ps)",500,0,50);
    TH1F* fHistEMCALMCSigBeta          = new TH1F("fHistEMCALMCSigBeta","#beta EMCAL MC Sigma;#beta;Counts/(0.01)",150,0,1.5);
    TH1F* fHistEMCALMCSigE             = new TH1F("fHistEMCALMCSigE","Energy EMCAL MC Sigma;E (GeV);Counts/MeV",12000,0,12);
    TH1F* fHistEMCALMCSigNTracks       = new TH1F("fHistEMCALMCSigNTracks","N Tracks matched to Cluster EMCAL MC Sigma;N Tracks;Counts",11,-0.5,10.5);
    TH1F* fHistEMCALMCSigNCells        = new TH1F("fHistEMCALMCSigNCells","N Cells of Cluster EMCAL MC Sigma;N Cells;Counts",10,0.5,10.5);
    TH1F* fHistEMCALMCSigDx            = new TH1F("fHistEMCALMCSigDx","Dx (Closest Track in #phi) EMCAL MC Sigma;Dx (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCSigDz            = new TH1F("fHistEMCALMCSigDz","Dz (Closest Track in z (#eta)) EMCAL MC Sigma;Dz (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCSigDy            = new TH1F("fHistEMCALMCSigDy","Dy (Closest Track in #eta#phi) EMCAL MC Sigma;Dy (cm);Counts/mm",500,0,50);
    TH1F* fHistEMCALMCSigDisp          = new TH1F("fHistEMCALMCSigDisp","Dispersion EMCAL MC Sigma;Dispersion;Counts/0.1",200,0,20);

    //Calculated photon
    TH1F* fHistCalcStatistics          = new TH1F("fHistCalcStatistics","Calculated Photon Counter;Stage;Counts",11,0.5,11.5);
    TH1F* fHistCalcStatisticsMC        = new TH1F("fHistCalcStatisticsMC","Calculated Photon Counter MC;Stage;Counts",11,0.5,11.5);
    TH1F* fHistEqPar1Calc              = new TH1F("fHistEqPar1Calc","Equation Parameter 1;par1;Counts/0.01",400,-2,2);
    TH1F* fHistEqPar1CalcMC            = new TH1F("fHistEqPar1CalcMC","MC Equation Parameter 1;par1;Counts/0.01",400,-2,2);
    TH1F* fHistNIterCalc               = new TH1F("fHistNIterCalc","Number of Iteraterations to fix NaN;N_{Iterations};Counts/Iteration",fMaxIter+2,-0.5,fMaxIter+1.5);
    TH1F* fHistNIterCalcMC             = new TH1F("fHistNIterCalcMC","MC Number of Iteraterations to fix NaN;N_{Iterations};Counts/Iteration",fMaxIter+2,-0.5,fMaxIter+1.5);
    TH1F* fHistDeltaThetaCalc          = new TH1F("fHistDeltaThetaCalc",Form("#Delta#Theta to fix NaN;#Delta#Theta (rad);Counts/(%f rad)",fThetaRange/200),400,-fThetaRange,fThetaRange);
    TH1F* fHistDeltaThetaCalcMC        = new TH1F("fHistDeltaThetaCalcMC",Form("MC #Delta#Theta to fix NaN;#Delta#Theta (rad);Counts/(%f rad)",fThetaRange/200),400,-fThetaRange,fThetaRange);
    TH1F* fHistDeltaPhiCalc            = new TH1F("fHistDeltaPhiCalc",Form("#Delta#Phi to fix NaN;#Delta#Phi (rad);Counts/(%f rad)",fPhiRange/200),400,-fPhiRange,fPhiRange);
    TH1F* fHistDeltaPhiCalcMC          = new TH1F("fHistDeltaPhiCalcMC",Form("MC #Delta#Phi to fix NaN;#Delta#Phi (rad);Counts/(%f rad)",fPhiRange/200),400,-fPhiRange,fPhiRange);
    TH1F* fHistInvSigmaMassCalc        = new TH1F("fHistInvSigmaMassCalc","Invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(2 MeV/#it{c}^{2})",1000,1,3);
    TH1F* fHistInvSigmaMassCalcMC      = new TH1F("fHistInvSigmaMassCalcMC","MC Invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})",1000,1,2);
    TH1F* fHistSigmaDCAtoPVCalc        = new TH1F("fHistSigmaDCAtoPVCalc","#Sigma^{+}/#bar#Sigma^{-} Dist to PV KF;DCA (cm);Counts/(mm)",500,0,50);
    TH1F* fHistSigmaDCAtoPVCalcMC      = new TH1F("fHistSigmaDCAtoPVCalcMC","MC #Sigma^{+}/#bar#Sigma^{-} Dist to PV KF;DCA (cm);Counts/(mm)",500,0,50);
    TH1F* fHistProtPhotonDCACalc       = new TH1F("fHistProtPhotonDCACalc","#Sigma^{+}/#bar#Sigma^{-} Photon Proton DCA KF;DCA (cm);Counts/(100 #mum)",500,0,5);
    TH1F* fHistProtPhotonDCACalcMC     = new TH1F("fHistProtPhotonDCACalcMC","MC #Sigma^{+}/#bar#Sigma^{-} Photon Proton DCA KF;DCA (cm);Counts/(100 #mum)",500,0,5);
    TH1F* fHistSigmaYCalc              = new TH1F("fHistSigmaYCalc","#Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistSigmaYCalcMC            = new TH1F("fHistSigmaYCalcMC","MC #Sigma^{+}/#bar#Sigma^{-} Y;Y;Counts/(0.01)",300,-1.5,1.5);
    TH1F* fHistSigmaAntiPACalc         = new TH1F("fHistSigmaAntiPACalc","#Sigma^{+}/#bar#Sigma^{-} Anti-PA;PA (rad);Counts/(0.005)",700,0,3.5);
    TH1F* fHistSigmaAntiPACalcMC       = new TH1F("fHistSigmaAntiPACalcMC","MC #Sigma^{+}/#bar#Sigma^{-} Anti-PA;PA (rad);Counts/(0.005)",700,0,3.5);

    //Miscellaneous
    /*** EMPTY ***/

    //Alphanumeric Histogram Labels

    fHistCutBookKeeper->GetXaxis()->SetBinLabel(1,"fRemoveGenPileup");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(2,"fMaxVertexZ");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(3,"fEvPoolSize");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(4,"fEvTrackSize");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(5,"fCentralityBins");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(6,"fMinCentBin");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(7,"fMaxCentBin");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(8,"fZvtxBins");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(9,"fMinZBin");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(10,"fMaxZBin");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(11,"fEvPoolSize2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(12,"fEvTrackSize2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(13,"fCentralityBins2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(14,"fMinCentBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(15,"fMaxCentBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(16,"fZvtxBins2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(17,"fMinZBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(18,"fMaxZBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(19,"fRequireSigma");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(20,"fRequireProton");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(21,"fRequireSigmaCand");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(22,"fUseAbsZ");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(23,"fUseAbsZCorr");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(24,"fRejectNegIDs");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(25,"fRejectZeroFilterBit");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(26,"fMaxProtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(27,"fMinTPCClustProt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(28,"fMaxNsigProtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(29,"fRequireProtonTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(30,"fRequireProtonTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(31,"fRequireProtonTOFforPairs");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(32,"fMaxNsigProtTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(33,"fMaxpOnlyTPCPID");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(34,"fMinProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(35,"fMaxProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(36,"fStrictMaxProtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(37,"fStrictMinTPCClustProt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(38,"fStrictMaxNsigProtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(39,"fStrictMaxNsigProtTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(40,"fStrictMaxpOnlyTPCPID");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(41,"fStrictMinProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(42,"fStrictMaxProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(43,"fMaxMCEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(44,"fMaxDaughtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(45,"fMinTPCClustDaught");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(46,"fMaxNsigDaughtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(47,"fMaxalpha");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(48,"fMaxqt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(49,"fMaxopenangle");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(50,"fMaxdeltatheta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(51,"fMinV0CPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(52,"fMinV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(53,"fMaxV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(54,"fMaxphotonmass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(55,"fRequirePHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(56,"fMinClusterBeta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(57,"fMinClusterDy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(58,"fMaxClusterM02");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(59,"fCleanAutoCorr");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(60,"fMinPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(61,"fMaxPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(62,"fMaxSigmaPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(63,"fMaxSigmaY");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(64,"fMaxSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(65,"fMinProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(66,"fMinProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(67,"fMaxProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(68,"fMaxProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(69,"fRequireDCACut");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(70,"fMinPi0MassPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(71,"fMaxPi0MassPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(72,"fMaxSigmaPAPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(73,"fMaxSigmaPAPHOSHM");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(74,"fMinSigmaAntiPAPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(75,"fMaxProtPhotDCA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(76,"fMinSigmaDCAtoPVPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(77,"fMaxSigmaDCAtoPVPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(78,"fMaxSigmaYPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(79,"fMaxSigmaMassPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(80,"fMinProtonDCAxyPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(81,"fMinProtonDCAzPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(82,"fMaxProtonDCAxyPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(83,"fMaxProtonDCAzPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(84,"fRequireDCACutPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(85,"fMinCorrPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(86,"fMaxCorrPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(87,"fMaxCorrSigmaPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(88,"fMinCorrSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(89,"fMaxCorrSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(90,"fMinCorrProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(91,"fMaxCorrPairProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(92,"fMaxCorrPairProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(93,"fMaxCorrkstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(94,"fMinCorrPi0MassPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(95,"fMaxCorrPi0MassPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(96,"fMaxCorrSigmaPAPHOS");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(97,"fMaxNsigPionTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(98,"fMaxNsigKaonTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(99,"fMaxNsigElecTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(100,"fMaxNsigPionTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(101,"fMaxNsigKaonTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(102,"fMaxNsigElecTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(103,"fRejOtherTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(104,"fCheckProtonV0IDs");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(105,"fMaxSigmaMassCalc");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(106,"fMinSigmaAntiPACalc");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(107,"fMaxSigmaAntiPACalc");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(108,"fMaxIter");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(109,"fRejectNegPar1");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(110,"fMaxSigmaYCalc");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(111,"fThetaPar0");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(112,"fThetaPar1");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(113,"fThetaRange");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(114,"fPhiPar0");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(115,"fPhiPar1");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(116,"fPhiRange");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(117,"Number of Fills");

    fHistCutBookKeeper->SetBinContent(1,fRemoveGenPileup);
    fHistCutBookKeeper->SetBinContent(2,fMaxVertexZ);
    fHistCutBookKeeper->SetBinContent(3,fEvPoolSize);
    fHistCutBookKeeper->SetBinContent(4,fEvTrackSize);
    fHistCutBookKeeper->SetBinContent(5,fCentralityBins);
    fHistCutBookKeeper->SetBinContent(6,fMinCentBin);
    fHistCutBookKeeper->SetBinContent(7,fMaxCentBin);
    fHistCutBookKeeper->SetBinContent(8,fZvtxBins);
    fHistCutBookKeeper->SetBinContent(9,fMinZBin);
    fHistCutBookKeeper->SetBinContent(10,fMaxZBin);
    fHistCutBookKeeper->SetBinContent(11,fEvPoolSize2);
    fHistCutBookKeeper->SetBinContent(12,fEvTrackSize2);
    fHistCutBookKeeper->SetBinContent(13,fCentralityBins2);
    fHistCutBookKeeper->SetBinContent(14,fMinCentBin2);
    fHistCutBookKeeper->SetBinContent(15,fMaxCentBin2);
    fHistCutBookKeeper->SetBinContent(16,fZvtxBins2);
    fHistCutBookKeeper->SetBinContent(17,fMinZBin2);
    fHistCutBookKeeper->SetBinContent(18,fMaxZBin2);
    fHistCutBookKeeper->SetBinContent(19,fRequireSigma);
    fHistCutBookKeeper->SetBinContent(20,fRequireProton);
    fHistCutBookKeeper->SetBinContent(21,fRequireSigmaCand);
    fHistCutBookKeeper->SetBinContent(22,fUseAbsZ);
    fHistCutBookKeeper->SetBinContent(23,fUseAbsZCorr);
    fHistCutBookKeeper->SetBinContent(24,fRejectNegIDs);
    fHistCutBookKeeper->SetBinContent(25,fRejectZeroFilterBit);
    fHistCutBookKeeper->SetBinContent(26,fMaxProtEta);
    fHistCutBookKeeper->SetBinContent(27,fMinTPCClustProt);
    fHistCutBookKeeper->SetBinContent(28,fMaxNsigProtTPC);
    fHistCutBookKeeper->SetBinContent(29,fRequireProtonTPC);
    fHistCutBookKeeper->SetBinContent(30,fRequireProtonTOF);
    fHistCutBookKeeper->SetBinContent(31,fRequireProtonTOFforPairs);
    fHistCutBookKeeper->SetBinContent(32,fMaxNsigProtTOF);
    fHistCutBookKeeper->SetBinContent(33,fMaxpOnlyTPCPID);
    fHistCutBookKeeper->SetBinContent(34,fMinProtpt);
    fHistCutBookKeeper->SetBinContent(35,fMaxProtpt);
    fHistCutBookKeeper->SetBinContent(36,fStrictMaxProtEta);
    fHistCutBookKeeper->SetBinContent(37,fStrictMinTPCClustProt);
    fHistCutBookKeeper->SetBinContent(38,fStrictMaxNsigProtTPC);
    fHistCutBookKeeper->SetBinContent(39,fStrictMaxNsigProtTOF);
    fHistCutBookKeeper->SetBinContent(40,fStrictMaxpOnlyTPCPID);
    fHistCutBookKeeper->SetBinContent(41,fStrictMinProtpt);
    fHistCutBookKeeper->SetBinContent(42,fStrictMaxProtpt);
    fHistCutBookKeeper->SetBinContent(43,fMaxMCEta);
    fHistCutBookKeeper->SetBinContent(44,fMaxDaughtEta);
    fHistCutBookKeeper->SetBinContent(45,fMinTPCClustDaught);
    fHistCutBookKeeper->SetBinContent(46,fMaxNsigDaughtTPC);
    fHistCutBookKeeper->SetBinContent(47,fMaxalpha);
    fHistCutBookKeeper->SetBinContent(48,fMaxqt);
    fHistCutBookKeeper->SetBinContent(49,fMaxopenangle);
    fHistCutBookKeeper->SetBinContent(50,fMaxdeltatheta);
    fHistCutBookKeeper->SetBinContent(51,fMinV0CPA);
    fHistCutBookKeeper->SetBinContent(52,fMinV0Radius);
    fHistCutBookKeeper->SetBinContent(53,fMaxV0Radius);
    fHistCutBookKeeper->SetBinContent(54,fMaxphotonmass);
    fHistCutBookKeeper->SetBinContent(55,fRequirePHOS);
    fHistCutBookKeeper->SetBinContent(56,fMinClusterBeta);
    fHistCutBookKeeper->SetBinContent(57,fMinClusterDy);
    fHistCutBookKeeper->SetBinContent(58,fMaxClusterM02);
    fHistCutBookKeeper->SetBinContent(59,fCleanAutoCorr); 
    fHistCutBookKeeper->SetBinContent(60,fMinPi0Mass);
    fHistCutBookKeeper->SetBinContent(61,fMaxPi0Mass);
    fHistCutBookKeeper->SetBinContent(62,fMaxSigmaPA);
    fHistCutBookKeeper->SetBinContent(63,fMaxSigmaY);
    fHistCutBookKeeper->SetBinContent(64,fMaxSigmaMass);
    fHistCutBookKeeper->SetBinContent(65,fMinProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(66,fMinProtonDCAz);
    fHistCutBookKeeper->SetBinContent(67,fMaxProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(68,fMaxProtonDCAz);
    fHistCutBookKeeper->SetBinContent(69,fRequireDCACut);
    fHistCutBookKeeper->SetBinContent(70,fMinPi0MassPHOS);
    fHistCutBookKeeper->SetBinContent(71,fMaxPi0MassPHOS);
    fHistCutBookKeeper->SetBinContent(72,fMaxSigmaPAPHOS);
    fHistCutBookKeeper->SetBinContent(73,fMaxSigmaPAPHOSHM);
    fHistCutBookKeeper->SetBinContent(74,fMinSigmaAntiPAPHOS);
    fHistCutBookKeeper->SetBinContent(75,fMaxProtPhotDCA);
    fHistCutBookKeeper->SetBinContent(76,fMinSigmaDCAtoPVPHOS);
    fHistCutBookKeeper->SetBinContent(77,fMaxSigmaDCAtoPVPHOS);
    fHistCutBookKeeper->SetBinContent(78,fMaxSigmaYPHOS);
    fHistCutBookKeeper->SetBinContent(79,fMaxSigmaMassPHOS);
    fHistCutBookKeeper->SetBinContent(80,fMinProtonDCAxyPHOS);
    fHistCutBookKeeper->SetBinContent(81,fMinProtonDCAzPHOS);
    fHistCutBookKeeper->SetBinContent(82,fMaxProtonDCAxyPHOS);
    fHistCutBookKeeper->SetBinContent(83,fMaxProtonDCAzPHOS);
    fHistCutBookKeeper->SetBinContent(84,fRequireDCACutPHOS);
    fHistCutBookKeeper->SetBinContent(85,fMinCorrPi0Mass);
    fHistCutBookKeeper->SetBinContent(86,fMaxCorrPi0Mass);
    fHistCutBookKeeper->SetBinContent(87,fMaxCorrSigmaPA);
    fHistCutBookKeeper->SetBinContent(88,fMinCorrSigmaMass);
    fHistCutBookKeeper->SetBinContent(89,fMaxCorrSigmaMass);
    fHistCutBookKeeper->SetBinContent(90,fMinCorrProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(91,fMaxCorrPairProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(92,fMaxCorrPairProtonDCAz);
    fHistCutBookKeeper->SetBinContent(93,fMaxCorrkstar);
    fHistCutBookKeeper->SetBinContent(94,fMinCorrPi0MassPHOS);
    fHistCutBookKeeper->SetBinContent(95,fMaxCorrPi0MassPHOS);
    fHistCutBookKeeper->SetBinContent(96,fMaxCorrSigmaPAPHOS);
    fHistCutBookKeeper->SetBinContent(97,fMaxNsigPionTPC);
    fHistCutBookKeeper->SetBinContent(98,fMaxNsigKaonTPC);
    fHistCutBookKeeper->SetBinContent(99,fMaxNsigElecTPC);
    fHistCutBookKeeper->SetBinContent(100,fMaxNsigPionTOF);
    fHistCutBookKeeper->SetBinContent(101,fMaxNsigKaonTOF);
    fHistCutBookKeeper->SetBinContent(102,fMaxNsigElecTOF);
    fHistCutBookKeeper->SetBinContent(103,fRejOtherTOF);
    fHistCutBookKeeper->SetBinContent(104,fCheckProtonV0IDs);
    fHistCutBookKeeper->SetBinContent(105,fMaxSigmaMassCalc);
    fHistCutBookKeeper->SetBinContent(106,fMinSigmaAntiPACalc);
    fHistCutBookKeeper->SetBinContent(107,fMaxSigmaAntiPACalc);
    fHistCutBookKeeper->SetBinContent(108,fMaxIter);
    fHistCutBookKeeper->SetBinContent(109,fRejectNegPar1);
    fHistCutBookKeeper->SetBinContent(110,fMaxSigmaYCalc);
    fHistCutBookKeeper->SetBinContent(111,fThetaPar0);
    fHistCutBookKeeper->SetBinContent(112,fThetaPar1);
    fHistCutBookKeeper->SetBinContent(113,fThetaRange);
    fHistCutBookKeeper->SetBinContent(114,fPhiPar0);
    fHistCutBookKeeper->SetBinContent(115,fPhiPar1);
    fHistCutBookKeeper->SetBinContent(116,fPhiRange);
    fHistCutBookKeeper->SetBinContent(117,1);

    fHistMCCounter->GetXaxis()->SetBinLabel(1,"Events");
    fHistMCCounter->GetXaxis()->SetBinLabel(2,"MC Particles");
    fHistMCCounter->GetXaxis()->SetBinLabel(3,"p");
    fHistMCCounter->GetXaxis()->SetBinLabel(4,"#barp");
    fHistMCCounter->GetXaxis()->SetBinLabel(5,"#Delta^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(6,"#bar#Delta^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(7,"#Sigma^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(8,"#bar#Sigma^{-}");    
    fHistMCCounter->GetXaxis()->SetBinLabel(9,"pi^{0}");
    fHistMCCounter->GetXaxis()->SetBinLabel(10,"pi^{0} from #Sigma^{+}/#bar#Sigma^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(11,"#gamma");
    fHistMCCounter->GetXaxis()->SetBinLabel(12,"#gamma->e^{-}+e^{+} (R_{conv}<180cm)"); 
    
    fHistEventCounter->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdouble->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterHM->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdoubleHM->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterHMV0->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdoubleHMV0->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterHMSPD->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdoubleHMSPD->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterINT7->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdoubleINT7->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterINT7LF->GetXaxis()->SetBinLabel(1,"Events");
    fHistEventCounterdoubleINT7LF->GetXaxis()->SetBinLabel(1,"Events");

    fHistEventCounter->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdouble->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterHM->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdoubleHM->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterHMV0->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdoubleHMV0->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterHMSPD->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdoubleHMSPD->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterINT7->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdoubleINT7->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterINT7LF->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");
    fHistEventCounterdoubleINT7LF->GetXaxis()->SetBinLabel(2,"Events, |z_{Vertex}|<10cm");

    fHistV0Statistics->GetXaxis()->SetBinLabel(1, "On-the-fly V0s");
    fHistV0Statistics->GetXaxis()->SetBinLabel(2, "Offline V0s");
    fHistV0Statistics->GetXaxis()->SetBinLabel(3, "Offline only V0s");
    fHistV0Statistics->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistV0Statistics->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistV0Statistics->GetXaxis()->SetBinLabel(6, "Passed alpha cut");
    fHistV0Statistics->GetXaxis()->SetBinLabel(7, "Passed qt cut");
    fHistV0Statistics->GetXaxis()->SetBinLabel(8, "Passed open angle");
    fHistV0Statistics->GetXaxis()->SetBinLabel(9, "Passed delta theta");
    fHistV0Statistics->GetXaxis()->SetBinLabel(10,"Passed CosPA cut");
    fHistV0Statistics->GetXaxis()->SetBinLabel(11,"Passed TPC PID");
    fHistV0Statistics->GetXaxis()->SetBinLabel(12,"Passed min radius");
    fHistV0Statistics->GetXaxis()->SetBinLabel(13,"Passed max radius");
    fHistV0Statistics->GetXaxis()->SetBinLabel(14,"Passed inv mass");

    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(1, "On-the-fly V0s");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(2, "Offline V0s");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(3, "Offline only V0s");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(6, "Passed alpha cut");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(7, "Passed qt cut");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(8, "Passed open angle");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(9, "Passed delta theta");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(10,"Passed CosPA cut");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(11,"Passed TPC PID");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(12,"Passed min radius");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(13,"Passed max radius");
    fHistV0StatisticsMC->GetXaxis()->SetBinLabel(14,"Passed inv mass");

    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(1, "On-the-fly V0s");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(2, "Offline V0s");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(3, "Offline only V0s");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(6, "Passed alpha cut");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(7, "Passed qt cut");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(8, "Passed open angle");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(9, "Passed delta theta");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(10,"Passed CosPA cut");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(11,"Passed TPC PID");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(12,"Passed min radius");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(13,"Passed max radius");
    fHistV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(14,"Passed inv mass");

    fHistProtonStatistics->GetXaxis()->SetBinLabel(1, "Tracks");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(2, "ESD ID >= 0");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(3, "Filterbit != 0");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(6, "Passed TPC cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(7, "Passed TOF cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(8, "Passed pt cut");

    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(1, "Tracks");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(2, "ESD ID >= 0");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(3, "Filterbit != 0");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(6, "Passed TPC cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(7, "Passed TOF cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(8, "Passed pt cut");

    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(1, "Tracks");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(2, "ESD ID >= 0");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(3, "Filterbit != 0");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(4, "Passed eta cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(5, "Passed cluster cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(6, "Passed TPC cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(7, "Passed TOF cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(8, "Passed pt cut");

    fHistAddV0Statistics->GetXaxis()->SetBinLabel(1, "Pairs considered");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(2, "Pairs not Onfly V0s");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(3, "Passed daughter dca");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(4, "Vertex within volume");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(5, "Passed radius cut");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(6, "Passed CosPA cut");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(7, "Passed alpha cut");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(8, "Passed qt cut");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(9, "Passed open angle");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(10,"Passed delta theta");
    fHistAddV0Statistics->GetXaxis()->SetBinLabel(11,"Passed inv mass");

    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(1, "Pairs considered");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(2, "Pairs not Onfly V0s");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(3, "Passed daughter dca");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(4, "Vertex within volume");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(5, "Passed radius cut");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(6, "Passed CosPA cut");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(7, "Passed alpha cut");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(8, "Passed qt cut");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(9, "Passed open angle");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(10,"Passed delta theta");
    fHistAddV0StatisticsMC->GetXaxis()->SetBinLabel(11,"Passed inv mass");

    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(1, "Pairs considered");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(2, "Pairs not Onfly V0s");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(3, "Passed daughter dca");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(4, "Vertex within volume");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(5, "Passed radius cut");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(6, "Passed CosPA cut");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(7, "Passed alpha cut");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(8, "Passed qt cut");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(9, "Passed open angle");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(10,"Passed delta theta");
    fHistAddV0StatisticsSigmaMC->GetXaxis()->SetBinLabel(11,"Passed inv mass");

    fHistSigmaCounter->GetXaxis()->SetBinLabel(1,"Both from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(2,"One from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(3,"Additional #gammas");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(4,"4 Particle Reco");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(5,"#Sigma->p#gamma from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(6,"#Sigma->p#gamma additional");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(6,"PHOS");

    fHistGammaPairStats->GetXaxis()->SetBinLabel(1,"Pairs");
    fHistGammaPairStats->GetXaxis()->SetBinLabel(2,"#pi^{0}");
    fHistGammaPairStats->GetXaxis()->SetBinLabel(3,"Same #gamma");
    fHistGammaPairStats->GetXaxis()->SetBinLabel(4,"Pairs, diff. Idx");
    fHistGammaPairStats->GetXaxis()->SetBinLabel(5,"#pi^{0}, diff. Idx");
    fHistGammaPairStats->GetXaxis()->SetBinLabel(6,"Same #gamma, diff. Idx");

    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(1,"Pairs");
    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(2,"#pi^{0}");
    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(3,"Same #gamma");
    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(4,"Pairs, diff. Idx");
    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(5,"#pi^{0}, diff. Idx");
    fHistGammaPairStatsOneadd->GetXaxis()->SetBinLabel(6,"Same #gamma, diff. Idx");

    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(1,"Pairs");
    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(2,"#pi^{0}");
    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(3,"Same #gamma");
    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(4,"Pairs, diff. Idx");
    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(5,"#pi^{0}, diff. Idx");
    fHistGammaPairStatsOnlyadd->GetXaxis()->SetBinLabel(6,"Same #gamma, diff. Idx");

    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(1,"Total");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(2,"Primary");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(5,"K^{0}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(6,"#barK^{0}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(7,"K^{-}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(8,"K^{+}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(9,"pi^{-}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(10,"pi^{+}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(11,"n");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(12,"#barn");    
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(13,"p");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(14,"#barp");    
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(15,"#Lambda");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(16,"#bar#Lambda");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(17,"#Sigma^{+}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(18,"#Sigma^{-}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(19,"#bar#Sigma^{-}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(20,"#bar#Sigma^{+}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(21,"#Xi^{0}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(22,"#bar#Xi^{0}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(23,"#Xi^{-}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(24,"#bar#Xi^{+}");
    fHistSigmaMotherPart->GetXaxis()->SetBinLabel(25,"Other");

    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(1,"Total");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(2,"Primary");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(3,"K^{0}_{L}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(4,"K^{0}_{S}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(5,"K^{0}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(6,"#barK^{0}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(7,"K^{-}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(8,"K^{+}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(9,"pi^{-}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(10,"pi^{+}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(11,"n");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(12,"#barn");    
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(13,"p");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(14,"#barp");    
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(15,"#Lambda");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(16,"#bar#Lambda");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(17,"#Sigma^{+}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(18,"#Sigma^{-}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(19,"#bar#Sigma^{-}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(20,"#bar#Sigma^{+}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(21,"#Xi^{0}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(22,"#bar#Xi^{0}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(23,"#Xi^{-}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(24,"#bar#Xi^{+}");
    fHistAntiSigmaMotherPart->GetXaxis()->SetBinLabel(25,"Other");

    fHistClusterStatistics->GetXaxis()->SetBinLabel(1, "Calo Clusters");
    fHistClusterStatistics->GetXaxis()->SetBinLabel(2, "Passed PHOS-only cut");
    fHistClusterStatistics->GetXaxis()->SetBinLabel(3, "Passed Beta cut");
    fHistClusterStatistics->GetXaxis()->SetBinLabel(4, "Passed Dy cut");
    fHistClusterStatistics->GetXaxis()->SetBinLabel(5, "Passed M02 cut");

    fHistClusterStatisticsMC->GetXaxis()->SetBinLabel(1, "Calo Clusters");
    fHistClusterStatisticsMC->GetXaxis()->SetBinLabel(2, "Passed PHOS-only cut");
    fHistClusterStatisticsMC->GetXaxis()->SetBinLabel(3, "Passed Beta cut");
    fHistClusterStatisticsMC->GetXaxis()->SetBinLabel(4, "Passed Dy cut");
    fHistClusterStatisticsMC->GetXaxis()->SetBinLabel(5, "Passed M02 cut");

    fHistClusterStatisticsMCSig->GetXaxis()->SetBinLabel(1, "Calo Clusters");
    fHistClusterStatisticsMCSig->GetXaxis()->SetBinLabel(2, "Passed PHOS-only cut");
    fHistClusterStatisticsMCSig->GetXaxis()->SetBinLabel(3, "Passed Beta cut");
    fHistClusterStatisticsMCSig->GetXaxis()->SetBinLabel(4, "Passed Dy cut");
    fHistClusterStatisticsMCSig->GetXaxis()->SetBinLabel(5, "Passed M02 cut");

    fHistCalcStatistics->GetXaxis()->SetBinLabel(1, "Proton Photon Pairs");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(2, "ID check ok");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(3, "Passed DCAtoPV Cut");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(4, "KF test ok");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(5, "Passed FlightDist Cut");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(6,"Passed ProtPhotDCA Cut");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(7,"Passed Anti-PA Cut");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(8, "Par2 not NaN");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(9, "Phase space ok");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(10, "Passed mass cut");
    fHistCalcStatistics->GetXaxis()->SetBinLabel(11,"Passed Rapidity Cut");

    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(1, "Proton Photon Pairs");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(2, "ID check ok");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(3, "Passed DCAtoPV Cut");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(4, "KF test ok");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(5, "Passed FlightDist Cut");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(6,"Passed ProtPhotDCA Cut");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(7,"Passed Anti-PA Cut");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(8, "Par2 not NaN");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(9, "Phase space ok");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(10, "Passed mass cut");
    fHistCalcStatisticsMC->GetXaxis()->SetBinLabel(11,"Passed Rapidity Cut");

    //Add Histograms to Output List
    
    //Book Keeper for used Cuts
    fOutputList->Add(fHistCutBookKeeper);
    //Event related
    fOutputList->Add(fHistMCGenPileup);
    fOutputList->Add(fHistVertexZ);
    fOutputList->Add(fHistVertexZMC);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistCentralityHMV0);
    fOutputList->Add(fHistCentralityHMSPD);
    fOutputList->Add(fHistCentralityINT7);
    fOutputList->Add(fHistCentralityINT7LF);

    fOutputList->Add(fHistEventCounter);
    fOutputList->Add(fHistEventCounterdouble);
    fOutputList->Add(fHistEventCounterHM);
    fOutputList->Add(fHistEventCounterdoubleHM);
    fOutputList->Add(fHistEventCounterHMV0);
    fOutputList->Add(fHistEventCounterdoubleHMV0);
    fOutputList->Add(fHistEventCounterHMSPD);
    fOutputList->Add(fHistEventCounterdoubleHMSPD);
    fOutputList->Add(fHistEventCounterINT7);
    fOutputList->Add(fHistEventCounterdoubleINT7);
    fOutputList->Add(fHistEventCounterINT7LF);
    fOutputList->Add(fHistEventCounterdoubleINT7LF);

    fOutputList->Add(fHistTrackMultiplicity);
    fOutputList->Add(fHistRefMultComb05);
    fOutputList->Add(fHistRefMultComb08);
    fOutputList->Add(fHistRefMultComb10);
    fOutputList->Add(fHistRefMultComb05HM);
    fOutputList->Add(fHistRefMultComb08HM);
    fOutputList->Add(fHistRefMultComb10HM);
    fOutputList->Add(fHistRefMultComb05HMV0);
    fOutputList->Add(fHistRefMultComb08HMV0);
    fOutputList->Add(fHistRefMultComb10HMV0);
    fOutputList->Add(fHistRefMultComb05HMSPD);
    fOutputList->Add(fHistRefMultComb08HMSPD);
    fOutputList->Add(fHistRefMultComb10HMSPD);
    fOutputList->Add(fHistRefMultComb05INT7);
    fOutputList->Add(fHistRefMultComb08INT7);
    fOutputList->Add(fHistRefMultComb10INT7);
    fOutputList->Add(fHistRefMultComb05INT7LF);
    fOutputList->Add(fHistRefMultComb08INT7LF);
    fOutputList->Add(fHistRefMultComb10INT7LF);

    //Counters
    fOutputList->Add(fHistMCCounter);
    fOutputList->Add(fHistV0Statistics);
    fOutputList->Add(fHistV0StatisticsMC);
    fOutputList->Add(fHistV0StatisticsSigmaMC);
    fOutputList->Add(fHistProtonStatistics);
    fOutputList->Add(fHistProtonStatisticsMC);
    fOutputList->Add(fHistProtonStatisticsSigmaMC);
    fOutputList->Add(fHistAddV0Statistics);
    fOutputList->Add(fHistAddV0StatisticsMC);
    fOutputList->Add(fHistAddV0StatisticsSigmaMC);
    fOutputList->Add(fHistGammaPairStats);    
    fOutputList->Add(fHistGammaPairStatsOneadd);
    fOutputList->Add(fHistGammaPairStatsOnlyadd);
    fOutputList->Add(fHistSigmaCounter);
    //MC Kinematics
    fOutputList->Add(fHistMCPrimSigmaPtvsRap);
    fOutputList->Add(fHistMCPrimAntiSigmaPtvsRap);
    fOutputList->Add(fHistMCPrimSigmaPtRap05);
    fOutputList->Add(fHistMCPrimAntiSigmaPtRap05);
    fOutputList->Add(fHistMCPrimSigmaPtRap08);
    fOutputList->Add(fHistMCPrimAntiSigmaPtRap08);
    fOutputList->Add(fHistMCPrimSigmaPtRap10);
    fOutputList->Add(fHistMCPrimAntiSigmaPtRap10);
    fOutputList->Add(fHistMCSigmaPtvsRap);
    fOutputList->Add(fHistMCAntiSigmaPtvsRap);
    fOutputList->Add(fHistMCSigmaPtRap05);
    fOutputList->Add(fHistMCAntiSigmaPtRap05);
    fOutputList->Add(fHistMCSigmaPtRap08);
    fOutputList->Add(fHistMCAntiSigmaPtRap08);
    fOutputList->Add(fHistMCSigmaPtRap10);
    fOutputList->Add(fHistMCAntiSigmaPtRap10);                  
    fOutputList->Add(fHistMCSigmaPt);
    fOutputList->Add(fHistMCAntiSigmaPt);
    fOutputList->Add(fHistMCPrimSigmaPt);
    fOutputList->Add(fHistMCPrimAntiSigmaPt);
    fOutputList->Add(fHistMCPrimSigmaPt08);
    fOutputList->Add(fHistMCPrimAntiSigmaPt08);
    fOutputList->Add(fHistMCPrimSigmaPt10);
    fOutputList->Add(fHistMCPrimAntiSigmaPt10);
    fOutputList->Add(fHistMCPrimSigmaPtvsEta);
    fOutputList->Add(fHistMCPrimAntiSigmaPtvsEta);
    fOutputList->Add(fHistMCSigmaPtvsEta);
    fOutputList->Add(fHistMCAntiSigmaPtvsEta);
    fOutputList->Add(fHistMCSigmaOrigin);
    fOutputList->Add(fHistMCAntiSigmaOrigin);
    fOutputList->Add(fHistMCDeltaPt);    
    fOutputList->Add(fHistMCAntiDeltaPt);
    fOutputList->Add(fHistMCPi0Pt);      
    fOutputList->Add(fHistMCProtonPt);         
    fOutputList->Add(fHistMCAntiProtonPt);
    fOutputList->Add(fHistMCPrimProtonPt);
    fOutputList->Add(fHistMCPrimAntiProtonPt);     
    fOutputList->Add(fHistMCPhotonPt);         
    fOutputList->Add(fHistMCSigmaPhotonPt);
    fOutputList->Add(fHistMCConvPhotonPt);     
    fOutputList->Add(fHistMCSigmaConvPhotonPt);
    fOutputList->Add(fHistMCConvRadius);    
    fOutputList->Add(fHistMCConvRadiusvspt);    
    fOutputList->Add(fHistSigmaMotherPart);    
    fOutputList->Add(fHistAntiSigmaMotherPart);    
    //Track Quality
    fOutputList->Add(fHistTrackEtaPhi);
    fOutputList->Add(fHistTrackChi2);
    fOutputList->Add(fHistTrackTPCCluster);
    fOutputList->Add(fHistTrackITSCluster);
    fOutputList->Add(fHistTrackpvsdEdx);
    fOutputList->Add(fHistTrackpvsbeta);
    fOutputList->Add(fHistTrackpt);
    //Proton QA
    fOutputList->Add(fHistProtonEtaPhiMC);
    fOutputList->Add(fHistProtonChi2MC);
    fOutputList->Add(fHistProtonTPCClusterMC);
    fOutputList->Add(fHistProtonITSClusterMC);
    fOutputList->Add(fHistProtonpvsNSigmaTPC);
    fOutputList->Add(fHistProtonpvsNSigmaTPCMC);
    fOutputList->Add(fHistProtonpvsNSigmaTOF);
    fOutputList->Add(fHistProtonpvsNSigmaTOFMC);
    fOutputList->Add(fHistProtonDCAxy);
    fOutputList->Add(fHistProtonDCAz);
    fOutputList->Add(fHistProtonDCAxyMC);
    fOutputList->Add(fHistProtonDCAzMC);
    fOutputList->Add(fHistProtonDCAxyMCSigma);
    fOutputList->Add(fHistProtonDCAzMCSigma);
    fOutputList->Add(fHistProtonDCAxyMCPrimSig);
    fOutputList->Add(fHistProtonDCAzMCPrimSig);
    fOutputList->Add(fHistPrimProtonDCAxyMC);
    fOutputList->Add(fHistPrimProtonDCAzMC);
    fOutputList->Add(fHistMaterialProtonDCAxyMC);
    fOutputList->Add(fHistMaterialProtonDCAzMC);
    fOutputList->Add(fHistWeakProtonDCAxyMC);
    fOutputList->Add(fHistWeakProtonDCAzMC);
    fOutputList->Add(fHistLambdaProtonDCAxyMC);
    fOutputList->Add(fHistLambdaProtonDCAzMC);
    fOutputList->Add(fHistProtonptMC);
    fOutputList->Add(fHistProtonptwCutsMC);
    fOutputList->Add(fHistProtonptwCuts);
    //V0 QA
    fOutputList->Add(fHistV0OnflyvsOffline);
    fOutputList->Add(fHistV0DaughtEtaPhi);
    fOutputList->Add(fHistV0DaughtEtaPhiMC);
    fOutputList->Add(fHistV0DaughtChi2);
    fOutputList->Add(fHistV0DaughtChi2MC);
    fOutputList->Add(fHistV0DaughtTPCClust);
    fOutputList->Add(fHistV0DaughtITSClust);
    fOutputList->Add(fHistV0DaughtTPCClustMC);
    fOutputList->Add(fHistV0DaughtITSClustMC);
    fOutputList->Add(fHistV0DaughtpvsNSigmaTPC);
    fOutputList->Add(fHistV0DaughtDCAtoPV);
    fOutputList->Add(fHistV0DaughtDCAtoPVMC);
    fOutputList->Add(fHistV0DaughtDCA);
    fOutputList->Add(fHistV0DaughtDCAMC);
    fOutputList->Add(fHistV0CPA);
    fOutputList->Add(fHistV0CPAMC);
    fOutputList->Add(fHistV0CPAMCSigma);
    fOutputList->Add(fHistV0Radius);
    fOutputList->Add(fHistV0RadiusMC);
    fOutputList->Add(fHistV0Position2D);
    fOutputList->Add(fHistV0Position2DMC);
    fOutputList->Add(fHistV0PhotonDCAPV);    
    fOutputList->Add(fHistV0PhotonDCAPVMC);
    fOutputList->Add(fHistV0PhotonDCAPVMCSigma);
    fOutputList->Add(fHistV0ArmPod);
    fOutputList->Add(fHistV0ArmPodMC);
    fOutputList->Add(fHistV0OpenAngle);
    fOutputList->Add(fHistV0OpenAngleMC);
    fOutputList->Add(fHistV0DeltaTheta);
    fOutputList->Add(fHistV0DeltaThetaMC);
    fOutputList->Add(fHistV0InvMass);
    fOutputList->Add(fHistV0InvMassMC);
    fOutputList->Add(fHistV0ptMC);
    fOutputList->Add(fHistV0ptwCutsMC);
    fOutputList->Add(fHistV0SigmaptMC); 
    fOutputList->Add(fHistV0SigmaptwCutsMC); 
    fOutputList->Add(fHistV0ptwCuts);
    //Gamma Gamma QA
    fOutputList->Add(fHistGammaPairInvMass);
    fOutputList->Add(fHistGammaPairInvMassMC);
    fOutputList->Add(fHistGammaPairInvMassOnfly);
    fOutputList->Add(fHistGammaPairInvMassMCOnfly);
    fOutputList->Add(fHistGammaPairInvMassOneAdd);  
    fOutputList->Add(fHistGammaPairInvMassOneAddMC);
    fOutputList->Add(fHistGammaPairInvMassOnlyAdd); 
    fOutputList->Add(fHistGammaPairInvMassOnlyAddMC);
    fOutputList->Add(fHistGammaPairDCA);
    fOutputList->Add(fHistGammaPairDCAMC);
    //Gamma Gamma QA - check auto correlations                      
    fOutputList->Add(fHistGammaPairInvMass2);
    fOutputList->Add(fHistGammaPairInvMassOnfly2);
    fOutputList->Add(fHistGammaPairInvMassOneAdd2);  
    fOutputList->Add(fHistGammaPairInvMassOnlyAdd2); 
    fOutputList->Add(fHistGammaPairDCA2);
    //Gamma Gamma QA - likely auto correlations                      
    fOutputList->Add(fHistGammaPairInvMass3);       
    fOutputList->Add(fHistGammaPairInvMassOnfly3);  
    fOutputList->Add(fHistGammaPairInvMassOneAdd3); 
    fOutputList->Add(fHistGammaPairInvMassOnlyAdd3);
    fOutputList->Add(fHistGammaPairDCA3);           
    //MC Sigma Topology
    fOutputList->Add(fHistKFSigmaVertexResX);
    fOutputList->Add(fHistKFSigmaVertexResY);
    fOutputList->Add(fHistKFSigmaVertexResZ);
    fOutputList->Add(fHistPi0VertexvsMC);     
    fOutputList->Add(fHistMCSigmaPA);
    fOutputList->Add(fHistMCSigmaPArot);
    fOutputList->Add(fHistMCPrimSigmaPA);
    fOutputList->Add(fHistMCPrimSigmaPArot);
    fOutputList->Add(fHistSigmaPA);
    fOutputList->Add(fHistSigmaPAmix);
    fOutputList->Add(fHistSigmaPArot);
    fOutputList->Add(fHistSigmaY);
    fOutputList->Add(fHistMCSigmaY);
    fOutputList->Add(fHistSigmaYrot);
    fOutputList->Add(fHistMCSigmaYrot);
    fOutputList->Add(fHistSigmaYmix);
    fOutputList->Add(fHistInvSigmaMass);
    fOutputList->Add(fHistMCInvSigmaMass);   
    fOutputList->Add(fHistMCInvSigmaMassrot);         
    fOutputList->Add(fHistInvSigmaMassmix);         
    fOutputList->Add(fHistMCOneGammaSigmaPA);
    fOutputList->Add(fHistMCPrimOneGammaSigmaPA);
    fOutputList->Add(fHistOneGammaSigmaPA);
    fOutputList->Add(fHistPi0VertexMC);
    //Sigma Momentum Resolution 
    fOutputList->Add(fHistSigmaPxResnoprop);
    fOutputList->Add(fHistSigmaPyResnoprop);
    fOutputList->Add(fHistSigmaPzResnoprop);
    fOutputList->Add(fHistSigmaPxResprop);
    fOutputList->Add(fHistSigmaPyResprop);
    fOutputList->Add(fHistSigmaPzResprop);
    //Event Mixing
    fOutputList->Add(fHistPairNMixedEvents);
    fOutputList->Add(fHistBkgNMixedEvents);
    fOutputList->Add(fHistPairNMixedEventsPHOS);
    fOutputList->Add(fHistBkgNMixedEventsPHOS);
    //KF Checks
    fOutputList->Add(fHistPhotonKFCheck);
    fOutputList->Add(fHistPi0KFCheck);
    fOutputList->Add(fHistSigmaKFCheck);
    //PHOS
    fOutputList->Add(fHistPhotonKFCheckPHOS);
    fOutputList->Add(fHistSigmaKFCheckPHOS);
    fOutputList->Add(fHistGammaPairInvMassPHOS);
    fOutputList->Add(fHistGammaPairInvMassMCPHOS);
    fOutputList->Add(fHistGammaPairInvMassOnflyPHOS);
    fOutputList->Add(fHistGammaPairInvMassMCOnflyPHOS);
    fOutputList->Add(fHistKFSigmaVertexResXPHOS);
    fOutputList->Add(fHistKFSigmaVertexResYPHOS);
    fOutputList->Add(fHistKFSigmaVertexResZPHOS);
    fOutputList->Add(fHistInvSigmaMassPHOS);
    fOutputList->Add(fHistMCInvSigmaMassPHOS);
    fOutputList->Add(fHistMCInvSigmaMassPHOSrot);
    fOutputList->Add(fHistInvSigmaMassPHOSUncorr);
    fOutputList->Add(fHistMCInvSigmaMassPHOSUncorr);
    fOutputList->Add(fHistSigmaYPHOS);
    fOutputList->Add(fHistMCSigmaYPHOS);
    fOutputList->Add(fHistMCSigmaPAPHOS);
    fOutputList->Add(fHistMCSigmaAntiPAPHOS);
    fOutputList->Add(fHistMCPrimSigmaPAPHOS);
    fOutputList->Add(fHistMCPrimSigmaPAPHOSrot);
    fOutputList->Add(fHistMCPrimSigmaAntiPAPHOS);
    fOutputList->Add(fHistSigmaPAPHOS);
    fOutputList->Add(fHistSigmaPAPHOSrot);
    fOutputList->Add(fHistSigmaAntiPAPHOS);
    fOutputList->Add(fHistSigmaPxResnopropPHOS);
    fOutputList->Add(fHistSigmaPyResnopropPHOS);
    fOutputList->Add(fHistSigmaPzResnopropPHOS);
    fOutputList->Add(fHistSigmaPxRespropPHOS);
    fOutputList->Add(fHistSigmaPyRespropPHOS);
    fOutputList->Add(fHistSigmaPzRespropPHOS);
    fOutputList->Add(fHistInvSigmaMassPHOSmix);
    fOutputList->Add(fHistSigmaYPHOSmix);
    fOutputList->Add(fHistSigmaPAPHOSmix);
    fOutputList->Add(fHistSigmaAntiPAPHOSmix);
    fOutputList->Add(fHistProtPhotonDCAPHOS);
    fOutputList->Add(fHistPhotonSecPAPHOS);
    fOutputList->Add(fHistPhotonSecPAPHOSMC);
    fOutputList->Add(fHistSigmaDCAtoPVPHOS);
    fOutputList->Add(fHistSigmaRadiusPHOS);
    fOutputList->Add(fHistSigmaDCAtoPVPHOSMC);
    fOutputList->Add(fHistProtPhotonDCAPHOSMC);
    fOutputList->Add(fHistNCaloPhotons);
    fOutputList->Add(fHistNClusters);
    fOutputList->Add(fHistClusterStatistics);
    fOutputList->Add(fHistClusterStatisticsMC);
    fOutputList->Add(fHistClusterStatisticsMCSig);
    fOutputList->Add(fHistClusterType);
    fOutputList->Add(fHistClusterTypeMC);
    fOutputList->Add(fHistClusterTypeMCSig);
    fOutputList->Add(fHistPHOSDisttoBC);
    fOutputList->Add(fHistPHOSM02);
    fOutputList->Add(fHistPHOSM20);
    fOutputList->Add(fHistPHOSTOF);
    fOutputList->Add(fHistPHOSBeta);
    fOutputList->Add(fHistPHOSE);
    fOutputList->Add(fHistPHOSNTracks);
    fOutputList->Add(fHistPHOSNCells);
    fOutputList->Add(fHistPHOSDx);
    fOutputList->Add(fHistPHOSDz);
    fOutputList->Add(fHistPHOSDy);
    fOutputList->Add(fHistPHOSDisp);
    fOutputList->Add(fHistPHOSMCDisttoBC);
    fOutputList->Add(fHistPHOSMCM02);
    fOutputList->Add(fHistPHOSMCM20);
    fOutputList->Add(fHistPHOSMCTOF);
    fOutputList->Add(fHistPHOSMCBeta);
    fOutputList->Add(fHistPHOSMCE);
    fOutputList->Add(fHistPHOSMCNTracks);
    fOutputList->Add(fHistPHOSMCNCells);
    fOutputList->Add(fHistPHOSMCDx);
    fOutputList->Add(fHistPHOSMCDz);
    fOutputList->Add(fHistPHOSMCDy);
    fOutputList->Add(fHistPHOSMCDisp);
    fOutputList->Add(fHistPHOSMCSigDisttoBC);
    fOutputList->Add(fHistPHOSMCSigM02);
    fOutputList->Add(fHistPHOSMCSigM20);
    fOutputList->Add(fHistPHOSMCSigTOF);
    fOutputList->Add(fHistPHOSMCSigBeta);
    fOutputList->Add(fHistPHOSMCSigE);
    fOutputList->Add(fHistPHOSMCSigNTracks);
    fOutputList->Add(fHistPHOSMCSigNCells);
    fOutputList->Add(fHistPHOSMCSigDx);
    fOutputList->Add(fHistPHOSMCSigDz);
    fOutputList->Add(fHistPHOSMCSigDy);
    fOutputList->Add(fHistPHOSMCSigDisp);
    fOutputList->Add(fHistEMCALDisttoBC);
    fOutputList->Add(fHistEMCALM02);
    fOutputList->Add(fHistEMCALM20);
    fOutputList->Add(fHistEMCALTOF);
    fOutputList->Add(fHistEMCALBeta);
    fOutputList->Add(fHistEMCALE);
    fOutputList->Add(fHistEMCALNTracks);
    fOutputList->Add(fHistEMCALNCells);
    fOutputList->Add(fHistEMCALDx);
    fOutputList->Add(fHistEMCALDz);
    fOutputList->Add(fHistEMCALDy);
    fOutputList->Add(fHistEMCALDisp);
    fOutputList->Add(fHistEMCALMCDisttoBC);
    fOutputList->Add(fHistEMCALMCM02);
    fOutputList->Add(fHistEMCALMCM20);
    fOutputList->Add(fHistEMCALMCTOF);
    fOutputList->Add(fHistEMCALMCBeta);
    fOutputList->Add(fHistEMCALMCE);
    fOutputList->Add(fHistEMCALMCNTracks);
    fOutputList->Add(fHistEMCALMCNCells);
    fOutputList->Add(fHistEMCALMCDx);
    fOutputList->Add(fHistEMCALMCDz);
    fOutputList->Add(fHistEMCALMCDy);
    fOutputList->Add(fHistEMCALMCDisp);
    fOutputList->Add(fHistEMCALMCSigDisttoBC);
    fOutputList->Add(fHistEMCALMCSigM02);
    fOutputList->Add(fHistEMCALMCSigM20);
    fOutputList->Add(fHistEMCALMCSigTOF);
    fOutputList->Add(fHistEMCALMCSigBeta);
    fOutputList->Add(fHistEMCALMCSigE);
    fOutputList->Add(fHistEMCALMCSigNTracks);
    fOutputList->Add(fHistEMCALMCSigNCells);
    fOutputList->Add(fHistEMCALMCSigDx);
    fOutputList->Add(fHistEMCALMCSigDz);
    fOutputList->Add(fHistEMCALMCSigDy);
    fOutputList->Add(fHistEMCALMCSigDisp);
    //Calculated photon
    fOutputList->Add(fHistCalcStatistics);
    fOutputList->Add(fHistCalcStatisticsMC);
    fOutputList->Add(fHistEqPar1Calc);
    fOutputList->Add(fHistEqPar1CalcMC);
    fOutputList->Add(fHistNIterCalc);
    fOutputList->Add(fHistNIterCalcMC);
    fOutputList->Add(fHistDeltaThetaCalc);
    fOutputList->Add(fHistDeltaThetaCalcMC);
    fOutputList->Add(fHistDeltaPhiCalc);
    fOutputList->Add(fHistDeltaPhiCalcMC);
    fOutputList->Add(fHistInvSigmaMassCalc);
    fOutputList->Add(fHistInvSigmaMassCalcMC);
    fOutputList->Add(fHistSigmaDCAtoPVCalc);
    fOutputList->Add(fHistSigmaDCAtoPVCalcMC);
    fOutputList->Add(fHistProtPhotonDCACalc);
    fOutputList->Add(fHistProtPhotonDCACalcMC);
    fOutputList->Add(fHistSigmaYCalc);
    fOutputList->Add(fHistSigmaYCalcMC);
    fOutputList->Add(fHistSigmaAntiPACalc);
    fOutputList->Add(fHistSigmaAntiPACalcMC);

/**************************************************************************/
    PostData(1, fOutputList);         
    PostData(2, fSigmaCandTree);        
    PostData(3, fSigmaPairTreeSE);      
    PostData(4, fSigmaPairTreeME);      
    PostData(5, fSigmaMEBackgroundTree);
    PostData(6, fSigmaPHOSCandTree);    
    PostData(7, fSigmaPairTreePHOSSE);  
    PostData(8, fSigmaPairTreePHOSME);  
    PostData(9, fSigmaPHOSMEBkgTree);   
    PostData(10,fSigmaCalcCandTree);   
    PostData(11,fProtonTree);   

}//end of UserCreateOutputObjects()

//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlus::UserExec(Option_t *)
{ 

  // Main loop. Called once for each event
  if(fDebug) cout << "Now in UserExec(). Start of Event\n";

  // Reset bools
  fEventhasSigma = kFALSE;
  fEventhasProton = kFALSE;
  fEventhasSigmaCand = kFALSE;

  // Load the Input Event and check it
  aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aodEvent) {
    AliWarning("ERROR: AOD Event not available!");
    return;
  }

  AliVHeader *aodEventHeader = aodEvent->GetHeader(); //Get the Header from the event
  if(!aodEventHeader){
    AliWarning("ERROR: Event Header not available!"); 
  }
  else fGlobalEventID = aodEventHeader->GetEventIdAsLong(); //Get global ID of the event

  mcEvent = MCEvent();           // Get MC event (called mcEvent) from the input file 
  if(mcEvent) isMonteCarlo = kTRUE; 
  else isMonteCarlo = kFALSE;

  if(isMonteCarlo) {
    Bool_t isSimPileup = kFALSE;
    AliAODMCHeader *aodMCHeader = (AliAODMCHeader*) (fInputEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(aodMCHeader){
      TString genname = (TString)aodMCHeader->GetGeneratorName();
      isSimPileup = AliAnalysisUtils::IsPileupInGeneratedEvent(aodMCHeader,genname);
    }
    FillHistogram("fHistMCGenPileup",isSimPileup);
    if(fRemoveGenPileup&&isSimPileup) return;

    AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!AODMCTrackArray)
	  { 
		  AliWarning("ERROR: MC Track Array not available!"); 
		  return; 
	  }
  }

  // Check the PID response
  if(!fPIDResponse) {
    AliError("ERROR: No pid response!");               
    return;
  }

  //Number of Tracks
  nTracks = aodEvent->GetNumberOfTracks();      
  FillHistogram("fHistTrackMultiplicity",nTracks);  
  if(nTracks==0) return; //No point in continuing if there are no tracks

  // Check primary vertex position
  Double_t primaryVtxPos[3] = {-999,-999,-999};
  const AliAODVertex *aodVtx = aodEvent->GetPrimaryVertex();
  
  if(!aodVtx) {
    AliWarning("No primary vertex in AOD!");
    return;
  }
  aodVtx->GetXYZ(primaryVtxPos);

  primaryVtxPosX=primaryVtxPos[0];
  primaryVtxPosY=primaryVtxPos[1];
  primaryVtxPosZ=primaryVtxPos[2];

  FillHistogram("fHistVertexZ",primaryVtxPosZ);

  if(isMonteCarlo){
    //Check MC primary vertex position
    Double_t primaryVtxPosMC[3] = {-999,-999,-999};
    const AliAODVertex *aodVtxMC = aodEvent->GetPrimaryVertex();
    if(!aodVtxMC) {AliWarning("No primary vertex in MC!");}
    else aodVtxMC->GetXYZ(primaryVtxPosMC);
    primaryVtxPosXMC=primaryVtxPosMC[0];
    primaryVtxPosYMC=primaryVtxPosMC[1];
    primaryVtxPosZMC=primaryVtxPosMC[2];
    FillHistogram("fHistVertexZMC",primaryVtxPosZMC);
  }

  if(isMonteCarlo){
    //Create Pseudo Event ID in case of MC
    Short_t ParA = aodEvent->GetNumberOfTracks();
    Short_t ParB = mcEvent->GetNumberOfTracks();
    UInt_t ParC = (*reinterpret_cast<unsigned int*>(&primaryVtxPosZMC));
    ULong64_t PseudoEventID = (((ULong64_t)ParA<<48)|((ULong64_t)ParB<<32)|((ULong64_t)ParC));
    fGlobalEventID = PseudoEventID;
  }

  //Magnetic Field
  Bz = aodEvent->GetMagneticField();    
  //Set Magnetic field for ALL KFParticles
  KFParticle::SetField(Bz);

  // Trigger
  AliVEventHandler* fInputHandler = (AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler();
  if(fInputHandler){
    EventTriggers = fInputHandler->IsEventSelected();
  }

  if(fDebug) cout << "Filling Event Counters\n";

  FillHistogram("fHistEventCounter",1); //Event Counter 
  TH1D* EventCounterdoub = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdouble"));
  if(!EventCounterdoub){AliWarning("Error: Histogram 'fHistEventCounterdouble' does not exist in TList 'fOutputList'!");} 
  else EventCounterdoub->Fill(1);

  //Check if Event has corresponding Trigger: Get Trigger Mask and AND it with requested trigger
  if((EventTriggers&65536)||(EventTriggers&8)){
    FillHistogram("fHistEventCounterHM",1); //Event Counter 
    TH1D* EventCounterdoubHM = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHM"));
    if(!EventCounterdoubHM){AliWarning("Error: Histogram 'fHistEventCounterdoubleHM' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHM->Fill(1);
  }
  if(EventTriggers&65536){
    FillHistogram("fHistEventCounterHMV0",1); //Event Counter 
    TH1D* EventCounterdoubHMV0 = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHMV0"));
    if(!EventCounterdoubHMV0){AliWarning("Error: Histogram 'fHistEventCounterdoubleHMV0' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHMV0->Fill(1);
  }
  if(EventTriggers&8){
    FillHistogram("fHistEventCounterHMSPD",1); //Event Counter 
    TH1D* EventCounterdoubHMSPD = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHMSPD"));
    if(!EventCounterdoubHMSPD){AliWarning("Error: Histogram 'fHistEventCounterdoubleHMSPD' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHMSPD->Fill(1);
  }
  if(EventTriggers&2){
    if(TMath::Abs(Bz)>4){
      FillHistogram("fHistEventCounterINT7",1); //Event Counter 
      TH1D* EventCounterdoubINT7 = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleINT7"));
      if(!EventCounterdoubINT7){AliWarning("Error: Histogram 'fHistEventCounterdoubleINT7' does not exist in TList 'fOutputList'!");} 
      else EventCounterdoubINT7->Fill(1);
    }
    if(TMath::Abs(Bz)<4){
      FillHistogram("fHistEventCounterINT7LF",1); //Event Counter 
      TH1D* EventCounterdoubINT7LF = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleINT7LF"));
      if(!EventCounterdoubINT7LF){AliWarning("Error: Histogram 'fHistEventCounterdoubleINT7LF' does not exist in TList 'fOutputList'!");} 
      else EventCounterdoubINT7LF->Fill(1);
    }
  }

  //Primary vertex cut
  if(std::abs(primaryVtxPosZ)>fMaxVertexZ) return;   //Return if vertex z position >10cm!                  

  FillHistogram("fHistMCCounter",1); //Event Counter 
  FillHistogram("fHistEventCounter",2); //Event Counter 
  if(!EventCounterdoub){AliWarning("Error: Histogram 'fHistEventCounterdouble' does not exist in TList 'fOutputList'!");} 
  else EventCounterdoub->Fill(2);

  if((EventTriggers&65536)||(EventTriggers&8)){
    FillHistogram("fHistEventCounterHM",2); //Event Counter 
    TH1D* EventCounterdoubHM = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHM"));
    if(!EventCounterdoubHM){AliWarning("Error: Histogram 'fHistEventCounterdoubleHM' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHM->Fill(2);
  }
  if(EventTriggers&65536){
    FillHistogram("fHistEventCounterHMV0",2); //Event Counter 
    TH1D* EventCounterdoubHMV0 = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHMV0"));
    if(!EventCounterdoubHMV0){AliWarning("Error: Histogram 'fHistEventCounterdoubleHMV0' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHMV0->Fill(2);
  }
  if(EventTriggers&8){
    FillHistogram("fHistEventCounterHMSPD",2); //Event Counter 
    TH1D* EventCounterdoubHMSPD = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleHMSPD"));
    if(!EventCounterdoubHMSPD){AliWarning("Error: Histogram 'fHistEventCounterdoubleHMSPD' does not exist in TList 'fOutputList'!");} 
    else EventCounterdoubHMSPD->Fill(2);
  }
  if(EventTriggers&2){
    if(TMath::Abs(Bz)>4){
      FillHistogram("fHistEventCounterINT7",2); //Event Counter 
      TH1D* EventCounterdoubINT7 = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleINT7"));
      if(!EventCounterdoubINT7){AliWarning("Error: Histogram 'fHistEventCounterdoubleINT7' does not exist in TList 'fOutputList'!");} 
      else EventCounterdoubINT7->Fill(2);
    }
    if(TMath::Abs(Bz)<4){
      FillHistogram("fHistEventCounterINT7LF",2); //Event Counter 
      TH1D* EventCounterdoubINT7LF = dynamic_cast<TH1D*>(fOutputList->FindObject("fHistEventCounterdoubleINT7LF"));
      if(!EventCounterdoubINT7LF){AliWarning("Error: Histogram 'fHistEventCounterdoubleINT7LF' does not exist in TList 'fOutputList'!");} 
      else EventCounterdoubINT7LF->Fill(2);
    }
  }

  // Centrality
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(aodEvent->FindListObject("MultSelection"));   //Get Mult Selection..
  if(multSelection) Centrality = multSelection->GetMultiplicityPercentile("V0M");                                //..and retrieve centrality
  FillHistogram("fHistCentrality",Centrality);

  if(EventTriggers&65536) FillHistogram("fHistCentralityHMV0",Centrality);
  if(EventTriggers&8) FillHistogram("fHistCentralityHMSPD",Centrality);
  if(EventTriggers&2){
    if(TMath::Abs(Bz)>4) FillHistogram("fHistCentralityINT7",Centrality);
    if(TMath::Abs(Bz)<4) FillHistogram("fHistCentralityINT7LF",Centrality);
  }

  AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  if(!aodHeader){
    AliWarning("ERROR: Event Header not available!"); 
  }
  else{
    fRefMultComb05 = aodHeader->GetRefMultiplicityComb05();
    fRefMultComb08 = aodHeader->GetRefMultiplicityComb08();
    fRefMultComb10 = aodHeader->GetRefMultiplicityComb10();
    FillHistogram("fHistRefMultComb05",fRefMultComb05);
    FillHistogram("fHistRefMultComb08",fRefMultComb08);
    FillHistogram("fHistRefMultComb10",fRefMultComb10);

    if(EventTriggers&65536){
      FillHistogram("fHistRefMultComb05HMV0",fRefMultComb05);
      FillHistogram("fHistRefMultComb08HMV0",fRefMultComb08);
      FillHistogram("fHistRefMultComb10HMV0",fRefMultComb10);
    }
    
    if(EventTriggers&8){
      FillHistogram("fHistRefMultComb05HMSPD",fRefMultComb05);
      FillHistogram("fHistRefMultComb08HMSPD",fRefMultComb08);
      FillHistogram("fHistRefMultComb10HMSPD",fRefMultComb10);
    }
    
    if(EventTriggers&65536||EventTriggers&8){
      FillHistogram("fHistRefMultComb05HM",fRefMultComb05);
      FillHistogram("fHistRefMultComb08HM",fRefMultComb08);
      FillHistogram("fHistRefMultComb10HM",fRefMultComb10);
    }

    if(EventTriggers&2){
      if(TMath::Abs(Bz)>4) FillHistogram("fHistRefMultComb05INT7",fRefMultComb05);
      if(TMath::Abs(Bz)>4) FillHistogram("fHistRefMultComb08INT7",fRefMultComb08);
      if(TMath::Abs(Bz)>4) FillHistogram("fHistRefMultComb10INT7",fRefMultComb10);
      if(TMath::Abs(Bz)<4) FillHistogram("fHistRefMultComb05INT7LF",fRefMultComb05);
      if(TMath::Abs(Bz)<4) FillHistogram("fHistRefMultComb08INT7LF",fRefMultComb08);
      if(TMath::Abs(Bz)<4) FillHistogram("fHistRefMultComb10INT7LF",fRefMultComb10);
    }
  }

  if(fDebug) cout << "Checking V0s\n";

  //Fill V0 arrays
  fOnFlyVector.clear(); //clear the arrays
  fFinderVector.clear();
  fV0ParticleIDArray.clear();

  Int_t nV0 = aodEvent->GetNumberOfV0s(); //Number of V0s in the event

  for (Int_t iV0=0; iV0<nV0; iV0++) {    //Loop over V0s in the event
      
    AliAODv0* aodV0 = (AliAODv0*)aodEvent->GetV0(iV0);  //Get V0 object
    if(!aodV0) continue;

    // Check basic V0 properties: 2 Daughters, opposite charge, total charge = 0
    if(aodV0->GetNDaughters() != 2)                    continue;
    if(aodV0->GetNProngs() != 2)                       continue;
    if(aodV0->GetCharge() != 0)                        continue;
    if(aodV0->ChargeProng(0) == aodV0->ChargeProng(1)) continue;

    // Get daughter tracks      
    AliAODTrack* trackN = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
    AliAODTrack* trackP = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
    if (!trackN||!trackP) continue;
    if (trackN->GetSign() == trackP->GetSign()) continue;

    fFinderVector.push_back(trackN->GetID());
    fFinderVector.push_back(trackP->GetID());

    if(aodV0->GetOnFlyStatus()){
      fOnFlyVector.push_back(trackN->GetID());
      fOnFlyVector.push_back(trackP->GetID());
    }
  }//End of V0 Loop. Finished preparing the maps

/************************Start Event Processing**************************************/

  //Process MC Particles
  if(fDebug&&isMonteCarlo&&fProcessMCParticles) cout << "Processing MC Particles\n";
  if(isMonteCarlo && fProcessMCParticles) ProcessMCParticles();

  //Process Protons
  if(fDebug&&fProcessProtons) cout << "Processing Protons\n";
  if(fProcessProtons) FillProtonArray();

  //Process V0s
  if(fDebug&&fProcessV0s) cout << "Processing V0s\n";
  if(fProcessV0s) FillV0PhotonArray();

  //Process Calo Clusters
  if(fDebug&&fProcessClusters) cout << "Processing Calo Clusters\n";
  if(fProcessClusters) FillCaloClusterArray();

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessReco) cout << "Reconstructing Pi0s and Sigmas\n";
  if(fProcessReco) ReconstructParticles();

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessRecoPHOS) cout << "Reconstructing Pi0s and Sigmas. One PHOS Photon\n";
  if(fProcessRecoPHOS) ReconstructParticlesPHOS();

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessRecoCalc) cout << "Reconstructing Pi0s and Sigmas. One calculated Photon\n";
  if(fProcessRecoCalc) ReconstructParticlesCalc();

  //Save Protons
  if(fDebug&&fSaveProtons) cout << "Saving Protons\n";
  if(fSaveProtons&&fEventhasSigmaCand) FillProtonTree();

  //Update Event Pools at the End of the Event if Event Mixing is enabled
  if(fSaveMixedBackground){
    if(!fRequireSigmaCand||fEventhasSigmaCand){
      if(fDebug) cout << "Updating Event Pool 1\n";
   
        //Get Pool from Pool Manager for given RefMult and Z Vertex values
        AliEventPool* Evpool = 0x0;
	      if(fEvPoolMgr&&fUseAbsZ)  Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
	      if(fEvPoolMgr&&!fUseAbsZ) Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
	      if(Evpool){
    
        //Create TObjArray of selected Protons
        TObjArray* ProtonObjArray = new TObjArray();    
  	    ProtonObjArray->SetOwner(kTRUE);
        Int_t nProtonforMixing = fProtonArray.size();
        for(Int_t k=0; k<nProtonforMixing; k++) {
          AliAODTrack *prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(k));
          if(!prot) continue;
          AliAODTrackreduced* redprot = new AliAODTrackreduced();
          if(!redprot) continue;
          redprot->InitfromTrack(prot, fPIDResponse, primaryVtxPosX, primaryVtxPosY, primaryVtxPosZ);
          ProtonObjArray->Add(redprot);
        }

        //Clone it and update the Pool    
        TObjArray* ProtonCloneArray = (TObjArray*)ProtonObjArray->Clone();    
        ProtonCloneArray->SetOwner(kTRUE);
        Evpool->UpdatePool(ProtonCloneArray);
      }
      else{AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ));}
    } //Require Sigma Candidate
  } //End of Pool updating

  if(fSavePHOSMixedBackground){
    if(!fRequireSigmaCand||fEventhasSigmaCand){

      if(fDebug) cout << "Updating Event Pool 3\n";

      //Get Pool from Pool Manager for given RefMult and Z Vertex values
      AliEventPool* Evpool = 0x0;
  	  if(fEvPoolMgr3&&fUseAbsZ)  Evpool = fEvPoolMgr3->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
  	  if(fEvPoolMgr3&&!fUseAbsZ) Evpool = fEvPoolMgr3->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
  	  if(Evpool){
      
        //Create TObjArray of selected Clusters
        TObjArray* ClusterObjArray = new TObjArray();    
    	  ClusterObjArray->SetOwner(kTRUE);
        Int_t nClusterforMixing = fCaloPhotonArray.size();
        for(Int_t k=0; k<nClusterforMixing; k++) {
          AliAODCaloCluster *clust = (AliAODCaloCluster*)aodEvent->GetCaloCluster(fCaloPhotonArray.at(k));
          if(!clust) continue;
          AliAODCaloCells* calocells;
          if(clust->IsPHOS()) calocells = (AliAODCaloCells*)aodEvent->GetPHOSCells();
          if(clust->IsEMCAL()) calocells = (AliAODCaloCells*)aodEvent->GetEMCALCells();
          AliAODClusterreduced* redclust = new AliAODClusterreduced();
          if(!redclust) continue;          
          redclust->InitfromCluster(clust,calocells);
          ClusterObjArray->Add(redclust);
        }

        //Clone it and update the Pool    
        TObjArray* ClusterCloneArray = (TObjArray*)ClusterObjArray->Clone();    
        ClusterCloneArray->SetOwner(kTRUE);
        Evpool->UpdatePool(ClusterCloneArray);
      }
      else{AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ));}
    } //Require Sigma Candidate
  } //End of Pool updating

  if(fFillPairTreeME||fFillPHOSPairTreeME){  
    if(!fRequireProton||fEventhasProton){
      if(!fRequireSigma||fEventhasSigma){
      
        if(fDebug) cout << "Updating Event Pool 2\n";

        //Get Pool from Pool Manager for given RefMult and Z Vertex values
        AliEventPool* Evpool = 0x0;
    	  if(fEvPoolMgr2&&fUseAbsZCorr)  Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
    	  if(fEvPoolMgr2&&!fUseAbsZCorr) Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
    	  if(Evpool){
        
          //Create TObjArray of selected Protons
          TObjArray* ProtonObjArray = new TObjArray();    
      	  ProtonObjArray->SetOwner(kTRUE);
          Int_t nProtonforMixing = fProtonArray2.size();
          for(Int_t k=0; k<nProtonforMixing; k++) {
            AliAODTrack *prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray2.at(k));
            if(!prot) continue;
            AliAODTrackcorrelation* redprot = new AliAODTrackcorrelation();
            if(!redprot) continue;
            redprot->InitfromTrack(prot, fPIDResponse);
            ProtonObjArray->Add(redprot);
          }

          //Clone it and update the Pool    
          TObjArray* ProtonCloneArray = (TObjArray*)ProtonObjArray->Clone();    
          ProtonCloneArray->SetOwner(kTRUE);
          Evpool->UpdatePool(ProtonCloneArray);
        }
        else{AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ));}
      } //Require Sigma
    } //Require Proton
  } //End of Pool updating

  if(fDebug) cout << "End of Event processing. Calling PostData()\n";

/************************End of Event Processing**************************************/

  PostData(1, fOutputList); // stream the analysis results of the current event to output manager               
  PostData(2, fSigmaCandTree);        
  PostData(3, fSigmaPairTreeSE);      
  PostData(4, fSigmaPairTreeME);      
  PostData(5, fSigmaMEBackgroundTree);
  PostData(6, fSigmaPHOSCandTree);    
  PostData(7, fSigmaPairTreePHOSSE);  
  PostData(8, fSigmaPairTreePHOSME);  
  PostData(9, fSigmaPHOSMEBkgTree);   
  PostData(10,fSigmaCalcCandTree);   
  PostData(11,fProtonTree);   

  if(fDebug) cout << "Returning from UserExec()\n";
  return;

}//end of UserExec()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillProtonArray(){

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fProtonArray.clear();
  fProtonArray2.clear();
  Int_t countProton = 0;

  //Loop for Proton Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {

    //Initialisation of local variables
    Bool_t isReallyProton        = kFALSE; //From MC
    Bool_t isProtonfromSigma     = kFALSE; //From MC
    Bool_t isPrimarySigma        = kFALSE; //From MC
    Bool_t isPrimaryProton       = kFALSE; //From MC
    Bool_t isProtonfromLambda    = kFALSE; //From MC
    Bool_t isProtonfromWeakDecay = kFALSE; //From MC
    Bool_t isProtonfromMaterial  = kFALSE; //From MC

    Bool_t   isTPCProton = kFALSE;
    Bool_t   isTOFProton = kFALSE;

    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack) {
      AliWarning("No AOD Track!");
      continue;
    }

    //Check MC Truth
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        if(mcPart->GetPdgCode()==2212 || mcPart->GetPdgCode()==-2212){
          isReallyProton = kTRUE;
          if(mcPart->IsPhysicalPrimary()) isPrimaryProton = kTRUE;
          if(mcPart->IsSecondaryFromWeakDecay()) isProtonfromWeakDecay = kTRUE;
          if(mcPart->IsSecondaryFromMaterial()) isProtonfromMaterial = kTRUE;
          if(mcPart->GetMother()!=-1){
            AliAODMCParticle* ProtonMotherPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mcPart->GetMother()));
            if(ProtonMotherPart){
              if(ProtonMotherPart->GetPdgCode()==3222 || ProtonMotherPart->GetPdgCode()==-3222){ 
                isProtonfromSigma = kTRUE;
                if(ProtonMotherPart->IsPhysicalPrimary()) isPrimarySigma = kTRUE;
              }//Is really Sigma
              if(ProtonMotherPart->GetPdgCode()==2114) isProtonfromLambda = kTRUE;
            }//MC Mother exists
          }//Proton has Mother
        }//is Proton or Anti-Proton
      }//MC Particle exists
    }//MC treatment

    FillHistogram("fHistProtonStatistics",0);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",0);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",0);

    //Reject Tracks with negative ESD ID (Rejected anyways, Global Tracks have positive IDs)
    if(aodTrack->GetID()<0&&fRejectNegIDs) continue; 

    FillHistogram("fHistProtonStatistics",1);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",1);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",1);

    //Reject Tracks with Filterbit 0 (These Tracks have low quality and are only stored because the are used by the V0 finder)
    if(!aodTrack->GetFilterMap()&&fRejectZeroFilterBit) continue; 

    FillHistogram("fHistProtonStatistics",2);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",2);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",2);

    Double_t eta    = aodTrack->Eta();
    Double_t phi    = aodTrack->Phi();
    Double_t dEdx   = aodTrack->GetTPCsignal();
    Double_t p      = aodTrack->GetTPCmomentum();
    Double_t pt     = aodTrack->Pt();
    Int_t    nTPCClustProt = aodTrack->GetTPCNcls();
    Int_t    nITSClustProt = aodTrack->GetITSNcls();
    Double_t chi2   = aodTrack->GetTPCchi2();

    Double_t len = aodTrack->GetIntegratedLength();
    Double_t tim;
    if(aodTrack->GetTPCmomentum()!=0) tim = aodTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(aodTrack->GetTPCmomentum());
    else tim = aodTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(TMath::Sqrt(aodTrack->Px()*aodTrack->Px()+aodTrack->Py()*aodTrack->Py()+aodTrack->Pz()*aodTrack->Pz()));
    Double_t beta = -1;
    if(tim != 0.) beta = len / (tim * c);

    Double_t nSigmaTPCProt = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton);
    Double_t nSigmaTOFProt = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton);

    if(TMath::Abs(nSigmaTPCProt)<fMaxNsigProtTPC) isTPCProton = kTRUE;
    if(TMath::Abs(nSigmaTOFProt)<fMaxNsigProtTOF) isTOFProton = kTRUE;

    Float_t  DCAxy = -999., DCAz = -999.;
    aodTrack->GetImpactParameters(DCAxy,DCAz);
    DCAxy = TMath::Abs(DCAxy); DCAz = TMath::Abs(DCAz);

    //QA Histograms
    FillHistogram("fHistTrackEtaPhi",eta,phi);
    if(isReallyProton) FillHistogram("fHistProtonEtaPhiMC",eta,phi);

    //Acceptance Cut
    if(TMath::Abs(eta) > fMaxProtEta) continue; 

    FillHistogram("fHistProtonStatistics",3);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",3);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",3);

    //Topology
    FillHistogram("fHistProtonDCAxy",DCAxy);
    FillHistogram("fHistProtonDCAz",DCAz);

    //Track Quality
    FillHistogram("fHistTrackChi2",chi2);
    FillHistogram("fHistTrackTPCCluster",nTPCClustProt);
    FillHistogram("fHistTrackITSCluster",nITSClustProt);
    FillHistogram("fHistTrackpvsdEdx",p,dEdx);
    FillHistogram("fHistTrackpvsbeta",p,beta);
    FillHistogram("fHistTrackpt",pt);

    //Proton PID
    FillHistogram("fHistProtonpvsNSigmaTPC",p,nSigmaTPCProt);
    FillHistogram("fHistProtonpvsNSigmaTOF",p,nSigmaTOFProt);

    //MC Histograms
    if(isReallyProton){
      FillHistogram("fHistProtonChi2MC",chi2);
      FillHistogram("fHistProtonDCAxyMC",DCAxy);
      FillHistogram("fHistProtonDCAzMC",DCAz);
      FillHistogram("fHistProtonTPCClusterMC",nTPCClustProt);
      FillHistogram("fHistProtonITSClusterMC",nITSClustProt);
      FillHistogram("fHistProtonpvsNSigmaTPCMC",p,nSigmaTPCProt);
      FillHistogram("fHistProtonpvsNSigmaTOFMC",p,nSigmaTOFProt);
      FillHistogram("fHistProtonptMC",pt);
    }

    if(isProtonfromSigma){
      FillHistogram("fHistProtonDCAxyMCSigma",DCAxy);
      FillHistogram("fHistProtonDCAzMCSigma",DCAz);
    }

    if(isPrimarySigma){
      FillHistogram("fHistProtonDCAxyMCPrimSig",DCAxy);
      FillHistogram("fHistProtonDCAzMCPrimSig",DCAz);
    }

    if(isPrimaryProton){
    FillHistogram("fHistPrimProtonDCAxyMC",DCAxy);
    FillHistogram("fHistPrimProtonDCAzMC",DCAz);
    }

    if(isProtonfromMaterial){
    FillHistogram("fHistMaterialProtonDCAxyMC",DCAxy);
    FillHistogram("fHistMaterialProtonDCAzMC",DCAz);
    }

    if(isProtonfromWeakDecay){
    FillHistogram("fHistWeakProtonDCAxyMC",DCAxy);
    FillHistogram("fHistWeakProtonDCAzMC",DCAz);
    }

    if(isProtonfromLambda){
    FillHistogram("fHistLambdaProtonDCAxyMC",DCAxy);
    FillHistogram("fHistLambdaProtonDCAzMC",DCAz);
    }

    // Cluster Cut
    if(nTPCClustProt<fMinTPCClustProt) continue;

    FillHistogram("fHistProtonStatistics",4);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",4);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",4);

    // TPC PID Cut. Continue if TPC PID is forced. Else require at least TOF PID
    if(!isTPCProton&&fRequireProtonTPC) continue;
    else if(!isTPCProton&&!isTOFProton) continue;

    FillHistogram("fHistProtonStatistics",5);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",5);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",5);

    // TOF PID Cut. Continue if TOF PID is forced and p > fMaxpOnlyTPCPID
    if(p > fMaxpOnlyTPCPID && !isTOFProton && fRequireProtonTOF) continue;
    if(p > fMaxpOnlyTPCPID && !isTOFProton && nSigmaTOFProt!=-999) continue;

    FillHistogram("fHistProtonStatistics",6);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",6);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",6);

    //Histograms for Purity and Efficiency Calcultation
    if(isReallyProton) FillHistogram("fHistProtonptwCutsMC",pt);
    FillHistogram("fHistProtonptwCuts",pt);
    
    // Kinematic Cut
    if(pt < fMinProtpt || pt > fMaxProtpt) continue;

    FillHistogram("fHistProtonStatistics",7);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",7);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",7);

    // Store (Anti-)Proton candidates after selection
    fProtonArray.push_back(iTrack);

    //DCA to PV Cut
    if(DCAxy>fMaxCorrPairProtonDCAxy) continue;
    if(DCAz>fMaxCorrPairProtonDCAz) continue;

    // TOF PID Cut. Continue if TOF PID is forced and p > fMaxpOnlyTPCPID
    if(p > fMaxpOnlyTPCPID && !isTOFProton && fRequireProtonTOFforPairs) continue;

    fProtonArray2.push_back(iTrack);

    //Acceptance Cut
    if(TMath::Abs(eta) > fStrictMaxProtEta) continue; 
    // Cluster Cut
    if(nTPCClustProt<fStrictMinTPCClustProt) continue;
    // TPC PID Cut. Continue if TPC PID is forced. Else require at least TOF PID
    if(TMath::Abs(nSigmaTPCProt)>fStrictMaxNsigProtTPC) continue;
    // TOF PID Cut. Continue if TOF PID is forced and p > fMaxpOnlyTPCPID
    if(p > fStrictMaxpOnlyTPCPID && fRequireProtonTOFforPairs && TMath::Abs(nSigmaTOFProt)>fStrictMaxNsigProtTOF) continue;
    // Kinematic Cut
    if(pt < fStrictMinProtpt || pt > fStrictMaxProtpt) continue;

    countProton++;

  } // End of proton track loop

  if(countProton!=0) fEventhasProton = kTRUE;

return;

}//End of FillProtonArray  

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillProtonTree(){

  //Loop for Proton Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {

    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack){
      AliWarning("No AOD Track!");
      continue;
    }

    //Reject Tracks with negative ESD ID (Rejected anyways, Global Tracks have positive IDs)
    if(aodTrack->GetID()<0&&fRejectNegIDs) continue; 

    //Reject Tracks with Filterbit 0 (These Tracks have low quality and are only stored because the are used by the V0 finder)
    if(!aodTrack->GetFilterMap()&&fRejectZeroFilterBit) continue; 

    //Acceptance Cut
    if(TMath::Abs(aodTrack->Eta()) > fMaxProtEta) continue; 

    // Cluster Cut
    if(aodTrack->GetTPCNcls()<fMinTPCClustProt) continue;

    Float_t  DCAxy = -999., DCAz = -999.;
    aodTrack->GetImpactParameters(DCAxy,DCAz);

    //DCA to PV Cut
    if(TMath::Abs(DCAxy)>fMaxCorrPairProtonDCAxy) continue;
    if(TMath::Abs(DCAz)>fMaxCorrPairProtonDCAz) continue;

    Bool_t isReallyProton  = kFALSE;
    Int_t  ParticlePdg = -999;
    Bool_t isPrimary = kFALSE;
    Float_t MCPx = -999;
    Float_t MCPy = -999;
    Float_t MCPz = -999;

    //Check MC Truth
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        ParticlePdg = mcPart->GetPdgCode();
        if(mcPart->IsPhysicalPrimary()) isPrimary = kTRUE;
        MCPx = mcPart->Px();
        MCPy = mcPart->Py();
        MCPz = mcPart->Pz();
        if(TMath::Abs(ParticlePdg)==2212){isReallyProton = kTRUE;}//is Proton or Anti-Proton
      }//MC Particle exists
    }//MC treatment

    Float_t nSigmaTPCProt = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton);
    Float_t nSigmaTPCPion = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion);
    Float_t nSigmaTPCKaon = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon);
    Float_t nSigmaTPCElec = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kElectron);
    Float_t nSigmaTOFProt = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton);
    Float_t nSigmaTOFPion = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kPion);
    Float_t nSigmaTOFKaon = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kKaon);
    Float_t nSigmaTOFElec = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kElectron);

    Bool_t isTPCProt = (TMath::Abs(nSigmaTPCProt)<fMaxNsigProtTPC);
    Bool_t isTPCPion = (TMath::Abs(nSigmaTPCPion)<fMaxNsigPionTPC);
    Bool_t isTPCKaon = (TMath::Abs(nSigmaTPCKaon)<fMaxNsigKaonTPC);
    Bool_t isTPCElec = (TMath::Abs(nSigmaTPCElec)<fMaxNsigElecTPC);
    Bool_t isTOFProt = (TMath::Abs(nSigmaTOFProt)<fMaxNsigProtTOF);
    Bool_t isTOFPion = (TMath::Abs(nSigmaTOFPion)<fMaxNsigPionTOF);
    Bool_t isTOFKaon = (TMath::Abs(nSigmaTOFKaon)<fMaxNsigKaonTOF);
    Bool_t isTOFElec = (TMath::Abs(nSigmaTOFElec)<fMaxNsigElecTOF);

    Bool_t isTPCOther = (isTPCPion||isTPCKaon||isTPCElec);
    Bool_t isTOFOther = (isTOFPion||isTOFKaon||isTOFElec);

    // TPC PID Cut
    if(!isTPCProt) continue;

    // if p > 0.x require TOF or reject other TPC bands
    if(aodTrack->GetTPCmomentum() > fMaxpOnlyTPCPID && !isTOFProt && isTPCOther) continue;

    // Reject other TOF bands
    if(aodTrack->GetTPCmomentum() > fMaxpOnlyTPCPID && isTOFOther && isTPCOther && fRejOtherTOF) continue;

    fSigRunnumber = aodEvent->GetRunNumber();        
    fSigTriggerMask = EventTriggers;      
    fSigEventID = fGlobalEventID;          
    fSigCentrality = Centrality;       
    fSigBField = Bz;           
    fSigRefMultComb05 = fRefMultComb05;    
    fSigRefMultComb08 = fRefMultComb08;    
    fSigRefMultComb10 = fRefMultComb10;    
    fPrimVertX = primaryVtxPosX;           
    fPrimVertY = primaryVtxPosY;           
    fPrimVertZ = primaryVtxPosZ;           
    fPairProtonID = aodTrack->GetID();        
    fPairProtonStatus = aodTrack->GetStatus();    
    fPairProtonFilterMap = aodTrack->GetFilterMap(); 
    fPairProtonCharge = aodTrack->Charge();    
    fPairProtonPx = aodTrack->Px();        
    fPairProtonPy = aodTrack->Py();        
    fPairProtonPz = aodTrack->Pz();        
    fPairProtonP = aodTrack->GetTPCmomentum();         
    fPairProtonEta = aodTrack->Eta();       
    fPairProtonDCAtoPVxy = DCAxy; 
    fPairProtonDCAtoPVz = DCAz;  
    fPairProtonNSigTPC = nSigmaTPCProt;   
    fPairProtNSigTPCKaon = nSigmaTPCPion; 
    fPairProtNSigTPCPion = nSigmaTPCKaon; 
    fPairProtNSigTPCElec = nSigmaTPCElec; 
    fPairProtonNSigTOF = nSigmaTOFProt;   
    fPairProtNSigTOFKaon = nSigmaTOFPion; 
    fPairProtNSigTOFPion = nSigmaTOFKaon; 
    fPairProtNSigTOFElec = nSigmaTOFElec; 
    fPairProtonCluster = aodTrack->GetTPCNcls();   
    fPairProtonITSCluster = aodTrack->GetITSNcls();
    fPairProtonChi2 = aodTrack->GetTPCchi2();
    fPrimVertXMC = primaryVtxPosXMC;         
    fPrimVertYMC = primaryVtxPosYMC;         
    fPrimVertZMC = primaryVtxPosZMC;         
    fPairProtonPxMC = MCPx;
    fPairProtonPyMC = MCPy;
    fPairProtonPzMC = MCPz;
    fPairProtonIsMC = isReallyProton;  
    fPairProtonPDG = ParticlePdg;    
    fPairProtonIsPrimary = isPrimary; 

    fProtonTree->Fill();

  } // End of proton track loop

return;

}//End of FillProtonTree  

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ProcessMCParticles() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  if(!isMonteCarlo) return;

  Int_t nMCTracks = mcEvent->GetNumberOfTracks();

  //loop over all MC tracks
  for(Int_t iMCtrack = 0; iMCtrack < nMCTracks; iMCtrack++){

    AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iMCtrack));
    if(!mcPart) {
      AliWarning("Could not retrieve AOD MC Particle!");
      continue;
    }

    Int_t MCPartPDGCode = mcPart->PdgCode(); 
    Double_t Rapidity = mcPart->Y();

    if(mcPart->IsPhysicalPrimary()){
        if(MCPartPDGCode == 3222){
            FillHistogram("fHistMCPrimSigmaPtvsRap",mcPart->Pt(),Rapidity);
            if(TMath::Abs(Rapidity)<0.5) FillHistogram("fHistMCPrimSigmaPtRap05",mcPart->Pt());
            if(TMath::Abs(Rapidity)<0.8) FillHistogram("fHistMCPrimSigmaPtRap08",mcPart->Pt());
            if(TMath::Abs(Rapidity)<1) FillHistogram("fHistMCPrimSigmaPtRap10",mcPart->Pt());
        }
        if(MCPartPDGCode == -3222){
            FillHistogram("fHistMCPrimAntiSigmaPtvsRap",mcPart->Pt(),Rapidity);
            if(TMath::Abs(Rapidity)<0.5) FillHistogram("fHistMCPrimAntiSigmaPtRap05",mcPart->Pt());
            if(TMath::Abs(Rapidity)<0.8) FillHistogram("fHistMCPrimAntiSigmaPtRap08",mcPart->Pt());
            if(TMath::Abs(Rapidity)<1) FillHistogram("fHistMCPrimAntiSigmaPtRap10",mcPart->Pt());
        }
    }

    if(MCPartPDGCode == 3222){
        FillHistogram("fHistMCSigmaPtvsRap",mcPart->Pt(),Rapidity);
        if(TMath::Abs(Rapidity)<0.5) FillHistogram("fHistMCSigmaPtRap05",mcPart->Pt());
        if(TMath::Abs(Rapidity)<0.8) FillHistogram("fHistMCSigmaPtRap08",mcPart->Pt());
        if(TMath::Abs(Rapidity)<1) FillHistogram("fHistMCSigmaPtRap10",mcPart->Pt());
    }
    if(MCPartPDGCode == -3222){
        FillHistogram("fHistMCAntiSigmaPtvsRap",mcPart->Pt(),Rapidity);
        if(TMath::Abs(Rapidity)<0.5) FillHistogram("fHistMCAntiSigmaPtRap05",mcPart->Pt());
        if(TMath::Abs(Rapidity)<0.8) FillHistogram("fHistMCAntiSigmaPtRap08",mcPart->Pt());
        if(TMath::Abs(Rapidity)<1) FillHistogram("fHistMCAntiSigmaPtRap10",mcPart->Pt());
    }

    if(MCPartPDGCode == 3222) FillHistogram("fHistMCSigmaPtvsEta",mcPart->Pt(),TMath::Abs(mcPart->Eta())); 
    if(MCPartPDGCode == -3222) FillHistogram("fHistMCAntiSigmaPtvsEta",mcPart->Pt(),TMath::Abs(mcPart->Eta()));

    if(mcPart->IsPhysicalPrimary()){
      if(MCPartPDGCode == 3222) FillHistogram("fHistMCPrimSigmaPtvsEta",mcPart->Pt(),TMath::Abs(mcPart->Eta())); 
      if(MCPartPDGCode == -3222) FillHistogram("fHistMCPrimAntiSigmaPtvsEta",mcPart->Pt(),TMath::Abs(mcPart->Eta()));
      if(MCPartPDGCode == 3222 && TMath::Abs(mcPart->Eta())<0.8) FillHistogram("fHistMCPrimSigmaPt08",mcPart->Pt()); 
      if(MCPartPDGCode == -3222 && TMath::Abs(mcPart->Eta())<0.8) FillHistogram("fHistMCPrimAntiSigmaPt08",mcPart->Pt());
      if(MCPartPDGCode == 3222 && TMath::Abs(mcPart->Eta())<1) FillHistogram("fHistMCPrimSigmaPt10",mcPart->Pt()); 
      if(MCPartPDGCode == -3222 && TMath::Abs(mcPart->Eta())<1) FillHistogram("fHistMCPrimAntiSigmaPt10",mcPart->Pt());
    }

    if(TMath::Abs(mcPart->Eta())>fMaxMCEta) continue; //Acceptance Cut

    //MC Particle Counter: Protons, Deltas, Sigma, Pi0, etc.
    FillHistogram("fHistMCCounter",2); //All Particles
    
    if(MCPartPDGCode == 2212)  {FillHistogram("fHistMCCounter",3); 
      FillHistogram("fHistMCProtonPt",mcPart->Pt()); 
      if(mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimProtonPt",mcPart->Pt());}
    }
    if(MCPartPDGCode == -2212) {FillHistogram("fHistMCCounter",4); 
      FillHistogram("fHistMCAntiProtonPt",mcPart->Pt()); 
      if(mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimAntiProtonPt",mcPart->Pt());}
    }

    if(MCPartPDGCode == 2214)  {FillHistogram("fHistMCCounter",5); FillHistogram("fHistMCDeltaPt",mcPart->Pt());}
    if(MCPartPDGCode == -2214) {FillHistogram("fHistMCCounter",6); FillHistogram("fHistMCAntiDeltaPt",mcPart->Pt());}
    
    if(MCPartPDGCode == 3222) {FillHistogram("fHistMCCounter",7); 
      FillHistogram("fHistSigmaMotherPart",1);
      FillHistogram("fHistMCSigmaPt",mcPart->Pt()); 
      FillHistogram("fHistMCSigmaOrigin",TMath::Sqrt(mcPart->Xv()*mcPart->Xv()+mcPart->Yv()*mcPart->Yv())); 
      if(mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimSigmaPt",mcPart->Pt()); FillHistogram("fHistSigmaMotherPart",2); continue;}
    
      AliAODMCParticle* SigmaMother = NULL;
      Int_t SigmaMotherPdg = 0;
      if(mcPart->GetMother()!=-1) SigmaMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      if(SigmaMother) SigmaMotherPdg = SigmaMother->GetPdgCode();

      if(SigmaMotherPdg==0) continue;

      if(SigmaMotherPdg==130) FillHistogram("fHistSigmaMotherPart",3);
      else if(SigmaMotherPdg==310) FillHistogram("fHistSigmaMotherPart",4);
      else if(SigmaMotherPdg==311) FillHistogram("fHistSigmaMotherPart",5);
      else if(SigmaMotherPdg==-311) FillHistogram("fHistSigmaMotherPart",6);
      else if(SigmaMotherPdg==-321) FillHistogram("fHistSigmaMotherPart",7);
      else if(SigmaMotherPdg==321) FillHistogram("fHistSigmaMotherPart",8);
      else if(SigmaMotherPdg==-211) FillHistogram("fHistSigmaMotherPart",9);
      else if(SigmaMotherPdg==211) FillHistogram("fHistSigmaMotherPart",10);
      else if(SigmaMotherPdg==2112) FillHistogram("fHistSigmaMotherPart",11);
      else if(SigmaMotherPdg==-2112) FillHistogram("fHistSigmaMotherPart",12);
      else if(SigmaMotherPdg==2212) FillHistogram("fHistSigmaMotherPart",13);
      else if(SigmaMotherPdg==-2212) FillHistogram("fHistSigmaMotherPart",14);
      else if(SigmaMotherPdg==3122) FillHistogram("fHistSigmaMotherPart",15);
      else if(SigmaMotherPdg==-3122) FillHistogram("fHistSigmaMotherPart",16);
      else if(SigmaMotherPdg==3222) FillHistogram("fHistSigmaMotherPart",17);
      else if(SigmaMotherPdg==3112) FillHistogram("fHistSigmaMotherPart",18);
      else if(SigmaMotherPdg==-3222) FillHistogram("fHistSigmaMotherPart",19);
      else if(SigmaMotherPdg==-3112) FillHistogram("fHistSigmaMotherPart",20);
      else if(SigmaMotherPdg==3322) FillHistogram("fHistSigmaMotherPart",21);
      else if(SigmaMotherPdg==-3322) FillHistogram("fHistSigmaMotherPart",22);
      else if(SigmaMotherPdg==3312) FillHistogram("fHistSigmaMotherPart",23);
      else if(SigmaMotherPdg==-3312) FillHistogram("fHistSigmaMotherPart",24);
      else FillHistogram("fHistSigmaMotherPart",25);
    
    } //End of SigmaPlus treatment
    
    if(MCPartPDGCode == -3222) {FillHistogram("fHistMCCounter",8); 
      FillHistogram("fHistAntiSigmaMotherPart",1); 
      FillHistogram("fHistMCAntiSigmaPt",mcPart->Pt()); 
      FillHistogram("fHistMCAntiSigmaOrigin",TMath::Sqrt(mcPart->Xv()*mcPart->Xv()+mcPart->Yv()*mcPart->Yv()));       
      if(mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimAntiSigmaPt",mcPart->Pt()); FillHistogram("fHistAntiSigmaMotherPart",2); continue;}

      AliAODMCParticle* SigmaMother = NULL;
      Int_t SigmaMotherPdg = 0;
      if(mcPart->GetMother()!=-1) SigmaMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      if(SigmaMother) SigmaMotherPdg = SigmaMother->GetPdgCode();

      if(SigmaMotherPdg==0) continue;

      if(SigmaMotherPdg==130) FillHistogram("fHistAntiSigmaMotherPart",3);
      else if(SigmaMotherPdg==310) FillHistogram("fHistAntiSigmaMotherPart",4);
      else if(SigmaMotherPdg==311) FillHistogram("fHistAntiSigmaMotherPart",5);
      else if(SigmaMotherPdg==-311) FillHistogram("fHistAntiSigmaMotherPart",6);
      else if(SigmaMotherPdg==-321) FillHistogram("fHistAntiSigmaMotherPart",7);
      else if(SigmaMotherPdg==321) FillHistogram("fHistAntiSigmaMotherPart",8);
      else if(SigmaMotherPdg==-211) FillHistogram("fHistAntiSigmaMotherPart",9);
      else if(SigmaMotherPdg==211) FillHistogram("fHistAntiSigmaMotherPart",10);
      else if(SigmaMotherPdg==2112) FillHistogram("fHistAntiSigmaMotherPart",11);
      else if(SigmaMotherPdg==-2112) FillHistogram("fHistAntiSigmaMotherPart",12);
      else if(SigmaMotherPdg==2212) FillHistogram("fHistAntiSigmaMotherPart",13);
      else if(SigmaMotherPdg==-2212) FillHistogram("fHistAntiSigmaMotherPart",14);
      else if(SigmaMotherPdg==3122) FillHistogram("fHistAntiSigmaMotherPart",15);
      else if(SigmaMotherPdg==-3122) FillHistogram("fHistAntiSigmaMotherPart",16);
      else if(SigmaMotherPdg==3222) FillHistogram("fHistAntiSigmaMotherPart",17);
      else if(SigmaMotherPdg==3112) FillHistogram("fHistAntiSigmaMotherPart",18);
      else if(SigmaMotherPdg==-3222) FillHistogram("fHistAntiSigmaMotherPart",19);
      else if(SigmaMotherPdg==-3112) FillHistogram("fHistAntiSigmaMotherPart",20);
      else if(SigmaMotherPdg==3322) FillHistogram("fHistAntiSigmaMotherPart",21);
      else if(SigmaMotherPdg==-3322) FillHistogram("fHistAntiSigmaMotherPart",22);
      else if(SigmaMotherPdg==3312) FillHistogram("fHistAntiSigmaMotherPart",23);
      else if(SigmaMotherPdg==-3312) FillHistogram("fHistAntiSigmaMotherPart",24);
      else FillHistogram("fHistAntiSigmaMotherPart",25);

    } //End of AntiSigmaMinus treatment

    if(MCPartPDGCode == 111)   {FillHistogram("fHistMCCounter",9); FillHistogram("fHistMCPi0Pt",mcPart->Pt());
      AliAODMCParticle* Pi0Mother = NULL;
      if(mcPart->GetMother()!=-1) Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      if(Pi0Mother){
        if(TMath::Abs(Pi0Mother->GetPdgCode())==3222){FillHistogram("fHistMCCounter",10);}
      }
    } //End of MCPartPDGCode == 111
    if(MCPartPDGCode == 22)    {
      
      FillHistogram("fHistMCCounter",11);
      FillHistogram("fHistMCPhotonPt",mcPart->Pt());

      if(mcPart->GetNDaughters()!=2) continue;  //Check if Photon has 2 Daughters (i.e. if it was converted to e+ e-)                    
      AliAODMCParticle* FirstDaughtParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetDaughterFirst())));
      AliAODMCParticle* LastDaughtParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetDaughterLast())));
      if(!FirstDaughtParticle || !LastDaughtParticle) continue;

      //Check if Photon Daughters are e+ and e-
      if(TMath::Abs(FirstDaughtParticle->GetPdgCode())!=11) continue; 
      if(TMath::Abs(LastDaughtParticle->GetPdgCode())!=11) continue; 
        
      //Check conversion vertex (aka production vertex of e+/e-)
      TVector3 prdvtx1(FirstDaughtParticle->Xv(),FirstDaughtParticle->Yv(),FirstDaughtParticle->Zv());   
      TVector3 prdvtx2(FirstDaughtParticle->Xv(),FirstDaughtParticle->Yv(),FirstDaughtParticle->Zv());   
      TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;

      FillHistogram("fHistMCConvRadius",prdvtx.Perp());
      FillHistogram("fHistMCConvRadiusvspt",prdvtx.Perp(),mcPart->Pt());
      
      if(prdvtx.Perp()<180){FillHistogram("fHistMCCounter",12); FillHistogram("fHistMCConvPhotonPt",mcPart->Pt());}

      AliAODMCParticle* Pi0Part = NULL;
      AliAODMCParticle* SigmaPart = NULL;
      if(mcPart->GetMother()!=-1){
        Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      }
      if(Pi0Part){
        if(Pi0Part->GetPdgCode()==111&&Pi0Part->GetMother()!=-1) SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
      }
      if(SigmaPart){
        if(TMath::Abs(SigmaPart->GetPdgCode())==3222){
          FillHistogram("fHistMCSigmaPhotonPt",mcPart->Pt());
          if(prdvtx.Perp()<180) FillHistogram("fHistMCSigmaConvPhotonPt",mcPart->Pt());
        }
      }
    } //End of MCPartPDGCode == 22

  } //End of MC Particle loop 

return;

}//End of ProcessMCParticles()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillV0PhotonArray() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear V0 Photon Array and reset counter
  fConvPhotonArray.clear();
  Int_t countPhotons = 0;

  Int_t nV0 = aodEvent->GetNumberOfV0s(); //Number of V0s in the event
  if(nV0 == 0) return;      //Return if there is no V0 to be processed

  Int_t non = 0, noff = 0;

  for (Int_t iV0=0; iV0<nV0; iV0++) {    //Loop over V0s in the event

    //Initialisation of the local Bools
    Bool_t   isElectronTPC = kFALSE;
    Bool_t   isPositronTPC = kFALSE;
    Bool_t   isPhotonTPC   = kFALSE;
    Bool_t   isRealV0       = kFALSE;
    Bool_t   isReallyPhoton = kFALSE;
    Bool_t   isPhotonfromSigma = kFALSE;

    TVector3 vecN, vecP, vecM;                             //Momentum Vectors for V0 tracks
    TLorentzVector electron, positron, photon;             //Lorentzvectors for invariant mass calculation

    // Daughter Track parameters for KF and ExtTrckPar initialization
    Double_t trackxyz[3];
    Double_t trackpxpypz[3];
    Double_t trackparams[6];
    Double_t covMatrix[21];
      
    AliAODv0* aodV0 = (AliAODv0*)aodEvent->GetV0(iV0);        //Get V0 object
    if(!aodV0){
      AliWarning("ERROR: Could not retrieve AOD V0 Object!");
      continue;
    }

    // Check basic V0 properties: 2 Daughters, opposite charge, total charge = 0
    if(aodV0->GetNDaughters() != 2)                    continue;
    if(aodV0->GetNProngs() != 2)                       continue;
    if(aodV0->GetCharge() != 0)                        continue;
    if(aodV0->ChargeProng(0) == aodV0->ChargeProng(1)) continue;

    // Get daughter tracks      
    AliAODTrack* trackN = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
    AliAODTrack* trackP = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
    if (trackN->GetSign() == trackP->GetSign()) continue;
    if (trackN->Charge() > 0){ //Check correct charge 
    trackN = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
    trackP = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
    }

    if(!trackP || !trackN) {
      AliWarning("ERROR: Could not retrieve both AOD daughter tracks of V0!");
      continue;
    }

    //If V0 is not On-fly, check the map if there is an equivalent On-fly V0
    if(!aodV0->GetOnFlyStatus()){ 
      noff++;
      FillHistogram("fHistV0Statistics",2);
      FillHistogram("fHistV0StatisticsMC",2);
      FillHistogram("fHistV0StatisticsSigmaMC",2);

      //Check if the tracks are used by the on-the-fly finder
      Int_t nFound = fOnFlyVector.size();
      Bool_t isused = kFALSE;
      for(Int_t iID = 0; iID<nFound; iID++){
        if(trackN->GetID()==fFinderVector[iID]) isused = kTRUE; 
        if(trackP->GetID()==fFinderVector[iID]) isused = kTRUE; 
      }
      if(isused) continue;

      FillHistogram("fHistV0Statistics",3);
      FillHistogram("fHistV0StatisticsMC",3);
      FillHistogram("fHistV0StatisticsSigmaMC",3);    
    }
    else{
      non++; 
      FillHistogram("fHistV0Statistics",1);
      FillHistogram("fHistV0StatisticsMC",1);
      FillHistogram("fHistV0StatisticsSigmaMC",1);
    }

    // Check track quality
    Int_t nTPCClustNeg = trackN->GetTPCNcls();
    Int_t nTPCClustPos = trackP->GetTPCNcls();
    Int_t nITSClustNeg = trackN->GetITSNcls();
    Int_t nITSClustPos = trackP->GetITSNcls();
  
    // Daughter track PID using TPC
    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kElectron);
    Double_t nSigmaTPCpositron = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kElectron);      
    if(TMath::Abs(nSigmaTPCelectron) < fMaxNsigDaughtTPC) isElectronTPC = kTRUE;    
    if(TMath::Abs(nSigmaTPCpositron) < fMaxNsigDaughtTPC) isPositronTPC = kTRUE; 
    if(isElectronTPC && isPositronTPC) isPhotonTPC = kTRUE;

    // Get topological values
    Double_t dcaV0Daughters  = TMath::Abs(aodV0->DcaV0Daughters());
    Double_t dcaPosToPrimVtx = TMath::Abs(aodV0->DcaPosToPrimVertex());
    Double_t dcaNegToPrimVtx = TMath::Abs(aodV0->DcaNegToPrimVertex());
    Double_t cosPointAngle   = aodV0->CosPointingAngle(primaryVtxPos);
    Double_t vtxPosV0[3];
    vtxPosV0[0]=aodV0->DecayVertexV0X();
    vtxPosV0[1]=aodV0->DecayVertexV0Y();
    vtxPosV0[2]=aodV0->DecayVertexV0Z();
    Double_t Vtxradius=TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0]+vtxPosV0[1]*vtxPosV0[1]);

    //Calculating DCA of Photon to PV 
    TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);                            //Prim. Vertex
    TVector3 CV(aodV0->DecayVertexV0X(),aodV0->DecayVertexV0Y(),aodV0->DecayVertexV0Z()); //Conv. Vertex
    TVector3 p(aodV0->Px(),aodV0->Py(),aodV0->Pz());                   //Momentum vectors of the photons
    Double_t DCAPV = (p.Cross((CV-PV))).Mag()/p.Mag();                            //DCA to PV of Photons

    //Get reconstructed cartesian momentum
    vecN.SetXYZ(aodV0->MomNegX(),aodV0->MomNegY(),aodV0->MomNegZ()); //negative daughter
    vecP.SetXYZ(aodV0->MomPosX(),aodV0->MomPosY(),aodV0->MomPosZ()); //positive daughter 
    vecM.SetXYZ(aodV0->MomV0X(),aodV0->MomV0Y(),aodV0->MomV0Z());    //mother 

    //Custom Armenteros Podolanski calculation since V0 member functions are not reliable!
    Double_t pLNeg = vecN.Dot(vecM)/vecM.Mag(); //Momentum longitudinal
    Double_t pLPos = vecP.Dot(vecM)/vecM.Mag(); //to V0 momentum
    Double_t alpha = (pLPos - pLNeg)/(pLPos + pLNeg);
    Double_t qt    = vecN.Perp(vecM);  

    // Get kinematic values
    Double_t ptV0   = aodV0->Pt();
    Double_t pPos   = trackP->P();
    Double_t pNeg   = trackN->P();
    Double_t thetaPos = trackP->Theta();
    Double_t thetaNeg = trackN->Theta();
    Double_t totangle = aodV0->OpenAngleV0();

    //*******************AOD V0 MC treatment************************//

    // AOD MC treatment
    if(isMonteCarlo){     
      AliAODMCParticle* V0Part = NULL;
      AliAODMCParticle* NPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(trackN->GetLabel())));
      AliAODMCParticle* PPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(trackP->GetLabel())));
      if(NPart&&PPart){
        if(NPart->GetMother()==PPart->GetMother()&&NPart->GetMother()!=-1){
          V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NPart->GetMother()));    
          if(V0Part){if(V0Part->GetPdgCode()==22){ 
            isReallyPhoton = kTRUE;
                  
            AliAODMCParticle* V0Mother = NULL;
            if(V0Part->GetMother()!=-1) V0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part->GetMother())));
            if(V0Mother){
              if(TMath::Abs(V0Mother->GetPdgCode())==3222) FillHistogram("fHistSigmaCounter",5);
              if(V0Mother->GetPdgCode()==111){
                      
              AliAODMCParticle* Pi0Mother = NULL;
              if(V0Mother->GetMother()!=-1) Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Mother->GetMother())));
              if(TMath::Abs(Pi0Mother?Pi0Mother->GetPdgCode():0)==3222) isPhotonfromSigma = kTRUE;
            
            }}//Mother of Photon exists and is a Pi0
          }}//Mother exists and is a Photon
        }//Both Tracks have a common Mother
      }//Both Tracks have matched MC Particle
    }//End of isMonteCarlo

    //***************End of AOD V0 MC treatment************************//

    // Reconstruct photon with TLorentzVector
    electron.SetXYZM(vecN(0),vecN(1),vecN(2),cElectronMass);
    positron.SetXYZM(vecP(0),vecP(1),vecP(2),cElectronMass);
    photon=electron+positron;

    // Calculate photon invariant mass with TL
    Double_t photonmass = photon.M();

    // Angle calculation
    Double_t deltatheta = thetaPos - thetaNeg;

    FillHistogram("fHistV0DaughtEtaPhi",trackN->Eta(),trackN->Phi());
    FillHistogram("fHistV0DaughtEtaPhi",trackP->Eta(),trackP->Phi());
    if(isReallyPhoton) FillHistogram("fHistV0DaughtEtaPhiMC",trackN->Eta(),trackN->Phi());
    if(isReallyPhoton) FillHistogram("fHistV0DaughtEtaPhiMC",trackP->Eta(),trackP->Phi());

    // Acceptance Cut
    if(TMath::Abs(trackN->Eta()) > fMaxDaughtEta) continue; 
    if(TMath::Abs(trackP->Eta()) > fMaxDaughtEta) continue;

    FillHistogram("fHistV0Statistics",4);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",4);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",4);

    // Fill QA Histograms
    FillHistogram("fHistV0DaughtChi2",trackN->GetTPCchi2());
    FillHistogram("fHistV0DaughtChi2",trackP->GetTPCchi2());
    FillHistogram("fHistV0DaughtTPCClust",nTPCClustPos);
    FillHistogram("fHistV0DaughtTPCClust",nTPCClustNeg);
    FillHistogram("fHistV0DaughtITSClust",nITSClustPos);
    FillHistogram("fHistV0DaughtITSClust",nITSClustNeg);
    FillHistogram("fHistV0DaughtpvsNSigmaTPC",pNeg,nSigmaTPCelectron);
    FillHistogram("fHistV0DaughtpvsNSigmaTPC",pPos,nSigmaTPCpositron);
    FillHistogram("fHistV0DaughtDCAtoPV",dcaPosToPrimVtx);
    FillHistogram("fHistV0DaughtDCAtoPV",dcaNegToPrimVtx);
    FillHistogram("fHistV0DaughtDCA",dcaV0Daughters);
    FillHistogram("fHistV0CPA",cosPointAngle);
    FillHistogram("fHistV0Radius",Vtxradius);
    FillHistogram("fHistV0Position2D",vtxPosV0[0],vtxPosV0[1]);
    FillHistogram("fHistV0PhotonDCAPV",DCAPV);
    FillHistogram("fHistV0ArmPod",alpha,qt);
    FillHistogram("fHistV0OpenAngle",totangle);
    FillHistogram("fHistV0DeltaTheta",deltatheta);
    FillHistogram("fHistV0InvMass",photonmass);

    if(isReallyPhoton){
      FillHistogram("fHistV0DaughtChi2MC",trackN->GetTPCchi2());
      FillHistogram("fHistV0DaughtChi2MC",trackP->GetTPCchi2());
      FillHistogram("fHistV0DaughtDCAtoPVMC",dcaPosToPrimVtx);
      FillHistogram("fHistV0DaughtDCAtoPVMC",dcaNegToPrimVtx);
      FillHistogram("fHistV0DaughtTPCClustMC",nTPCClustPos);
      FillHistogram("fHistV0DaughtTPCClustMC",nTPCClustNeg);
      FillHistogram("fHistV0DaughtITSClustMC",nITSClustPos);
      FillHistogram("fHistV0DaughtITSClustMC",nITSClustNeg);
      FillHistogram("fHistV0DaughtDCAMC",dcaV0Daughters);
      FillHistogram("fHistV0CPAMC",cosPointAngle);
      FillHistogram("fHistV0RadiusMC",Vtxradius);
      FillHistogram("fHistV0Position2DMC",vtxPosV0[0],vtxPosV0[1]);
      FillHistogram("fHistV0PhotonDCAPVMC",DCAPV);
      FillHistogram("fHistV0ArmPodMC",alpha,qt);
      FillHistogram("fHistV0OpenAngleMC",totangle);
      FillHistogram("fHistV0DeltaThetaMC",deltatheta);
      FillHistogram("fHistV0InvMassMC",photonmass);
      FillHistogram("fHistV0ptMC",ptV0);
    }

    if(isPhotonfromSigma){
      FillHistogram("fHistV0CPAMCSigma",cosPointAngle);
      FillHistogram("fHistV0SigmaptMC",ptV0);
      FillHistogram("fHistV0PhotonDCAPVMCSigma",DCAPV);
    }

    // Check Track quality and reject poor qualty tracks
    if(nTPCClustNeg < fMinTPCClustDaught) continue;
    if(nTPCClustPos < fMinTPCClustDaught) continue;

    FillHistogram("fHistV0Statistics",5);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",5);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",5);

    // Armenteros-Podolanski Cuts
    if(TMath::Abs(alpha) > fMaxalpha) continue;

    FillHistogram("fHistV0Statistics",6);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",6);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",6);

    if(TMath::Abs(qt) > fMaxqt) continue;

    FillHistogram("fHistV0Statistics",7);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",7);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",7);

    // Angle Cut
    if(TMath::Abs(totangle) > fMaxopenangle) continue;

    FillHistogram("fHistV0Statistics",8);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",8);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",8);

    if(TMath::Abs(deltatheta) > fMaxdeltatheta) continue;

    FillHistogram("fHistV0Statistics",9);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",9);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",9);

    // CPA Cut
    if (cosPointAngle < fMinV0CPA) continue;

    FillHistogram("fHistV0Statistics",10);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",10);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",10);

    // PID Cut
    if(!isPhotonTPC) continue;    

    FillHistogram("fHistV0Statistics",11);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",11);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",11);

    //Radius Cut
    if(Vtxradius < fMinV0Radius) continue;

    FillHistogram("fHistV0Statistics",12);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",12);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",12);

    if(Vtxradius > fMaxV0Radius) continue;

    FillHistogram("fHistV0Statistics",13);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",13);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",13);

    // Inv. Mass Cut
    if(photonmass > fMaxphotonmass) continue;
    FillHistogram("fHistV0Statistics",14);
    if(isReallyPhoton) FillHistogram("fHistV0StatisticsMC",14);
    if(isPhotonfromSigma) FillHistogram("fHistV0StatisticsSigmaMC",14);

    // Histograms for Efficiency and Purity calculation
    if(isReallyPhoton) FillHistogram("fHistV0ptwCutsMC",ptV0);
    if(isPhotonfromSigma) FillHistogram("fHistV0SigmaptwCutsMC",ptV0);      
    FillHistogram("fHistV0ptwCuts",ptV0); 

    // Store Photon candidates after selection
    fV0ParticleIDArray.push_back(trackN->GetID());
    fV0ParticleIDArray.push_back(trackP->GetID());
        
    fConvPhotonArray.push_back(iV0);
    countPhotons++;

  }//End of V0 Loop

FillHistogram("fHistV0OnflyvsOffline",non,noff);

return;

}//End of FillV0PhotonArray()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillCaloClusterArray() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Calo Photon Array and reset counter
  fCaloPhotonArray.clear();
  Int_t countPhotons = 0;

  Int_t nCluster = aodEvent->GetNumberOfCaloClusters(); //Number of Clusters in the event
  FillHistogram("fHistNClusters",nCluster);
  if(nCluster == 0) return;               //Return if there is no Cluster to be processed

  //Loop for Proton Selection
  for(Int_t iCluster=0; iCluster < nCluster; iCluster++) { //Loop over Clusters in the event

    AliAODCaloCluster *aodCluster = dynamic_cast<AliAODCaloCluster*>(aodEvent->GetCaloCluster(iCluster));
    if(!aodCluster) {
      AliWarning("No AOD Cluster!");
      continue;
    }

    Float_t x[3];
    aodCluster->GetPosition(x);

    Float_t E = aodCluster->E();
    Int_t NLabels = aodCluster->GetNLabels();
    Float_t Dispersion = TMath::Abs(aodCluster->GetDispersion());
    Float_t M20 = TMath::Abs(aodCluster->GetM20());
    Float_t M02 = TMath::Abs(aodCluster->GetM02());
    Int_t NTracksMatched = aodCluster->GetNTracksMatched();
    Float_t TrackDx = TMath::Abs(aodCluster->GetTrackDx());
    Float_t TrackDz = TMath::Abs(aodCluster->GetTrackDz());
    Double_t TOF = 1000000000*aodCluster->GetTOF();
    Int_t NCells = aodCluster->GetNCells();
    Float_t DisttoBC = aodCluster->GetDistanceToBadChannel();

    Double_t FlightDist = TMath::Sqrt((x[0]-primaryVtxPosX)*(x[0]-primaryVtxPosX)+(x[1]-primaryVtxPosY)*(x[1]-primaryVtxPosY)+(x[2]-primaryVtxPosZ)*(x[2]-primaryVtxPosZ));
    Double_t vTOF = 0;
    if(TOF!=0) vTOF = FlightDist/(1000*TOF);
    Float_t beta = vTOF/c;

    Bool_t isEmcal = aodCluster->IsEMCAL();
    Bool_t isPhos = aodCluster->IsPHOS();
    Int_t type = -1; 
    if(isPhos) type = 0;
    else if(isEmcal) type = 1;
    else type = 2;

    if(TrackDx==0) TrackDx=-999;
    if(TrackDz==0) TrackDz=-999;

    //Check closest track if its EMCAL
    if(fMapClusterstoTracks&&isEmcal){
      //Radial position of EMCAL Cluster
      Double_t Rcalo = TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
      //Point of track at radius Rcalo
      Double_t pnt[3];
      //Absolute distance of closest track
      Double_t minD = 999;
      //Loop over Tracks
      for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {
        AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
        if(!aodTrack) continue;
        //Propagate track to radius Rcalo
        if(!aodTrack->GetXYZatR(Rcalo,Bz,pnt,0)) continue; //Continue if propagation fails
        //Calculate distances
        Double_t dx = TMath::Abs(pnt[0]-x[0]);
        Double_t dy = TMath::Abs(pnt[1]-x[1]);
        Double_t dz = TMath::Abs(pnt[2]-x[2]);
        Double_t d = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
        if(d>minD) continue; //Continue if there is already a closer track
        minD=d;
        TrackDz=dz;
        TrackDx=TMath::Sqrt(dx*dx+dy*dy);
        //Overwrite track distances for later use
        aodCluster->SetTrackDistance(TrackDx,TrackDz);
      }
    }//End of track mapping

    //Absolute distance of closest track 
    Float_t TrackDy = -999;
    if(TrackDx!=-999) TrackDy = TMath::Sqrt(TrackDx*TrackDx+TrackDz*TrackDz);

    //Check MC Truth
    Bool_t isPhoton = kFALSE;
    Bool_t isfromSigma = kFALSE;

    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodCluster->GetLabel())));
      AliAODMCParticle* Pi0Part = NULL;
      AliAODMCParticle* SigmaPart = NULL;
      if(mcPart){
        if(mcPart->GetPdgCode()==22) isPhoton = kTRUE;        
        if(mcPart->GetMother()!=-1&&isPhoton){
          Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
        }
        if(Pi0Part){
          if(Pi0Part->GetPdgCode()==111&&Pi0Part->GetMother()!=-1) SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
        }
        if(SigmaPart){
          if(TMath::Abs(SigmaPart->GetPdgCode())==3222) isfromSigma = kTRUE;
        }
      }//MC Particle exists
    }//MC treatment

    if(isPhos){
      FillHistogram("fHistClusterType",type);
      FillHistogram("fHistClusterTypeMC",type);
      FillHistogram("fHistClusterTypeMCSig",type);
      FillHistogram("fHistPHOSDisttoBC",DisttoBC);
      FillHistogram("fHistPHOSM02",M02);
      FillHistogram("fHistPHOSM20",M20);
      FillHistogram("fHistPHOSTOF",TOF);
      FillHistogram("fHistPHOSBeta",beta);
      FillHistogram("fHistPHOSE",E);
      FillHistogram("fHistPHOSNTracks",NTracksMatched);
      FillHistogram("fHistPHOSNCells",NCells);
      FillHistogram("fHistPHOSDx",TrackDx);
      FillHistogram("fHistPHOSDz",TrackDz);
      FillHistogram("fHistPHOSDy",TrackDy);
      FillHistogram("fHistPHOSDisp",Dispersion);
      if(isPhoton){
        FillHistogram("fHistPHOSMCDisttoBC",DisttoBC);
        FillHistogram("fHistPHOSMCM02",M02);
        FillHistogram("fHistPHOSMCM20",M20);
        FillHistogram("fHistPHOSMCTOF",TOF);
        FillHistogram("fHistPHOSMCBeta",beta);
        FillHistogram("fHistPHOSMCE",E);
        FillHistogram("fHistPHOSMCNTracks",NTracksMatched);
        FillHistogram("fHistPHOSMCNCells",NCells);
        FillHistogram("fHistPHOSMCDx",TrackDx);
        FillHistogram("fHistPHOSMCDz",TrackDz);
        FillHistogram("fHistPHOSMCDy",TrackDy);
        FillHistogram("fHistPHOSMCDisp",Dispersion);
      }
      if(isfromSigma){
        FillHistogram("fHistPHOSMCSigDisttoBC",DisttoBC);
        FillHistogram("fHistPHOSMCSigM02",M02);
        FillHistogram("fHistPHOSMCSigM20",M20);
        FillHistogram("fHistPHOSMCSigTOF",TOF);
        FillHistogram("fHistPHOSMCSigBeta",beta);
        FillHistogram("fHistPHOSMCSigE",E);
        FillHistogram("fHistPHOSMCSigNTracks",NTracksMatched);
        FillHistogram("fHistPHOSMCSigNCells",NCells);
        FillHistogram("fHistPHOSMCSigDx",TrackDx);
        FillHistogram("fHistPHOSMCSigDz",TrackDz);
        FillHistogram("fHistPHOSMCSigDy",TrackDy);
        FillHistogram("fHistPHOSMCSigDisp",Dispersion);
      }
    }
    if(isEmcal){
      FillHistogram("fHistEMCALDisttoBC",DisttoBC);
      FillHistogram("fHistEMCALM02",M02);
      FillHistogram("fHistEMCALM20",M20);
      FillHistogram("fHistEMCALTOF",TOF);
      FillHistogram("fHistEMCALBeta",beta);
      FillHistogram("fHistEMCALE",E);
      FillHistogram("fHistEMCALNTracks",NTracksMatched);
      FillHistogram("fHistEMCALNCells",NCells);
      FillHistogram("fHistEMCALDx",TrackDx);
      FillHistogram("fHistEMCALDz",TrackDz);
      FillHistogram("fHistEMCALDy",TrackDy);
      FillHistogram("fHistEMCALDisp",Dispersion);
      if(isPhoton){
        FillHistogram("fHistEMCALMCDisttoBC",DisttoBC);
        FillHistogram("fHistEMCALMCM02",M02);
        FillHistogram("fHistEMCALMCM20",M20);
        FillHistogram("fHistEMCALMCTOF",TOF);
        FillHistogram("fHistEMCALMCBeta",beta);
        FillHistogram("fHistEMCALMCE",E);
        FillHistogram("fHistEMCALMCNTracks",NTracksMatched);
        FillHistogram("fHistEMCALMCNCells",NCells);
        FillHistogram("fHistEMCALMCDx",TrackDx);
        FillHistogram("fHistEMCALMCDz",TrackDz);
        FillHistogram("fHistEMCALMCDy",TrackDy);
        FillHistogram("fHistEMCALMCDisp",Dispersion);
      }
      if(isfromSigma){
        FillHistogram("fHistEMCALMCSigDisttoBC",DisttoBC);
        FillHistogram("fHistEMCALMCSigM02",M02);
        FillHistogram("fHistEMCALMCSigM20",M20);
        FillHistogram("fHistEMCALMCSigTOF",TOF);
        FillHistogram("fHistEMCALMCSigBeta",beta);
        FillHistogram("fHistEMCALMCSigE",E);
        FillHistogram("fHistEMCALMCSigNTracks",NTracksMatched);
        FillHistogram("fHistEMCALMCSigNCells",NCells);
        FillHistogram("fHistEMCALMCSigDx",TrackDx);
        FillHistogram("fHistEMCALMCSigDz",TrackDz);
        FillHistogram("fHistEMCALMCSigDy",TrackDy);
        FillHistogram("fHistEMCALMCSigDisp",Dispersion);
      }
    }

    if(type<0||type>1) continue;

    FillHistogram("fHistClusterStatistics",1);
    if(isPhoton) FillHistogram("fHistClusterStatisticsMC",1);
    if(isfromSigma) FillHistogram("fHistClusterStatisticsMCSig",1);

    // PHOS Cut
    if(!isPhos&&fRequirePHOS) continue;

    FillHistogram("fHistClusterStatistics",2);
    if(isPhoton) FillHistogram("fHistClusterStatisticsMC",2);
    if(isfromSigma) FillHistogram("fHistClusterStatisticsMCSig",2);

    // Beta Cut
    if(beta<fMinClusterBeta) continue;

    FillHistogram("fHistClusterStatistics",3);
    if(isPhoton) FillHistogram("fHistClusterStatisticsMC",3);
    if(isfromSigma) FillHistogram("fHistClusterStatisticsMCSig",3);

    // Dy Cut
    if(TrackDy<fMinClusterDy) continue;

    FillHistogram("fHistClusterStatistics",4);
    if(isPhoton) FillHistogram("fHistClusterStatisticsMC",4);
    if(isfromSigma) FillHistogram("fHistClusterStatisticsMCSig",4);

    // M02 Cut
    if(M02>fMaxClusterM02) continue;

    FillHistogram("fHistClusterStatistics",5);
    if(isPhoton) FillHistogram("fHistClusterStatisticsMC",5);
    if(isfromSigma) FillHistogram("fHistClusterStatisticsMCSig",5);

    fCaloPhotonArray.push_back(iCluster);
    countPhotons++;

  }//End of Cluster Loop

  FillHistogram("fHistNCaloPhotons",countPhotons);

  return;

}//End of FillCaloClusterArray()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticles() {

  KFParticle::SetField(Bz);

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TLorentzVector trackPhoton1, trackPhoton2, trackPi0, trackProton, trackSigmaplus;
  KFParticle KFElectron1, KFElectron2, KFPositron1, KFPositron2, KFProton; 

  const Int_t nConvPhoton = fConvPhotonArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nConvPhoton-1; i++) {

    AliAODv0 *v0_1 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(i));
    if(!v0_1) continue;

    for(Int_t j=i+1; j<nConvPhoton; j++) {

      AliAODv0 *v0_2 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(j));
      if(!v0_2) continue;

      // Get daughter tracks      
      AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(0));
      AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(1));
      AliAODTrack* track3 = dynamic_cast<AliAODTrack*>(v0_2->GetDaughter(0));
      AliAODTrack* track4 = dynamic_cast<AliAODTrack*>(v0_2->GetDaughter(1));

      if(!track1 || !track2 || !track3 || !track4) {
        AliWarning("ERROR: Could not retrieve all AOD tracks in Pi0 reconstruction!");
        continue;
      }

      // AOD MC treatment
      Double_t MCPi0DCAPV; //MC Pi0 Decay Vertex;
      Double_t MCPi0decayX; //MC Pi0 Decay Vertex;
      Double_t MCPi0decayY; //MC Pi0 Decay Vertex;
      Double_t MCPi0decayZ; //MC Pi0 Decay Vertex;
      Bool_t isReallyPi0 = kFALSE; 
      Bool_t isSameGamma = kFALSE;
      Bool_t isReallyPi0fromSigma = kFALSE;
      Bool_t isReallyPi0fromDelta = kFALSE;
      Bool_t isPrimary = kFALSE;
      Bool_t isDalitzDecay = kFALSE;
      Int_t Pi0MotherLabel=-1;
      TLorentzVector MCSigmaMom; 

      if(isMonteCarlo){
        AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
        AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
        AliAODMCParticle* V02Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track3->GetLabel())));
        AliAODMCParticle* V02Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track4->GetLabel())));
        if(V01Daught1&&V01Daught2&&V02Daught1&&V02Daught2){
          if(TMath::Abs(V01Daught1->GetPdgCode())==11&&TMath::Abs(V01Daught2->GetPdgCode())==11&&TMath::Abs(V02Daught1->GetPdgCode())==11&&TMath::Abs(V02Daught2->GetPdgCode())==11){
            if(V01Daught1->GetMother()!=-1&&V02Daught1->GetMother()!=-1&&V01Daught1->GetMother()==V01Daught2->GetMother()&&V02Daught1->GetMother()==V02Daught2->GetMother()){
              AliAODMCParticle* V0Part1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
              AliAODMCParticle* V0Part2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V02Daught1->GetMother())));
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22&&TMath::Abs(V0Part1->GetLabel())==TMath::Abs(V0Part2->GetLabel())) isSameGamma = kTRUE;}
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22&&TMath::Abs(V0Part1->GetLabel())!=TMath::Abs(V0Part2->GetLabel())){
                if(V0Part1->GetMother()!=-1&&V0Part1->GetMother()==V0Part2->GetMother()){

                  AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                  if(Pi0Part){if(Pi0Part->GetPdgCode()==111){
                  
                    isReallyPi0=kTRUE;

                    TVector3 prdvtx1(V0Part1->Xv(),V0Part1->Yv(),V0Part1->Zv());   
                    TVector3 prdvtx2(V0Part2->Xv(),V0Part2->Yv(),V0Part2->Zv());   
                    TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
                    MCPi0DCAPV = prdvtx.Perp();

                    MCPi0decayX = prdvtx.X();
                    MCPi0decayY = prdvtx.Y();
                    MCPi0decayZ = prdvtx.Z();

                    AliAODMCParticle* Pi0Mother = NULL;
                    if(Pi0Part->GetMother()!=-1){
                      Pi0MotherLabel=Pi0Part->GetMother(); 
                      Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                    }
                    if(Pi0Mother){
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) {
                        isReallyPi0fromSigma = kTRUE; 
                        if(Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE; 
                        MCSigmaMom.SetXYZM(Pi0Mother->Px(),Pi0Mother->Py(),Pi0Mother->Pz(),cSigmaMass);
                      }
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==2214) isReallyPi0fromDelta = kTRUE;
                    }

                  }}//Photon Mother exists and is Pi0
                }//Photon and Single electron have common mother
              }}//MC Photon exists
              //Consider dalitz decays
              if(V0Part1&&V0Part2){if((V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==111)||(V0Part1->GetPdgCode()==111&&V0Part2->GetPdgCode()==22)){
                if(TMath::Abs(V0Part1->GetMother())==TMath::Abs(V0Part2->GetLabel())||TMath::Abs(V0Part2->GetMother())==TMath::Abs(V0Part1->GetLabel())){
                  isReallyPi0=kTRUE;
                  isDalitzDecay = kTRUE;

                  TVector3 prdvtx;
                  if(V0Part1->GetPdgCode()==22) prdvtx.SetXYZ(V0Part1->Xv(),V0Part1->Yv(),V0Part1->Zv());   
                  if(V0Part2->GetPdgCode()==22) prdvtx.SetXYZ(V0Part2->Xv(),V0Part2->Yv(),V0Part2->Zv());   
                  MCPi0DCAPV = prdvtx.Perp();
                  MCPi0decayX = prdvtx.X();
                  MCPi0decayY = prdvtx.Y();
                  MCPi0decayZ = prdvtx.Z();

                  AliAODMCParticle* Pi0Part = NULL;
                  if(V0Part1->GetPdgCode()==111) Pi0Part = V0Part1;
                  else Pi0Part = V0Part2;
                  AliAODMCParticle* Pi0Mother = NULL;

                  if(Pi0Part->GetMother()!=-1){
                    Pi0MotherLabel=Pi0Part->GetMother(); 
                    Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                  }

                  if(Pi0Mother){
                    if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) {
                      isReallyPi0fromSigma = kTRUE; 
                      if(Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE; 
                      MCSigmaMom.SetXYZM(Pi0Mother->Px(),Pi0Mother->Py(),Pi0Mother->Pz(),cSigmaMass);
                    }
                    if(TMath::Abs(Pi0Mother->GetPdgCode())==2214) isReallyPi0fromDelta = kTRUE;
                  }

                } //V0s belong together
              }} //V0s exist and one is Photon and one is from dalitz decay
            }//Electrons have common mother
          }//Daughters are all electrons
        }//V0 Daughters exist
      }//End of isMonteCarlo

      trackPhoton1.SetXYZM(v0_1->Px(),v0_1->Py(),v0_1->Pz(),0);
      trackPhoton2.SetXYZM(v0_2->Px(),v0_2->Py(),v0_2->Pz(),0);
      trackPi0 = trackPhoton1 + trackPhoton2;
      Double_t pi0mass = trackPi0.M();

      //KF Pi0 calculations
      // Set up KFParticle
      trackparams[0] = v0_1->DecayVertexV0X();
      trackparams[1] = v0_1->DecayVertexV0Y();
      trackparams[2] = v0_1->DecayVertexV0Z();
      trackparams[3] = v0_1->MomPosX();
      trackparams[4] = v0_1->MomPosY();
      trackparams[5] = v0_1->MomPosZ();
      if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron1.Create(trackparams,covMatrix,1,cElectronMass);

      // Repeat for all other particles
      trackparams[0] = v0_1->DecayVertexV0X();
      trackparams[1] = v0_1->DecayVertexV0Y();
      trackparams[2] = v0_1->DecayVertexV0Z();
      trackparams[3] = v0_1->MomNegX();
      trackparams[4] = v0_1->MomNegY();
      trackparams[5] = v0_1->MomNegZ();
      if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron1.Create(trackparams,covMatrix,-1,cElectronMass);

      trackparams[0] = v0_2->DecayVertexV0X();
      trackparams[1] = v0_2->DecayVertexV0Y();
      trackparams[2] = v0_2->DecayVertexV0Z();
      trackparams[3] = v0_2->MomPosX();
      trackparams[4] = v0_2->MomPosY();
      trackparams[5] = v0_2->MomPosZ();
      if(track3->Charge()>0) track3->GetCovarianceXYZPxPyPz(covMatrix);
      else track4->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron2.Create(trackparams,covMatrix,1,cElectronMass);

      trackparams[0] = v0_2->DecayVertexV0X();
      trackparams[1] = v0_2->DecayVertexV0Y();
      trackparams[2] = v0_2->DecayVertexV0Z();
      trackparams[3] = v0_2->MomNegX();
      trackparams[4] = v0_2->MomNegY();
      trackparams[5] = v0_2->MomNegZ();
      if(track3->Charge()<0) track3->GetCovarianceXYZPxPyPz(covMatrix);
      else track4->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron2.Create(trackparams,covMatrix,-1,cElectronMass);

      KFParticleCD KFPhoton1CD; //Check Daughters to avoid floating point exceptions. See .h-file
      KFPhoton1CD.AddDaughter(KFElectron1);
      if(!KFPhoton1CD.CheckDaughter(KFPositron1)) {FillHistogram("fHistPhotonKFCheck",1); if(isReallyPi0) FillHistogram("fHistPhotonKFCheck",2); continue;}
      else{ FillHistogram("fHistPhotonKFCheck",3); if(isReallyPi0) FillHistogram("fHistPhotonKFCheck",4);}

      KFParticleCD KFPhoton2CD; //Check Daughters to avoid floating point exceptions. See .h-file
      KFPhoton2CD.AddDaughter(KFElectron2);
      if(!KFPhoton2CD.CheckDaughter(KFPositron2)) {FillHistogram("fHistPhotonKFCheck",1); if(isReallyPi0) FillHistogram("fHistPhotonKFCheck",2); continue;}
      else{ FillHistogram("fHistPhotonKFCheck",3); if(isReallyPi0) FillHistogram("fHistPhotonKFCheck",4);}

      //Reconstruct the Photons with two Methods: Standard and special gamma reconstruction
      KFParticle KFPhoton1(KFElectron1,KFPositron1);
      KFParticle KFPhoton2(KFElectron2,KFPositron2);

      //Transport Photons to Conversion Points
      KFPhoton1.TransportToDecayVertex();
      KFPhoton2.TransportToDecayVertex();

      Double_t PhotPhotDCA = TMath::Abs(KFPhoton1.GetDistanceFromParticle(KFPhoton2));

      //Calculating DCA of Photons to PV 
      TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);                            //Prim. Vertex
      TVector3 CV1(v0_1->DecayVertexV0X(),v0_1->DecayVertexV0Y(),v0_1->DecayVertexV0Z()); //Conv. Vertices
      TVector3 CV2(v0_2->DecayVertexV0X(),v0_2->DecayVertexV0Y(),v0_2->DecayVertexV0Z());
      TVector3 p1(v0_1->Px(),v0_1->Py(),v0_1->Pz());                     //Momentum vectors of the photons
      TVector3 p2(v0_2->Px(),v0_2->Py(),v0_2->Pz());   
      Double_t DCAPV1 = (p1.Cross((CV1-PV))).Mag()/p1.Mag(); //DCA to PV of Photons
      Double_t DCAPV2 = (p2.Cross((CV2-PV))).Mag()/p2.Mag(); //using line-point distance equation

  	  FillHistogram("fHistGammaPairInvMass",pi0mass);
  	  FillHistogram("fHistGammaPairDCA",PhotPhotDCA);
  	  if(v0_1->GetOnFlyStatus()&&v0_2->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassOnfly",pi0mass);

      Bool_t hasdiffindices = kFALSE;  
      if(track1->GetID()!=track3->GetID()&&track1->GetID()!=track4->GetID()&&track2->GetID()!=track3->GetID()&&track2->GetID()!=track4->GetID()) hasdiffindices = kTRUE;

      FillHistogram("fHistGammaPairStats",1);
      if(isReallyPi0) FillHistogram("fHistGammaPairStats",2);
      if(isSameGamma) FillHistogram("fHistGammaPairStats",3);
      if(hasdiffindices) FillHistogram("fHistGammaPairStats",4);
      if(isReallyPi0&&hasdiffindices) FillHistogram("fHistGammaPairStats",5);
      if(isSameGamma&&hasdiffindices) FillHistogram("fHistGammaPairStats",6);

      if(hasdiffindices){
  	    FillHistogram("fHistGammaPairInvMass2",pi0mass);
  	    FillHistogram("fHistGammaPairDCA2",PhotPhotDCA);
  	    if(v0_1->GetOnFlyStatus()&&v0_2->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassOnfly2",pi0mass);
      }
      else{
  	    FillHistogram("fHistGammaPairInvMass3",pi0mass);
  	    FillHistogram("fHistGammaPairDCA3",PhotPhotDCA);
  	    if(v0_1->GetOnFlyStatus()&&v0_2->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassOnfly3",pi0mass);
      }

      if(isReallyPi0){
  	    FillHistogram("fHistGammaPairInvMassMC",pi0mass);
  	    FillHistogram("fHistGammaPairDCAMC",PhotPhotDCA);
  	    if(v0_1->GetOnFlyStatus()&&v0_2->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassMCOnfly",pi0mass);
      }
      
      if(!hasdiffindices&&fCleanAutoCorr) continue;  //Discard Photon Pairs that share at least one Track. Is most probably auto-correlation!
      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

      // Save Photon quality
      Int_t nMinTPCClustDaught = track1->GetTPCNcls();
      if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();
      if(track3->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track3->GetTPCNcls();
      if(track4->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track4->GetTPCNcls();

      Int_t nMinITSClustDaught = track1->GetITSNcls();
      if(track2->GetITSNcls()<nMinITSClustDaught) nMinITSClustDaught = track2->GetITSNcls();
      if(track3->GetITSNcls()<nMinITSClustDaught) nMinITSClustDaught = track3->GetITSNcls();
      if(track4->GetITSNcls()<nMinITSClustDaught) nMinITSClustDaught = track4->GetITSNcls();

      Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron));

      Double_t MaxAlpha = TMath::Abs(v0_1->AlphaV0());
      if(TMath::Abs(v0_2->AlphaV0())>MaxAlpha) MaxAlpha = TMath::Abs(v0_2->AlphaV0());

      Double_t MaxQt = TMath::Abs(v0_1->PtArmV0());
      if(TMath::Abs(v0_2->PtArmV0())>MaxQt) MaxQt = TMath::Abs(v0_2->PtArmV0());

      Double_t MaxOpenAngle = v0_1->OpenAngleV0();
      if(v0_2->OpenAngleV0()>MaxOpenAngle) MaxOpenAngle = v0_2->OpenAngleV0();

      TLorentzVector tl1, tl2, tlp;
      tl1.SetXYZM(v0_1->MomNegX(),v0_1->MomNegY(),v0_1->MomNegZ(),cElectronMass);
      tl2.SetXYZM(v0_1->MomPosX(),v0_1->MomPosY(),v0_1->MomPosZ(),cElectronMass);
      tlp=tl1+tl2;

      Double_t Maxphotonmass = tlp.M();

      tl1.SetXYZM(v0_2->MomNegX(),v0_2->MomNegY(),v0_2->MomNegZ(),cElectronMass);
      tl2.SetXYZM(v0_2->MomPosX(),v0_2->MomPosY(),v0_2->MomPosZ(),cElectronMass);
      tlp=tl1+tl2;

      if(tlp.M()>Maxphotonmass) Maxphotonmass = tlp.M();

      Double_t nMaxPhotchi2 = v0_1->Chi2V0();
      if(v0_2->Chi2V0()>nMaxPhotchi2) nMaxPhotchi2 = v0_2->Chi2V0();

      Double_t MaxDaughtEta = TMath::Abs(track1->Eta());
      if(TMath::Abs(track2->Eta())>MaxDaughtEta) MaxDaughtEta = TMath::Abs(track2->Eta());
      if(TMath::Abs(track3->Eta())>MaxDaughtEta) MaxDaughtEta = TMath::Abs(track3->Eta());
      if(TMath::Abs(track4->Eta())>MaxDaughtEta) MaxDaughtEta = TMath::Abs(track4->Eta());

      Double_t MaxDeltaTheta = TMath::Abs(track1->Theta()-track2->Theta());
      if(TMath::Abs(track3->Theta()-track4->Theta())>MaxDeltaTheta) MaxDeltaTheta = TMath::Abs(track3->Theta()-track4->Theta());

      /************************Sigma+ reconstruction************************************/

        for(Int_t k=0; k<nProton; k++) {

          AliAODTrack *prot;
          prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(k));
          if(!prot) continue;

          trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);
          trackSigmaplus = trackPi0 + trackProton;

          Float_t sigmaplusmass = trackSigmaplus.M();
          FillHistogram("fHistInvSigmaMass",sigmaplusmass);
          if(sigmaplusmass>fMaxSigmaMass) continue;   //Limit the mass range to reduce tree size

          prot->GetXYZ(trackxyz);      
          prot->GetPxPyPz(trackpxpypz);
          for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
          for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
          prot->GetCovarianceXYZPxPyPz(covMatrix);
          if(prot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
          else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

          //Reconstruct the Pi0 with the gammas
          KFParticleCD KFPi0CD; //Check Daughters to avoid floating point exceptions. See .h-file
          KFPi0CD.AddDaughter(KFPhoton1);
          if(!KFPi0CD.CheckDaughter(KFPhoton2)) {FillHistogram("fHistPi0KFCheck",1); if(isReallyPi0) FillHistogram("fHistPi0KFCheck",2); continue;}
          else{ FillHistogram("fHistPi0KFCheck",3); if(isReallyPi0) FillHistogram("fHistPi0KFCheck",4);}

          KFParticle KFPi0(KFPhoton1,KFPhoton2);
          KFPi0.TransportToDecayVertex();

          Double_t KFPi0DCAPV = TMath::Sqrt(KFPi0.GetX()*KFPi0.GetX()+KFPi0.GetY()*KFPi0.GetY());  
          if(isReallyPi0){ 
            FillHistogram("fHistPi0VertexMC",MCPi0DCAPV);
            FillHistogram("fHistPi0VertexvsMC",KFPi0DCAPV-MCPi0DCAPV);
          }

          KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
          KFSigmaPlusCD.AddDaughter(KFProton);
          if(!KFSigmaPlusCD.CheckDaughter(KFPi0)) {FillHistogram("fHistSigmaKFCheck",1); if(isReallyPi0) FillHistogram("fHistSigmaKFCheck",2); continue;}
          else{ FillHistogram("fHistSigmaKFCheck",3); if(isReallyPi0) FillHistogram("fHistSigmaKFCheck",4);}

          KFParticle KFSigmaPlus(KFProton,KFPi0);
          KFSigmaPlus.TransportToDecayVertex();
          Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));

          Bool_t isReallySigma = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){if(TMath::Abs(ProtonPart->GetPdgCode())==2212){
              if(TMath::Abs(ProtonPart->GetMother())==TMath::Abs(Pi0MotherLabel)){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",1);}
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          if(isReallySigma){
            FillHistogram("fHistKFSigmaVertexResX",KFSigmaPlus.GetX()-MCPi0decayX);            
            FillHistogram("fHistKFSigmaVertexResY",KFSigmaPlus.GetY()-MCPi0decayY);            
            FillHistogram("fHistKFSigmaVertexResZ",KFSigmaPlus.GetZ()-MCPi0decayZ);              
          }

          //Calculating DCA of Photon and Proton to Decay Vertex 
          TVector3 SV(KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()); //Sec. Vertex
          Double_t DCASV1 = (p1.Cross((CV1-SV))).Mag()/p1.Mag();                 //DCA to SV of Photons
          Double_t DCASV2 = (p1.Cross((CV2-SV))).Mag()/p2.Mag();                 
          Double_t Trackpnt[3];
          Double_t SigmaRadius = TMath::Sqrt((KFSigmaPlus.GetX()*KFSigmaPlus.GetX())+(KFSigmaPlus.GetY()*KFSigmaPlus.GetY()));
          prot->GetXYZatR(SigmaRadius,Bz,Trackpnt,0);
          TVector3 TrackpntVec(Trackpnt[0],Trackpnt[1],Trackpnt[2]);  
          Double_t DCATrack = (TrackpntVec-SV).Mag();
          Double_t DCATrackKF = -999;
          //DCATrackKF = TMath::Abs(KFProton.GetDistanceFromParticle(KFSigmaPlus));

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

          TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
          Float_t Rapidity = sig.Rapidity();
          if(isReallySigma) FillHistogram("fHistMCSigmaY",Rapidity);
          FillHistogram("fHistSigmaY",Rapidity);
          if(TMath::Abs(Rapidity)>fMaxSigmaY) continue;  //Cut on the rapidity

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);

          TLorentzVector trackProtonRot = trackProton;
          TVector3 protonpath(KFSigmaPlus.GetX()-trackxyz[0],KFSigmaPlus.GetY()-trackxyz[1],KFSigmaPlus.GetZ()-trackxyz[2]);
          Double_t propdir = 1; 
          if(TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY))
          <TMath::Sqrt((trackxyz[0]-primaryVtxPosX)*(trackxyz[0]-primaryVtxPosX)+(trackxyz[1]-primaryVtxPosY)*(trackxyz[1]-primaryVtxPosY))) propdir = -1;
          Double_t qprot = 1; if(prot->Charge()<0) qprot=-1;
          Double_t Lprop = protonpath.Mag();
          Double_t Rcurve = trackProton.Vect().Mag()*1000/(0.2998*TMath::Abs(Bz));
          Double_t Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
          if(Lprop/(2*Rcurve)<1) Alpha = -2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
          trackProtonRot.RotateZ(Alpha);
          TLorentzVector trackSigmaplusRot = trackPi0 + trackProtonRot;
          Float_t sigmaplusmassrot = trackSigmaplusRot.M();

          TVector3 sigmamomentumrot(trackSigmaplusRot.Px(),trackSigmaplusRot.Py(),trackSigmaplusRot.Pz());
          TVector3 sigmavertexrot(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
          Lprop = sigmavertex.Mag();
          Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
          if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
          sigmavertexrot.RotateZ(Alpha);          
          Float_t SigmaPointingAnglerot = sigmamomentumrot.Angle(sigmavertexrot);

          Lprop = protonpath.Mag();
          Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
          Alpha = propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
          if(Lprop/(2*Rcurve)<1) Alpha = 2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
          trackSigmaplusRot.RotateZ(Alpha);

          Float_t AntiSigmaPointingAngle = trackProtonRot.Angle(sigmavertexrot);

          if(isReallySigma){
            FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
            FillHistogram("fHistMCInvSigmaMass",sigmaplusmass);
            FillHistogram("fHistMCInvSigmaMassrot",sigmaplusmassrot);
          }
          if(isReallySigma&&isPrimary){
            FillHistogram("fHistMCPrimSigmaPA",SigmaPointingAngle);
            FillHistogram("fHistMCPrimSigmaPArot",SigmaPointingAnglerot);
          }

          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          //Check the deviation from MC
          if(isReallySigma){
            FillHistogram("fHistSigmaPxResnoprop",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResnoprop",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResnoprop",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResprop",1000*(trackSigmaplusRot.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResprop",1000*(trackSigmaplusRot.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResprop",1000*(trackSigmaplusRot.Pz()-MCSigmaMom.Pz()));              
          }

          //MC for Correlations
          fTrackLabel = -999;
          fTrackPDGCode = -999;
          fTrackMotherID = -999;
          fTrackMotherPDGCode = -999;
          fTrackMotherMCPx = -999;
          fTrackMotherMCPy = -999;
          fTrackMotherMCPz = -999;
          fTrackMCPx = -999;
          fTrackMCPy = -999;
          fTrackMCPz = -999;
          fPhoton1Label = -1;
          fPhoton1PDGCode = -999;
          fPhoton1MotherID = -999;
          fPhoton1MotherPDGCode = -999;
          fPhoton1GMotherID = -999;
          fPhoton1GMotherPDGCode = -999;
          fPhoton1MCPx = -999;
          fPhoton1MCPy = -999;
          fPhoton1MCPz = -999;
          fPhoton2Label = -1;
          fPhoton2PDGCode = -999;
          fPhoton2MotherID = -999;
          fPhoton2MotherPDGCode = -999;
          fPhoton2GMotherID = -999;
          fPhoton2GMotherPDGCode = -999;
          fPhoton2MCPx = -999;
          fPhoton2MCPy = -999;
          fPhoton2MCPz = -999;

          if(isMonteCarlo&&fSaveAddMCBranches){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){
              fTrackLabel = ProtonPart->GetLabel();
              fTrackPDGCode = ProtonPart->GetPdgCode();
              fTrackMotherID = ProtonPart->GetMother();
              if(ProtonPart->IsPhysicalPrimary()) fTrackMotherID = -1;
              fTrackMCPx = ProtonPart->Px();
              fTrackMCPy = ProtonPart->Py();
              fTrackMCPz = ProtonPart->Pz();
              AliAODMCParticle* ProtonPartM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ProtonPart->GetMother())));
              if(ProtonPartM){
                fTrackMotherPDGCode = ProtonPartM->GetPdgCode();
                fTrackMotherMCPx = ProtonPartM->Px();
                fTrackMotherMCPy = ProtonPartM->Py();
                fTrackMotherMCPz = ProtonPartM->Pz();
              }//Proton Mother exists
            }//Track exists
            AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
            AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
            AliAODMCParticle* V02Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track3->GetLabel())));
            AliAODMCParticle* V02Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track4->GetLabel())));
            if(V01Daught1&&V01Daught2){
              if(V01Daught1->GetMother()!=-1&&V01Daught1->GetMother()==V01Daught2->GetMother()){
                AliAODMCParticle* V0Part1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
                if(V0Part1){
                  fPhoton1Label = V0Part1->GetLabel();
                  fPhoton1PDGCode = V0Part1->GetPdgCode();
                  fPhoton1MotherID = V0Part1->GetMother();
                  if(V0Part1->IsPhysicalPrimary()) fPhoton1MotherID = -1;
                  fPhoton1MCPx = V0Part1->Px();
                  fPhoton1MCPy = V0Part1->Py();
                  fPhoton1MCPz = V0Part1->Pz();
                  AliAODMCParticle* V0Part1M = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                  if(V0Part1M){
                    fPhoton1MotherPDGCode = V0Part1M->GetPdgCode();
                    fPhoton1GMotherID = V0Part1M->GetMother();
                    if(V0Part1M->IsPhysicalPrimary()) fPhoton1GMotherID = -1;
                    AliAODMCParticle* V0Part1GM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1M->GetMother())));
                    if(V0Part1GM){
                      fPhoton1GMotherPDGCode = V0Part1GM->GetPdgCode();        
                    }//V0 G Mother exists
                  }//V0 Mother exits
                }//MC V0 exists
              }//Real V0
            }//V01s Daughters exist
            if(V02Daught1&&V02Daught2){
              if(V02Daught1->GetMother()!=-1&&V02Daught1->GetMother()==V02Daught2->GetMother()){
                AliAODMCParticle* V0Part2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V02Daught1->GetMother())));
                if(V0Part2){
                  fPhoton2Label = V0Part2->GetLabel();
                  fPhoton2PDGCode = V0Part2->GetPdgCode();
                  fPhoton2MotherID = V0Part2->GetMother();
                  if(V0Part2->IsPhysicalPrimary()) fPhoton2MotherID = -1;
                  fPhoton2MCPx = V0Part2->Px();
                  fPhoton2MCPy = V0Part2->Py();
                  fPhoton2MCPz = V0Part2->Pz();
                  AliAODMCParticle* V0Part2M = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part2->GetMother())));
                  if(V0Part2M){
                    fPhoton2MotherPDGCode = V0Part2M->GetPdgCode();
                    fPhoton2GMotherID = V0Part2M->GetMother();
                    if(V0Part2M->IsPhysicalPrimary()) fPhoton2GMotherID = -1;
                    AliAODMCParticle* V0Part2GM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part2M->GetMother())));
                    if(V0Part2GM){
                     fPhoton2GMotherPDGCode = V0Part2GM->GetPdgCode();        
                    }//V0 G Mother exists
                  }//V0 Mother exits
                }//MC V0 exists
              }//Real V0
            }//V02s Daughters exist
          }//IsMonteCarlo

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          fIsMCDalitz = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          if(isDalitzDecay) fIsMCDalitz = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          fIsV01Onthefly = v0_1->GetOnFlyStatus();
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fHas4DiffIDs = hasdiffindices;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = EventTriggers;
          fSigMCLabel = Pi0MotherLabel;
          fSigProtonID = prot->GetID();
          fSigProtonStatus = prot->GetStatus();
          fSigProtonFilterMap = prot->GetFilterMap();
          fSigEventID = fGlobalEventID;
          fSigCentrality = Centrality;
          fSigRefMultComb05 = fRefMultComb05;
          fSigRefMultComb08 = fRefMultComb08;
          fSigRefMultComb10 = fRefMultComb10;
          fSigBField = Bz;
          fInvSigMass = sigmaplusmass;
          fInvSigpropMass = sigmaplusmassrot;
          fSigY = trackSigmaplus.Rapidity();               
          fSigYprop = trackSigmaplusRot.Rapidity();               
          fSigPA = SigmaPointingAngle; 
          fSigPAprop = SigmaPointingAnglerot; 
          fSigAntiPA = AntiSigmaPointingAngle; 
          fSigCharge = prot->Charge(); 
          fSigPx = trackSigmaplus.Px(); 
          fSigPy = trackSigmaplus.Py(); 
          fSigPz = trackSigmaplus.Pz(); 
          fSigPt = trackSigmaplus.Pt(); 
          fSigPxprop = trackSigmaplusRot.Px(); 
          fSigPyprop = trackSigmaplusRot.Py(); 
          fSigPzprop = trackSigmaplusRot.Pz(); 
          fPrimVertX = primaryVtxPosX; 
          fPrimVertY = primaryVtxPosY; 
          fPrimVertZ = primaryVtxPosZ; 
          fSigDecayVertX = KFSigmaPlus.GetX(); 
          fSigDecayVertY = KFSigmaPlus.GetY(); 
          fSigDecayVertZ = KFSigmaPlus.GetZ();
          fSigFlightDist = TMath::Sqrt((fSigDecayVertX-primaryVtxPosX)*(fSigDecayVertX-primaryVtxPosX)+(fSigDecayVertY-primaryVtxPosY)*(fSigDecayVertY-primaryVtxPosY)+(fSigDecayVertZ-primaryVtxPosZ)*(fSigDecayVertZ-primaryVtxPosZ));          
          fSigDecayVertXMC = MCPi0decayX;
          fSigDecayVertYMC = MCPi0decayY;
          fSigDecayVertZMC = MCPi0decayZ;
          fSigPxMC = MCSigmaMom.Px();        
          fSigPyMC = MCSigmaMom.Py();        
          fSigPzMC = MCSigmaMom.Pz();
          fPrimVertXMC = primaryVtxPosXMC;
          fPrimVertYMC = primaryVtxPosYMC;
          fPrimVertZMC = primaryVtxPosZMC;
          fPhoton1Px = v0_1->Px();
          fPhoton1Py = v0_1->Py();
          fPhoton1Pz = v0_1->Pz();
          fPhoton2Px = v0_2->Px();
          fPhoton2Py = v0_2->Py();
          fPhoton2Pz = v0_2->Pz();
          fPhotonDaughtMaxEta = MaxDaughtEta;
          fPhotonsMaxDeltaTheta = MaxDeltaTheta;
          fPhoton1CPA = v0_1->CosPointingAngle(primaryVtxPos);
          fPhoton2CPA = v0_2->CosPointingAngle(primaryVtxPos);
          fPhoton1Radius = TMath::Sqrt(v0_1->DecayVertexV0X()*v0_1->DecayVertexV0X()+v0_1->DecayVertexV0Y()*v0_1->DecayVertexV0Y());
          fPhoton2Radius = TMath::Sqrt(v0_2->DecayVertexV0X()*v0_2->DecayVertexV0X()+v0_2->DecayVertexV0Y()*v0_2->DecayVertexV0Y()); 
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2;
          fPhoton1DCASV = DCASV1;
          fPhoton2DCASV = DCASV2;
          fTrackDCASV = DCATrack;
          fTrackDCASVKF = DCATrackKF;
          fKFChi2 = KFSigmaPlus.GetChi2();
          fPhotonsMinCluster   = nMinTPCClustDaught;
          fPhotonsMinITSCluster= nMinITSClustDaught;
          fPhotonsMaxalpha     = MaxAlpha;
          fPhotonsMaxqt        = MaxQt;
          fPhotonsMaxOpenAngle = MaxOpenAngle;
          fPhotonsMaxinvmass   = Maxphotonmass;
          fPhotonsMaxNSigTPC   = nMaxNsigTPCDaught;
          fPhotonsMaxChi2      = nMaxPhotchi2;
          fInvPi0Mass = pi0mass; 
          fPi0Px = trackPi0.Px(); 
          fPi0Py = trackPi0.Py(); 
          fPi0Pz = trackPi0.Pz(); 
          fPi0DecayVertX = KFPi0.GetX(); 
          fPi0DecayVertY = KFPi0.GetY(); 
          fPi0DecayVertZ = KFPi0.GetZ();
          fPi0PhotPhotDCA = PhotPhotDCA;
          fProtonPx = prot->Px(); 
          fProtonPy = prot->Py(); 
          fProtonPz = prot->Pz();
          fProtonX = trackxyz[0];
          fProtonY = trackxyz[1];
          fProtonZ = trackxyz[2];
          fProtonEta = prot->Eta(); 
          fProtonpropPx = trackProtonRot.Px(); 
          fProtonpropPy = trackProtonRot.Py(); 
          fProtonpropPz = trackProtonRot.Pz();
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
          fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
          fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
          fProtonNCluster = prot->GetTPCNcls();
          fProtonNITSCluster = prot->GetITSNcls();
          fProtonChi2 = prot->GetTPCchi2();  
          fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
          fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
          fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
          fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
          fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
          fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);

          fEventhasSigmaCand = kTRUE;
          if(fSavePartCand) fSigmaCandTree->Fill();

          //Now apply some Selections to filter out potential Sigma Candidates.
          //If it passes the criteria, Sigma Proton pairs in SE and ME can be written out. 
          //Also all Protons in the Event can be written for offline Event-Mixing.
          if(!fFillPairTreeSE&&!fFillPairTreeME) continue;
          if(pi0mass<fMinCorrPi0Mass) continue;
          if(pi0mass>fMaxCorrPi0Mass) continue;
          if(SigmaPointingAngle>fMaxCorrSigmaPA) continue;
          if(TMath::Abs(DCAxy)<fMinCorrProtonDCAxy) continue;
          if(sigmaplusmass<fMinCorrSigmaMass) continue; 
          if(sigmaplusmass>fMaxCorrSigmaMass) continue;

          if(sigmaplusmass>1.17&&sigmaplusmass<1.21) fEventhasSigma = kTRUE;

          for(Int_t q=0; q<nProton; q++) { 
            AliAODTrack *pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            //if(pairprot->Charge()!=prot->Charge()) continue;

            //If the Sigma made it here, check properties of all Protons in the Event and write them to a tree

            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 deltapvec=sigmamomentum-protonmomentum;
            Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
            Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
            fSigmaProtonkstar = 1000*kstar;

            if(fSigmaProtonkstar>fMaxCorrkstar) continue;

            deltapvec=trackSigmaplusRot.Vect()-protonmomentum;
            SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentumrot.Mag()*sigmamomentumrot.Mag());
            ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
            fSigmaProtonpropkstar = 1000*kstar;

            fPairProtonIsMC = kFALSE;
            fPairProtonIsPrimary = kFALSE;
            if(isMonteCarlo){
              AliAODMCParticle* PairProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(pairprot->GetLabel())));
              if(PairProtonPart){
                if(TMath::Abs(PairProtonPart->GetPdgCode())==2212){
                  fPairProtonIsMC = kTRUE;
                  if(PairProtonPart->IsPhysicalPrimary()) fPairProtonIsPrimary = kTRUE;
                }//MC Particle is a Proton
              }//MC Particle exists 
            }//End of isMonteCarlo

            pairprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

            fPairProtonPx = pairprot->Px();        
            fPairProtonPy = pairprot->Py();        
            fPairProtonPz = pairprot->Pz();
            fPairProtonP = pairprot->P();
            fPairProtonEta = pairprot->Eta();            
            fPairProtonCharge = pairprot->Charge();        
            fPairProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kProton);        
            fPairProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kProton);
            fPairProtNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kPion);
            fPairProtNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kKaon);
            fPairProtNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kElectron);
            fPairProtNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kPion);
            fPairProtNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kKaon);
            fPairProtNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kElectron);
            fPairProtonChi2 = pairprot->GetTPCchi2();    
            fPairProtonCluster = pairprot->GetTPCNcls();
            fPairProtonITSCluster = pairprot->GetITSNcls();
            fPairProtonID = pairprot->GetID();
            fPairProtonStatus = pairprot->GetStatus();
            fPairProtonFilterMap = pairprot->GetFilterMap();

            //Force TOF PID for Pair Proton if requested
            if(pairprot->P()>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

            if(fPairProtonDCAtoPVxy>fMaxCorrPairProtonDCAxy) continue;
            if(fPairProtonDCAtoPVz>fMaxCorrPairProtonDCAz) continue;

            if(q!=k&&fFillPairTreeSE) fSigmaPairTreeSE->Fill();
          }        

          //Continue here if no Event Mixing is requested
          if(!fFillPairTreeME) continue;

          //Get Pool from Pool Manager for given RefMult and Z Vertex values
        	AliEventPool* Evpool = 0x0;
	        if(fEvPoolMgr2&&fUseAbsZCorr)  Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
	        if(fEvPoolMgr2&&!fUseAbsZCorr) Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
		      if(!Evpool){  AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ)); continue;}
          if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << " is empty!\n";*/ continue;}

          //Get Number of Events in Pool
    			Int_t nMixEvents = Evpool->GetCurrentNEvents();
          if(fUseAbsZCorr) FillHistogram("fHistPairNMixedEvents",fRefMultComb08,TMath::Abs(primaryVtxPosZ),nMixEvents);
          else FillHistogram("fHistPairNMixedEvents",fRefMultComb08,primaryVtxPosZ,nMixEvents);

          //Now Loop over the mixed Events
			    for (Int_t iMixEvent = 0; iMixEvent < nMixEvents; iMixEvent++){
            //Retrieve Array of Protons for each mixed event 
				    TObjArray* MixedProtons = (TObjArray*)Evpool->GetEvent(iMixEvent);
      			Int_t nMixProtons = MixedProtons->GetEntriesFast();  			  
            for(Int_t iMixProton = 0; iMixProton < nMixProtons; iMixProton++){

              AliAODTrackcorrelation *mixprot = (AliAODTrackcorrelation*)MixedProtons->At(iMixProton);
              if(!mixprot) continue;

              TVector3 protonmomentum(mixprot->Px(),mixprot->Py(),mixprot->Pz());
              TVector3 deltapvec=sigmamomentum-protonmomentum;
              Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
              Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
              Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
              Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
              Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
              fSigmaProtonkstar = 1000*kstar;

              if(fSigmaProtonkstar>fMaxCorrkstar) continue;

              deltapvec=trackSigmaplusRot.Vect()-protonmomentum;
              SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentumrot.Mag()*sigmamomentumrot.Mag());
              ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
              qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
              vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
              kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
              fSigmaProtonpropkstar = 1000*kstar;

              fPairProtonIsMC = kFALSE;
              fPairProtonIsPrimary = kFALSE;

              mixprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

              fPairProtonPx = mixprot->Px();        
              fPairProtonPy = mixprot->Py();        
              fPairProtonPz = mixprot->Pz();
              fPairProtonCharge = mixprot->Charge();        
              fPairProtonNSigTPC = mixprot->NumberOfSigmasTPCProton();        
              fPairProtonNSigTOF = mixprot->NumberOfSigmasTOFProton();
              fPairProtNSigTPCKaon = mixprot->NumberOfSigmasTPCKaon();        
              fPairProtNSigTOFKaon = mixprot->NumberOfSigmasTOFKaon();
              fPairProtNSigTPCPion = mixprot->NumberOfSigmasTPCPion();        
              fPairProtNSigTOFPion = mixprot->NumberOfSigmasTOFPion();
              fPairProtonChi2 = mixprot->GetTPCchi2();    
              fPairProtonCluster = mixprot->GetTPCNcls();
              fPairProtonITSCluster = mixprot->GetITSNcls();
              fPairProtonID = mixprot->GetID();
              fPairProtonStatus = mixprot->GetStatus();
              fPairProtonFilterMap = mixprot->GetFilterMap();

              TVector3 mprt(fPairProtonPx,fPairProtonPy,fProtonPz);
              fPairProtonP = mprt.Mag();
              fProtonEta = mprt.Eta();

              //Force TOF PID for Pair Proton if requested
              if(fPairProtonP>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

              if(fPairProtonDCAtoPVxy>fMaxCorrPairProtonDCAxy) continue;
              if(fPairProtonDCAtoPVz>fMaxCorrPairProtonDCAz) continue;

              fSigmaPairTreeME->Fill();

            }//End of Loop over Mixed Protons
          }//End of Loop of Mixed Events

        }//End of Proton loop

        //Do Event Mixing for the Background if requested
        if(fSaveMixedBackground){

          Int_t nMixEvents = 0;

          //Get Pool from Pool Manager for given RefMult and Z Vertex values
          AliEventPool* Evpool = 0x0;
  	      if(fEvPoolMgr&&fUseAbsZ)  Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
  	      if(fEvPoolMgr&&!fUseAbsZ) Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
  		    if(!Evpool){AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ)); nMixEvents = -1;}
          else if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << " is empty!\n";*/ nMixEvents = -1;}

          //Get Number of Events in Pool. Number can be reduced with Setter Function to reduce Tree Size.
      		if(nMixEvents!=-1) nMixEvents = Evpool->GetCurrentNEvents();
          else nMixEvents = 0;
          if(fUseAbsZ) FillHistogram("fHistBkgNMixedEvents",fRefMultComb08,TMath::Abs(primaryVtxPosZ),nMixEvents);
          else FillHistogram("fHistBkgNMixedEvents",fRefMultComb08,primaryVtxPosZ,nMixEvents);

          //Now Loop over the mixed Events
  			  for (Int_t iMixEvent = 0; iMixEvent < nMixEvents; iMixEvent++){

            //Retrieve Array of Protons for each mixed event 
  				  TObjArray* MixedProtons = (TObjArray*)Evpool->GetEvent(iMixEvent);
        		Int_t nMixProtons = MixedProtons->GetEntriesFast();  			  
            for(Int_t iMixProton = 0; iMixProton < nMixProtons; iMixProton++){

            AliAODTrackreduced *mixprot = (AliAODTrackreduced*)MixedProtons->At(iMixProton);
            if(!mixprot) continue;

            //Now the hole Sigma reconstruction is repeated with Protons from a different Event.

            trackProton.SetXYZM(mixprot->Px(),mixprot->Py(),mixprot->Pz(),cProtonMass);
            trackSigmaplus = trackPi0 + trackProton;
            Float_t sigmaplusmass = trackSigmaplus.M();
            FillHistogram("fHistInvSigmaMassmix",sigmaplusmass);
            if(sigmaplusmass>fMaxSigmaMass) continue;   //Limit the mass range to reduce tree size

            mixprot->GetXYZ(trackxyz);
            //Move Track to PV of the Pion      
            trackxyz[0]+=primaryVtxPosX;
            trackxyz[1]+=primaryVtxPosY;
            trackxyz[2]+=primaryVtxPosZ;      
            mixprot->GetPxPyPz(trackpxpypz);
            for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
            for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
            mixprot->GetCovarianceXYZPxPyPz(covMatrix);
            if(mixprot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
            else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

            //Reconstruct the Pi0 with the gammas
            KFParticleCD KFPi0CD; //Check Daughters to avoid floating point exceptions. See .h-file
            KFPi0CD.AddDaughter(KFPhoton1);
            if(!KFPi0CD.CheckDaughter(KFPhoton2)) continue;

            KFParticle KFPi0(KFPhoton1,KFPhoton2);
            KFPi0.TransportToDecayVertex();

            KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
            KFSigmaPlusCD.AddDaughter(KFProton);
            if(!KFSigmaPlusCD.CheckDaughter(KFPi0)) continue;

            KFParticle KFSigmaPlus(KFProton,KFPi0);
            KFSigmaPlus.TransportToDecayVertex();
            Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));

            //Calculating DCA of Photon and Proton to Decay Vertex 
            TVector3 SV(KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()); //Sec. Vertex
            Double_t DCASV1 = (p1.Cross((CV1-SV))).Mag()/p1.Mag();                 //DCA to SV of Photons
            Double_t DCASV2 = (p1.Cross((CV2-SV))).Mag()/p2.Mag();                 
            Double_t DCATrackKF = -999; 
            //DCATrackKF = TMath::Abs(KFProton.GetDistanceFromParticle(KFSigmaPlus));

            Float_t  DCAxy = -999., DCAz = -999.;
            mixprot->GetImpactParameters(DCAxy,DCAz);

            if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
            if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
            if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
            if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
            if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
            if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

            TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
            Float_t Rapidity = sig.Rapidity();
            FillHistogram("fHistSigmaYmix",Rapidity);
            if(TMath::Abs(Rapidity)>fMaxSigmaY) continue;  //Cut on the rapidity

            TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
            TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
            Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);

            FillHistogram("fHistSigmaPAmix",SigmaPointingAngle);
            if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

            TLorentzVector trackProtonRot = trackProton;
            TVector3 protonpath(KFSigmaPlus.GetX()-trackxyz[0],KFSigmaPlus.GetY()-trackxyz[1],KFSigmaPlus.GetZ()-trackxyz[2]);
            Double_t propdir = 1; 
            if(TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY))
            <TMath::Sqrt((trackxyz[0]-primaryVtxPosX)*(trackxyz[0]-primaryVtxPosX)+(trackxyz[1]-primaryVtxPosY)*(trackxyz[1]-primaryVtxPosY))) propdir = -1;
            Double_t qprot = 1; if(mixprot->Charge()<0) qprot=-1;
            Double_t Lprop = protonpath.Mag();
            Double_t Rcurve = trackProton.Vect().Mag()*1000/(0.2998*TMath::Abs(Bz));
            Double_t Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
            if(Lprop/(2*Rcurve)<1) Alpha = -2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
            trackProtonRot.RotateZ(Alpha);
            TLorentzVector trackSigmaplusRot = trackPi0 + trackProtonRot;
            Float_t sigmaplusmassrot = trackSigmaplusRot.M();

            TVector3 sigmamomentumrot(trackSigmaplusRot.Px(),trackSigmaplusRot.Py(),trackSigmaplusRot.Pz());
            TVector3 sigmavertexrot(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
            Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
            Lprop = sigmavertex.Mag();
            Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
            if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
            sigmavertexrot.RotateZ(Alpha);          
            Float_t SigmaPointingAnglerot = sigmamomentumrot.Angle(sigmavertexrot);

            Lprop = protonpath.Mag();
            Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
            Alpha = propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
            if(Lprop/(2*Rcurve)<1) Alpha = 2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
            trackSigmaplusRot.RotateZ(Alpha);

            Float_t AntiSigmaPointingAngle = trackProtonRot.Angle(sigmavertexrot);

            // Fill the ME Sigma Trees
            fIsV01Onthefly = v0_1->GetOnFlyStatus();
            fIsV02Onthefly = v0_2->GetOnFlyStatus();
            fHas4DiffIDs = hasdiffindices;
            fSigRunnumber = aodEvent->GetRunNumber();
            fSigTriggerMask = EventTriggers;
            fSigProtonID = mixprot->GetID();
            fSigProtonStatus = mixprot->GetStatus();
            fSigProtonFilterMap = mixprot->GetFilterMap();
            fSigEventID = fGlobalEventID;
            fSigCentrality = Centrality;
            fSigRefMultComb05 = fRefMultComb05;
            fSigRefMultComb08 = fRefMultComb08;
            fSigRefMultComb10 = fRefMultComb10;
            fSigBField = Bz;
            fInvSigMass = sigmaplusmass;
            fInvSigpropMass = sigmaplusmassrot;
            fSigY = trackSigmaplus.Rapidity();               
            fSigYprop = trackSigmaplusRot.Rapidity();                 
            fSigPA = SigmaPointingAngle; 
            fSigPAprop = SigmaPointingAnglerot; 
            fSigAntiPA = AntiSigmaPointingAngle; 
            fSigCharge = mixprot->Charge(); 
            fSigPx = trackSigmaplus.Px(); 
            fSigPy = trackSigmaplus.Py(); 
            fSigPz = trackSigmaplus.Pz(); 
            fSigPt = trackSigmaplus.Pt(); 
            fSigPxprop = trackSigmaplusRot.Px(); 
            fSigPyprop = trackSigmaplusRot.Py(); 
            fSigPzprop = trackSigmaplusRot.Pz(); 
            fPrimVertX = primaryVtxPosX; 
            fPrimVertY = primaryVtxPosY; 
            fPrimVertZ = primaryVtxPosZ; 
            fSigDecayVertX = KFSigmaPlus.GetX(); 
            fSigDecayVertY = KFSigmaPlus.GetY(); 
            fSigDecayVertZ = KFSigmaPlus.GetZ();
            fSigFlightDist = TMath::Sqrt((fSigDecayVertX-primaryVtxPosX)*(fSigDecayVertX-primaryVtxPosX)+(fSigDecayVertY-primaryVtxPosY)*(fSigDecayVertY-primaryVtxPosY)+(fSigDecayVertZ-primaryVtxPosZ)*(fSigDecayVertZ-primaryVtxPosZ));        
            fSigDecayVertXMC = MCPi0decayX;
            fSigDecayVertYMC = MCPi0decayY;
            fSigDecayVertZMC = MCPi0decayZ;
            fSigPxMC = MCSigmaMom.Px();        
            fSigPyMC = MCSigmaMom.Py();        
            fSigPzMC = MCSigmaMom.Pz();        
            fPrimVertXMC = primaryVtxPosXMC;
            fPrimVertYMC = primaryVtxPosYMC;
            fPrimVertZMC = primaryVtxPosZMC;
            fPhoton1Px = v0_1->Px();
            fPhoton1Py = v0_1->Py();
            fPhoton1Pz = v0_1->Pz();
            fPhoton2Px = v0_2->Px();
            fPhoton2Py = v0_2->Py();
            fPhoton2Pz = v0_2->Pz();
            fPhotonDaughtMaxEta = MaxDaughtEta;
            fPhotonsMaxDeltaTheta = MaxDeltaTheta;
            fPhoton1CPA = v0_1->CosPointingAngle(primaryVtxPos);
            fPhoton2CPA = v0_2->CosPointingAngle(primaryVtxPos);
            fPhoton1Radius = TMath::Sqrt(v0_1->DecayVertexV0X()*v0_1->DecayVertexV0X()+v0_1->DecayVertexV0Y()*v0_1->DecayVertexV0Y());
            fPhoton2Radius = TMath::Sqrt(v0_2->DecayVertexV0X()*v0_2->DecayVertexV0X()+v0_2->DecayVertexV0Y()*v0_2->DecayVertexV0Y()); 
            fPhoton1DCAPV = DCAPV1;
            fPhoton2DCAPV = DCAPV2;
            fPhoton1DCASV = DCASV1;
            fPhoton2DCASV = DCASV2;
            fTrackDCASV = -999;
            fTrackDCASVKF = DCATrackKF;
            fKFChi2 = KFSigmaPlus.GetChi2();
            fPhotonsMinCluster   = nMinTPCClustDaught;
            fPhotonsMinITSCluster   = nMinITSClustDaught;
            fPhotonsMaxalpha     = MaxAlpha;
            fPhotonsMaxqt        = MaxQt;
            fPhotonsMaxOpenAngle = MaxOpenAngle;
            fPhotonsMaxinvmass   = Maxphotonmass;
            fPhotonsMaxNSigTPC   = nMaxNsigTPCDaught;
            fPhotonsMaxChi2      = nMaxPhotchi2;
            fInvPi0Mass = pi0mass; 
            fPi0Px = trackPi0.Px(); 
            fPi0Py = trackPi0.Py(); 
            fPi0Pz = trackPi0.Pz(); 
            fPi0DecayVertX = KFPi0.GetX(); 
            fPi0DecayVertY = KFPi0.GetY(); 
            fPi0DecayVertZ = KFPi0.GetZ();
            fPi0PhotPhotDCA = PhotPhotDCA;
            fProtonPx = mixprot->Px(); 
            fProtonPy = mixprot->Py(); 
            fProtonPz = mixprot->Pz();
            fProtonX = trackxyz[0];
            fProtonY = trackxyz[1];
            fProtonZ = trackxyz[2];
            TVector3 prt(fProtonPx,fProtonPy,fProtonPz);
            fProtonEta = prt.Eta();
            fProtonpropPx = trackProtonRot.Px(); 
            fProtonpropPy = trackProtonRot.Py(); 
            fProtonpropPz = trackProtonRot.Pz();
            fProtonDCAtoPVxy = DCAxy; 
            fProtonDCAtoPVz = DCAz; 
            fProtonPi0DCA = ProtPi0DCA;
            fProtonNSigTPC = mixprot->NumberOfSigmasTPCProton();
            fProtonNSigTOF = mixprot->NumberOfSigmasTOFProton();
            fProtonNCluster = mixprot->GetTPCNcls();
            fProtonNITSCluster = mixprot->GetITSNcls();
            fProtonChi2 = mixprot->GetTPCchi2();  
            fProtonNSigTPCPion = mixprot->NumberOfSigmasTPCPion();
            fProtonNSigTPCKaon = mixprot->NumberOfSigmasTPCKaon();
            fProtonNSigTPCElec = mixprot->NumberOfSigmasTPCElectron();
            fProtonNSigTOFPion = mixprot->NumberOfSigmasTOFPion();
            fProtonNSigTOFKaon = mixprot->NumberOfSigmasTOFKaon();
            fProtonNSigTOFElec = mixprot->NumberOfSigmasTOFElectron();
            fSigmaMEBackgroundTree->Fill();

            }//End of Mixed Proton loop
          }//End of Mixed Event loop

        }////End of Mixed Event Background

      /************************End of Sigma+ reconstruction*****************************/        

    }//End of Photon 2 loop
  }//End of Photon 1 loop
  
return;

} //End of ReconstructParticles()

//________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesPHOS() {

  KFParticle::SetField(Bz);

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TVector3 ConvPhoton, ClusterPhoton; 
  TLorentzVector trackPhoton1, trackPhoton2, trackPi0, trackProton, trackSigmaplus;
  KFParticle KFElectron, KFPositron, KFProton; 

  const Int_t nCaloPhoton = fCaloPhotonArray.size();
  const Int_t nConvPhoton = fConvPhotonArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nConvPhoton; i++) {

    AliAODv0 *v0_1 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(i));
    if(!v0_1) continue;

    // Get daughter tracks      
    AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(0));
    AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(1));

    if(!track1 || !track2) {
      AliWarning("ERROR: Could not retrieve all AOD tracks in Pi0 reconstruction!");
      continue;
    }

    //KF calculations
    // Set up KFParticle
    trackparams[0] = v0_1->DecayVertexV0X();
    trackparams[1] = v0_1->DecayVertexV0Y();
    trackparams[2] = v0_1->DecayVertexV0Z();
    trackparams[3] = v0_1->MomPosX();
    trackparams[4] = v0_1->MomPosY();
    trackparams[5] = v0_1->MomPosZ();
    if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
    else track2->GetCovarianceXYZPxPyPz(covMatrix);
    KFPositron.Create(trackparams,covMatrix,1,cElectronMass);

    // Repeat for all other particles
    trackparams[0] = v0_1->DecayVertexV0X();
    trackparams[1] = v0_1->DecayVertexV0Y();
    trackparams[2] = v0_1->DecayVertexV0Z();
    trackparams[3] = v0_1->MomNegX();
    trackparams[4] = v0_1->MomNegY();
    trackparams[5] = v0_1->MomNegZ();
    if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
    else track2->GetCovarianceXYZPxPyPz(covMatrix);
    KFElectron.Create(trackparams,covMatrix,-1,cElectronMass);

    KFParticleCD KFPhotonCD; //Check Daughters to avoid floating point exceptions. See .h-file
    KFPhotonCD.AddDaughter(KFElectron);
    if(!KFPhotonCD.CheckDaughter(KFPositron)) {FillHistogram("fHistPhotonKFCheckPHOS",1); continue;}
    else{ FillHistogram("fHistPhotonKFCheckPHOS",3);}

    //Reconstruct the Photons with two Methods: Standard and special gamma reconstruction
    KFParticle KFPhoton(KFElectron,KFPositron);

    //Transport Photon to Conversion Point
    KFPhoton.TransportToDecayVertex();

    //Calculating DCA of Photon to PV 
    TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);                            //Prim. Vertex
    TVector3 CV1(v0_1->DecayVertexV0X(),v0_1->DecayVertexV0Y(),v0_1->DecayVertexV0Z()); //Conv. Vertices
    TVector3 p1(v0_1->Px(),v0_1->Py(),v0_1->Pz());                     //Momentum vectors of the photons
    Double_t DCAPV1 = (p1.Cross((CV1-PV))).Mag()/p1.Mag();                        //DCA to PV of Photons

    // Save Photon quality
    Int_t nMinTPCClustDaught = track1->GetTPCNcls();
    if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();

    Int_t nMinITSClustDaught = track1->GetITSNcls();
    if(track2->GetITSNcls()<nMinITSClustDaught) nMinITSClustDaught = track2->GetITSNcls();

    Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));

    Double_t MaxAlpha = TMath::Abs(v0_1->AlphaV0());
    Double_t MaxQt = TMath::Abs(v0_1->PtArmV0());
    Double_t MaxOpenAngle = v0_1->OpenAngleV0();

    TLorentzVector tl1, tl2, tlp;
    tl1.SetXYZM(v0_1->MomNegX(),v0_1->MomNegY(),v0_1->MomNegZ(),cElectronMass);
    tl2.SetXYZM(v0_1->MomPosX(),v0_1->MomPosY(),v0_1->MomPosZ(),cElectronMass);
    tlp=tl1+tl2;
    Double_t Maxphotonmass = tlp.M();

    Double_t nMaxPhotchi2 = v0_1->Chi2V0();

    Double_t MaxDaughtEta = TMath::Abs(track1->Eta());
    if(TMath::Abs(track2->Eta())>MaxDaughtEta) MaxDaughtEta = TMath::Abs(track2->Eta());

    Double_t MaxDeltaTheta = TMath::Abs(track1->Theta()-track2->Theta());

    //Set TLorentzvector for further calculations  
    trackPhoton1.SetXYZM(v0_1->Px(),v0_1->Py(),v0_1->Pz(),0);    
    //Set also TVector3 (better suited for some applications)
    ConvPhoton.SetXYZ(v0_1->Px(),v0_1->Py(),v0_1->Pz());

    for(Int_t k=0; k<nProton; k++) {

      AliAODTrack *prot;
      prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(k));
      if(!prot) continue;

      if(fCheckProtonV0IDs&&(prot->GetID()==track1->GetID()||prot->GetID()==track2->GetID())) continue;

      Float_t  DCAxy = -999., DCAz = -999.;
      prot->GetImpactParameters(DCAxy,DCAz);

      if(TMath::Abs(DCAxy)<fMinProtonDCAxyPHOS) continue; //Coarse cut on DCA to PV to reduce the tree size
      if(TMath::Abs(DCAz)<fMinProtonDCAzPHOS) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxyPHOS&&DCAxy!=-999) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAzPHOS&&DCAz!=-999) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxyPHOS&&fRequireDCACutPHOS) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAzPHOS&&fRequireDCACutPHOS) continue;

      //Set TLorentzvector for further calculations  
      trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);

      prot->GetXYZ(trackxyz);      
      prot->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      prot->GetCovarianceXYZPxPyPz(covMatrix);
      if(prot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
      else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

      //Reconstruct pseudo sigma plus from proton and photon to access the decay vertex.
      KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
      KFSigmaPlusCD.AddDaughter(KFProton);
      if(!KFSigmaPlusCD.CheckDaughter(KFPhoton)) {FillHistogram("fHistSigmaKFCheckPHOS",1); continue;}
      else{ FillHistogram("fHistSigmaKFCheckPHOS",3);}

      KFParticle KFSigmaPlus(KFProton,KFPhoton);
      KFSigmaPlus.TransportToDecayVertex();

      TVector3 PhotonXDir(v0_1->DecayVertexV0X()-KFSigmaPlus.GetX(),v0_1->DecayVertexV0Y()-KFSigmaPlus.GetY(),v0_1->DecayVertexV0Z()-KFSigmaPlus.GetZ());
      TVector3 PhotonPDir(v0_1->Px(),v0_1->Py(),v0_1->Pz());
      Double_t PhotonXPAngle = PhotonXDir.Angle(PhotonPDir); 

      //Calculate distance between Photon and Proton
      Double_t ProtPhotonDCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPhoton));
      //Calculate distance to primary vertex and reject candidates which are too close to the PV
      Double_t SigmaDisttoPV = TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY)+(KFSigmaPlus.GetZ()-primaryVtxPosZ)*(KFSigmaPlus.GetZ()-primaryVtxPosZ));
      Double_t SigmaRadius = TMath::Sqrt((KFSigmaPlus.GetX()*KFSigmaPlus.GetX())+(KFSigmaPlus.GetY()*KFSigmaPlus.GetY()));

      //Calculating DCA of Photon and Proton to Decay Vertex 
      TVector3 SV(KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()); //Sec. Vertex
      Double_t DCASV1 = (p1.Cross((CV1-SV))).Mag()/p1.Mag();                 //DCA to SV of Photons
      Double_t Trackpnt[3];
      prot->GetXYZatR(SigmaRadius,Bz,Trackpnt,0);
      TVector3 TrackpntVec(Trackpnt[0],Trackpnt[1],Trackpnt[2]);  
      Double_t DCATrack = (TrackpntVec-SV).Mag();
      Double_t DCATrackKF = -999; 
      //DCATrackKF = TMath::Abs(KFProton.GetDistanceFromParticle(KFSigmaPlus));

      //Clone proton and propagate it to the secondary vertex
      TLorentzVector trackProtonRot = trackProton;
      TVector3 protonpath(KFSigmaPlus.GetX()-trackxyz[0],KFSigmaPlus.GetY()-trackxyz[1],KFSigmaPlus.GetZ()-trackxyz[2]);
      Double_t propdir = 1; 
      if(TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY))
      <TMath::Sqrt((trackxyz[0]-primaryVtxPosX)*(trackxyz[0]-primaryVtxPosX)+(trackxyz[1]-primaryVtxPosY)*(trackxyz[1]-primaryVtxPosY))) propdir = -1;
      Double_t qprot = 1; if(prot->Charge()<0) qprot=-1;
      Double_t Lprop = protonpath.Mag();
      Double_t Rcurve = trackProton.Vect().Mag()*1000/(0.2998*TMath::Abs(Bz));
      Double_t Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
      if(Lprop/(2*Rcurve)<1) Alpha = -2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
      trackProtonRot.RotateZ(Alpha);

      for(Int_t j=0; j<nCaloPhoton; j++) {

        AliAODCaloCluster *aodCluster = dynamic_cast<AliAODCaloCluster*>(aodEvent->GetCaloCluster(fCaloPhotonArray.at(j)));
        if(!aodCluster) continue;

        //MC treatment
        Bool_t   isReallyPi0 = kFALSE; 
        Double_t MCPi0Vtx[3];
        Int_t    clustPhotPDGcode;
        Double_t ClustEMC;
        Int_t    SigmaMCLabel=-1;
        Bool_t   isReallySigma = kFALSE;
        Bool_t   isDalitzDecay = kFALSE;
        Bool_t   isPrimary = kFALSE;
        TLorentzVector MCSigmaMom, MCCaloMom; 

        if(isMonteCarlo){
          AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
          AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
          AliAODMCParticle* ClustPart  = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodCluster->GetLabel())));
          AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));

          Int_t Pi0ID = -1;
          if(V01Daught1&&V01Daught2){
            if(V01Daught1->GetMother()!=-1&&TMath::Abs(V01Daught1->GetMother())==TMath::Abs(V01Daught2->GetMother())){
              AliAODMCParticle* V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
              if(V0Part){
                if(V0Part->GetPdgCode()==22){
                  Pi0ID = V0Part->GetMother();
                }//V0 is photon
                if(V0Part->GetPdgCode()==111){
                  Pi0ID = V0Part->Label();
                  isDalitzDecay = kTRUE;
                }//V0 is from dalitz decay
              }//V0 Particle exists
            }//V0 Daughters have common mother
          }//V0 Daughters exist

          Int_t SigmaID = -1;
          if(ClustPart){
            clustPhotPDGcode = ClustPart->GetPdgCode();
            ClustEMC = ClustPart->E();
            MCCaloMom.SetXYZM(ClustPart->Px(),ClustPart->Py(),ClustPart->Pz(),ClustPart->M());
            if(clustPhotPDGcode==22){
              if(ClustPart->GetMother()!=-1&&ClustPart->GetMother()==Pi0ID){
                AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0ID)));
                if(Pi0Part){
                  if(Pi0Part->GetPdgCode()==111){
                    isReallyPi0 = kTRUE; 
                    MCPi0Vtx[0]=Pi0Part->Xv();
                    MCPi0Vtx[1]=Pi0Part->Yv();
                    MCPi0Vtx[2]=Pi0Part->Zv();
                    SigmaID = Pi0Part->GetMother();
                  }//Particle is Pi0
                }//Particle exists
              }//Cluster is from same Particle
            }//Cluster is photon
          }//Cluster exists

          if(ProtonPart){
            if(TMath::Abs(ProtonPart->GetPdgCode())==2212&&ProtonPart->GetMother()!=-1&&ProtonPart->GetMother()==SigmaID){
              AliAODMCParticle* SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(SigmaID)));
              if(SigmaPart){
                if(TMath::Abs(SigmaPart->GetPdgCode())==3222){
                  SigmaMCLabel = SigmaID;
                  isReallySigma = kTRUE;
                  if(SigmaPart->IsPhysicalPrimary()) isPrimary = kTRUE; 
                  MCSigmaMom.SetXYZM(SigmaPart->Px(),SigmaPart->Py(),SigmaPart->Pz(),cSigmaMass);
                }//Particle is Sigma+
              }//Sigma exists
            }//Particle is Proton and has common mother with pion 
          }//Proton exists

        }//End of isMonteCarlo

        Float_t clustE = aodCluster->E();
        Float_t clustpos[3];
        aodCluster->GetPosition(clustpos);

        //Calculating beta of cluster 
        ClusterPhoton.SetXYZ(clustpos[0],clustpos[1],clustpos[2]);
        Double_t FlightDist = ClusterPhoton.Mag();
        Double_t vTOF = 0;
        if(aodCluster->GetTOF()!=0) vTOF = FlightDist/(1e12*aodCluster->GetTOF());
        Float_t clustbeta = vTOF/c;

        //Set Calo Photon. Best assumtion: Origin = Secondary Vertex
        ClusterPhoton.SetXYZ(clustpos[0]-KFSigmaPlus.GetX(),clustpos[1]-KFSigmaPlus.GetY(),clustpos[2]-KFSigmaPlus.GetZ());
        //Calculate pion mass from available information  
        Double_t pi0mass = TMath::Sqrt(2*clustE*ConvPhoton.Mag()*(1-TMath::Cos(ConvPhoton.Angle(ClusterPhoton))));

        //Normalize Calo Photon and scale to cluster energy
        ClusterPhoton*=clustE/ClusterPhoton.Mag();
        trackPhoton2.SetXYZM(ClusterPhoton.X(),ClusterPhoton.Y(),ClusterPhoton.Z(),0);            
        trackSigmaplus = trackPhoton1 + trackPhoton2 + trackProton;
        Float_t sigmamassuncorr = trackSigmaplus.M();  //Save uncorrected mass

    	  FillHistogram("fHistGammaPairInvMassPHOS",pi0mass);
    	  if(v0_1->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassOnflyPHOS",pi0mass);
        if(isReallyPi0){
    	    FillHistogram("fHistGammaPairInvMassMCPHOS",pi0mass);
    	    if(v0_1->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassMCOnflyPHOS",pi0mass);
        }
        if(pi0mass<fMinPi0MassPHOS || pi0mass>fMaxPi0MassPHOS) continue;  //Coarse mass cut to reduce combinatorics!

        //Reset Calo Photon.
        ClusterPhoton.SetXYZ(clustpos[0]-KFSigmaPlus.GetX(),clustpos[1]-KFSigmaPlus.GetY(),clustpos[2]-KFSigmaPlus.GetZ());
        //Calculate cluster energy using nominal pion mass
        Double_t ClustECorr = cPi0Mass*cPi0Mass/(2*ConvPhoton.Mag()*(1-TMath::Cos(ConvPhoton.Angle(ClusterPhoton))));       
        //Normalize Calo Photon and scale to corrected energy
        ClusterPhoton*=ClustECorr/ClusterPhoton.Mag();
        //Set TLorentzvector for further calculations  
        trackPhoton2.SetXYZM(ClusterPhoton.X(),ClusterPhoton.Y(),ClusterPhoton.Z(),0);            

        trackPi0 = trackPhoton1 + trackPhoton2;  
        trackSigmaplus = trackPi0 + trackProton;

        //Calculate sigma mass
        Float_t sigmaplusmass = trackSigmaplus.M();
        FillHistogram("fHistInvSigmaMassPHOS",sigmaplusmass);
        FillHistogram("fHistInvSigmaMassPHOSUncorr",sigmamassuncorr);
        if(sigmaplusmass>fMaxSigmaMassPHOS) continue;   //Limit the mass range to reduce tree size

        TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
        Float_t Rapidity = sig.Rapidity();

        TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
        TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
        Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);

        TLorentzVector trackSigmaplusRot = trackPi0 + trackProtonRot;
        TVector3 sigmamomentumrot(trackSigmaplusRot.Px(),trackSigmaplusRot.Py(),trackSigmaplusRot.Pz());
        TVector3 sigmavertexrot(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
        Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
        Lprop = sigmavertex.Mag();
        Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
        if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
        sigmavertexrot.RotateZ(Alpha);          
        Float_t AntiSigmaPointingAngle = trackProtonRot.Angle(sigmavertexrot); //Fake Pointing Angle

        Lprop = protonpath.Mag();
        Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
        Alpha = propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
        if(Lprop/(2*Rcurve)<1) Alpha = 2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
        trackSigmaplusRot.RotateZ(Alpha);
        Float_t SigmaPointingAnglerot = sigmamomentumrot.Angle(sigmavertexrot); 
        Float_t sigmaplusmassrot = trackSigmaplusRot.M();

        //Fill some histograms
        FillHistogram("fHistSigmaDCAtoPVPHOS",SigmaDisttoPV);          
        FillHistogram("fHistSigmaRadiusPHOS",SigmaRadius);          
        FillHistogram("fHistProtPhotonDCAPHOS",ProtPhotonDCA);          
        FillHistogram("fHistSigmaYPHOS",Rapidity);
        FillHistogram("fHistSigmaPAPHOS",SigmaPointingAngle);
        FillHistogram("fHistSigmaPAPHOSrot",SigmaPointingAnglerot);
        FillHistogram("fHistSigmaAntiPAPHOS",AntiSigmaPointingAngle);
        FillHistogram("fHistPhotonSecPAPHOS",PhotonXPAngle);
        if(isReallySigma){
          FillHistogram("fHistSigmaDCAtoPVPHOSMC",SigmaDisttoPV);          
          FillHistogram("fHistProtPhotonDCAPHOSMC",ProtPhotonDCA);          
          FillHistogram("fHistMCSigmaYPHOS",Rapidity);
          FillHistogram("fHistMCSigmaPAPHOS",SigmaPointingAngle);
          FillHistogram("fHistMCSigmaAntiPAPHOS",AntiSigmaPointingAngle);
          FillHistogram("fHistPhotonSecPAPHOSMC",PhotonXPAngle);
          FillHistogram("fHistMCInvSigmaMassPHOSUncorr",sigmamassuncorr);
          FillHistogram("fHistMCInvSigmaMassPHOS",sigmaplusmass);
          FillHistogram("fHistMCInvSigmaMassPHOSrot",sigmaplusmassrot);
          FillHistogram("fHistKFSigmaVertexResXPHOS",KFSigmaPlus.GetX()-MCPi0Vtx[0]);            
          FillHistogram("fHistKFSigmaVertexResYPHOS",KFSigmaPlus.GetY()-MCPi0Vtx[1]);            
          FillHistogram("fHistKFSigmaVertexResZPHOS",KFSigmaPlus.GetZ()-MCPi0Vtx[2]);              
          FillHistogram("fHistSigmaPxResnopropPHOS",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
          FillHistogram("fHistSigmaPyResnopropPHOS",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
          FillHistogram("fHistSigmaPzResnopropPHOS",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
          FillHistogram("fHistSigmaPxRespropPHOS",1000*(trackSigmaplusRot.Px()-MCSigmaMom.Px()));            
          FillHistogram("fHistSigmaPyRespropPHOS",1000*(trackSigmaplusRot.Py()-MCSigmaMom.Py()));            
          FillHistogram("fHistSigmaPzRespropPHOS",1000*(trackSigmaplusRot.Pz()-MCSigmaMom.Pz()));              
          if(isPrimary){
            FillHistogram("fHistMCPrimSigmaPAPHOS",SigmaPointingAngle);
            FillHistogram("fHistMCPrimSigmaPAPHOSrot",SigmaPointingAnglerot);
            FillHistogram("fHistMCPrimSigmaAntiPAPHOS",AntiSigmaPointingAngle);          
          }
        }

        //Coarse topological cuts to reduce the tree size
        if(SigmaDisttoPV<fMinSigmaDCAtoPVPHOS||SigmaDisttoPV>fMaxSigmaDCAtoPVPHOS) continue;
        if(TMath::Abs(Rapidity)>fMaxSigmaYPHOS) continue;  //Cut on the rapidity
        if(SigmaPointingAngle>fMaxSigmaPAPHOS) continue;  
        if(EventTriggers&65536||EventTriggers&8){
          if(!isMonteCarlo&&SigmaPointingAngle>fMaxSigmaPAPHOSHM) continue;  //Apply stricter cut in case of high multiplicity          
        }
        if(AntiSigmaPointingAngle<fMinSigmaAntiPAPHOS) continue; 
        if(ProtPhotonDCA>fMaxProtPhotDCA) continue;

        //MC for Correlations
        fTrackLabel = -999;
        fTrackPDGCode = -999;
        fTrackMotherID = -999;
        fTrackMotherPDGCode = -999;
        fTrackMotherMCPx = -999;
        fTrackMotherMCPy = -999;
        fTrackMotherMCPz = -999;
        fTrackMCPx = -999;
        fTrackMCPy = -999;
        fTrackMCPz = -999;
        fPhoton1Label = -1;
        fPhoton1PDGCode = -999;
        fPhoton1MotherID = -999;
        fPhoton1MotherPDGCode = -999;
        fPhoton1GMotherID = -999;
        fPhoton1GMotherPDGCode = -999;
        fPhoton1MCPx = -999;
        fPhoton1MCPy = -999;
        fPhoton1MCPz = -999;
        fPhoton2Label = -1;
        fPhoton2PDGCode = -999;
        fPhoton2MotherID = -999;
        fPhoton2MotherPDGCode = -999;
        fPhoton2GMotherID = -999;
        fPhoton2GMotherPDGCode = -999;
        fPhoton2MCPx = -999;
        fPhoton2MCPy = -999;
        fPhoton2MCPz = -999;

        if(isMonteCarlo&&fSaveAddMCBranches){
          AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
          if(ProtonPart){
            fTrackLabel = ProtonPart->GetLabel();
            fTrackPDGCode = ProtonPart->GetPdgCode();
            fTrackMotherID = ProtonPart->GetMother();
            if(ProtonPart->IsPhysicalPrimary()) fTrackMotherID = -1;
            fTrackMCPx = ProtonPart->Px();
            fTrackMCPy = ProtonPart->Py();
            fTrackMCPz = ProtonPart->Pz();
            AliAODMCParticle* ProtonPartM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ProtonPart->GetMother())));
            if(ProtonPartM){
              fTrackMotherPDGCode = ProtonPartM->GetPdgCode();
              fTrackMotherMCPx = ProtonPartM->Px();
              fTrackMotherMCPy = ProtonPartM->Py();
              fTrackMotherMCPz = ProtonPartM->Pz();
            }//Proton Mother exists
          }//Track exists
          AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
          AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
          AliAODMCParticle* ClustPart  = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodCluster->GetLabel())));
          if(V01Daught1&&V01Daught2){
            if(V01Daught1->GetMother()!=-1&&V01Daught1->GetMother()==V01Daught2->GetMother()){
              AliAODMCParticle* V0Part1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
              if(V0Part1){
                fPhoton1Label = V0Part1->GetLabel();
                fPhoton1PDGCode = V0Part1->GetPdgCode();
                fPhoton1MotherID = V0Part1->GetMother();
                if(V0Part1->IsPhysicalPrimary()) fPhoton1MotherID = -1;
                fPhoton1MCPx = V0Part1->Px();
                fPhoton1MCPy = V0Part1->Py();
                fPhoton1MCPz = V0Part1->Pz();
                AliAODMCParticle* V0Part1M = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                if(V0Part1M){
                  fPhoton1MotherPDGCode = V0Part1M->GetPdgCode();
                  fPhoton1GMotherID = V0Part1M->GetMother();
                  if(V0Part1M->IsPhysicalPrimary()) fPhoton1GMotherID = -1;
                  AliAODMCParticle* V0Part1GM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1M->GetMother())));
                  if(V0Part1GM){
                    fPhoton1GMotherPDGCode = V0Part1GM->GetPdgCode();        
                  }//V0 G Mother exists
                }//V0 Mother exits
              }//MC V0 exists
            }//Real V0
          }//V01s Daughters exist
          if(ClustPart){
            fPhoton2Label = ClustPart->GetLabel();
            fPhoton2PDGCode = ClustPart->GetPdgCode();
            fPhoton2MotherID = ClustPart->GetMother();
            if(ClustPart->IsPhysicalPrimary()) fPhoton2MotherID = -1;
            fPhoton2MCPx = ClustPart->Px();
            fPhoton2MCPy = ClustPart->Py();
            fPhoton2MCPz = ClustPart->Pz();
            AliAODMCParticle* ClustPartM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ClustPart->GetMother())));
            if(ClustPartM){
              fPhoton2MotherPDGCode = ClustPartM->GetPdgCode();
              fPhoton2GMotherID = ClustPartM->GetMother();
              if(ClustPartM->IsPhysicalPrimary()) fPhoton2GMotherID = -1;
              AliAODMCParticle* ClustPartGM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ClustPartM->GetMother())));
              if(ClustPartGM){
                fPhoton2GMotherPDGCode = ClustPartGM->GetPdgCode();        
              }//Cluster G Mother exists
            }//Cluster Mother exits
          }//MC Cluster exists
        }//IsMonteCarlo

        // Fill the Sigma Candidate Trees
        fIsMCSigma = kFALSE; 
        fIsMCDalitz = kFALSE; 
        if(isReallySigma) fIsMCSigma = kTRUE;
        if(isDalitzDecay) fIsMCDalitz = kTRUE;
        fIsMCPrimary = kFALSE;
        if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
        fIsV01Onthefly = v0_1->GetOnFlyStatus();
        fIsClusterEMCAL = aodCluster->IsEMCAL();
        fIsClusterPHOS = aodCluster->IsPHOS();
        fSigRunnumber = aodEvent->GetRunNumber();
        fSigTriggerMask = EventTriggers;
        fSigMCLabel = SigmaMCLabel;
        fSigProtonID = prot->GetID();
        fSigProtonStatus = prot->GetStatus();
        fSigProtonFilterMap = prot->GetFilterMap();
        fSigEventID = fGlobalEventID;
        fSigCentrality = Centrality;
        fSigRefMultComb05 = fRefMultComb05;
        fSigRefMultComb08 = fRefMultComb08;
        fSigRefMultComb10 = fRefMultComb10;
        fSigBField = Bz;
        fInvSigMass = sigmaplusmass;
        fInvSigpropMass = sigmaplusmassrot;
        fInvSigMassUncorr = sigmamassuncorr;
        fSigY = trackSigmaplus.Rapidity();               
        fSigYprop = trackSigmaplusRot.Rapidity();                
        fSigPA = SigmaPointingAngle; 
        fSigPAprop = SigmaPointingAnglerot; 
        fSigAntiPA = AntiSigmaPointingAngle; 
        fSigCharge = prot->Charge(); 
        fSigPx = trackSigmaplus.Px(); 
        fSigPy = trackSigmaplus.Py(); 
        fSigPz = trackSigmaplus.Pz(); 
        fSigPt = trackSigmaplus.Pt(); 
        fSigPxprop = trackSigmaplusRot.Px(); 
        fSigPyprop = trackSigmaplusRot.Py(); 
        fSigPzprop = trackSigmaplusRot.Pz(); 
        fPrimVertX = primaryVtxPosX; 
        fPrimVertY = primaryVtxPosY; 
        fPrimVertZ = primaryVtxPosZ; 
        fSigDecayVertX = KFSigmaPlus.GetX(); 
        fSigDecayVertY = KFSigmaPlus.GetY(); 
        fSigDecayVertZ = KFSigmaPlus.GetZ();
        fSigFlightDist = SigmaDisttoPV;
        fSigDecayVertXMC = MCPi0Vtx[0];
        fSigDecayVertYMC = MCPi0Vtx[1];
        fSigDecayVertZMC = MCPi0Vtx[2];
        fSigPxMC = MCSigmaMom.Px();        
        fSigPyMC = MCSigmaMom.Py();        
        fSigPzMC = MCSigmaMom.Pz();
        fPrimVertXMC = primaryVtxPosXMC;
        fPrimVertYMC = primaryVtxPosYMC;
        fPrimVertZMC = primaryVtxPosZMC;
        fPhoton1Px = v0_1->Px();
        fPhoton1Py = v0_1->Py();
        fPhoton1Pz = v0_1->Pz();
        fPhoton2Px = trackPhoton2.Px();
        fPhoton2Py = trackPhoton2.Py();
        fPhoton2Pz = trackPhoton2.Pz();
        fCaloPhotonX = clustpos[0];
        fCaloPhotonY = clustpos[1];
        fCaloPhotonZ = clustpos[2];
        fConvPhotonX = v0_1->DecayVertexV0X();
        fConvPhotonY = v0_1->DecayVertexV0Y();
        fConvPhotonZ = v0_1->DecayVertexV0Z();
        fConvPhotonSecPA = PhotonXPAngle;
        fCaloPhotonPxMC = MCCaloMom.Px();
        fCaloPhotonPyMC = MCCaloMom.Py();
        fCaloPhotonPzMC = MCCaloMom.Pz();
        fCaloPhotonE = clustE;
        fCaloPhotonEcorr = ClustECorr;
        fCaloPhotonEMC = ClustEMC;
        fClustNLabels = aodCluster->GetNLabels();
        fClustPDGCode = clustPhotPDGcode;
        fClustDispersion = aodCluster->GetDispersion();
        fClustM20 = aodCluster->GetM20();
        fClustM02 = aodCluster->GetM02();
        fClustNTracksMatched = aodCluster->GetNTracksMatched();
        fClustTrackDx = aodCluster->GetTrackDx();
        fClustTrackDz = aodCluster->GetTrackDz();
        fClustTrackD = TMath::Sqrt(fClustTrackDx*fClustTrackDx+fClustTrackDz*fClustTrackDz);
        fClustTOF = 1000000000*aodCluster->GetTOF();
        fClustBeta = clustbeta;
        fClustNCells = aodCluster->GetNCells();
        fClustDisttoBC = aodCluster->GetDistanceToBadChannel();
        for(Int_t iCell=0; iCell<20; iCell++){fCellsAbsId[iCell] = (Short_t)aodCluster->GetCellAbsId(iCell);}
        fPhotonDaughtMaxEta = MaxDaughtEta;
        fPhotonsMaxDeltaTheta = MaxDeltaTheta;
        fPhoton1CPA = v0_1->CosPointingAngle(primaryVtxPos);
        fPhoton1Radius = TMath::Sqrt(v0_1->DecayVertexV0X()*v0_1->DecayVertexV0X()+v0_1->DecayVertexV0Y()*v0_1->DecayVertexV0Y());
        fPhoton1DCAPV = DCAPV1;
        fPhoton1DCASV = DCASV1;
        fTrackDCASV = DCATrack;
        fTrackDCASVKF = DCATrackKF;
        fKFChi2 = KFSigmaPlus.GetChi2();
        fPhotonsMinCluster   = nMinTPCClustDaught;
        fPhotonsMinITSCluster= nMinITSClustDaught;
        fPhotonsMaxalpha     = MaxAlpha;
        fPhotonsMaxqt        = MaxQt;
        fPhotonsMaxOpenAngle = MaxOpenAngle;
        fPhotonsMaxinvmass   = Maxphotonmass;
        fPhotonsMaxNSigTPC   = nMaxNsigTPCDaught;
        fPhotonsMaxChi2      = nMaxPhotchi2;
        fInvPi0Mass = pi0mass; 
        fPi0Px = trackPi0.Px(); 
        fPi0Py = trackPi0.Py(); 
        fPi0Pz = trackPi0.Pz(); 
        fProtonPx = prot->Px(); 
        fProtonPy = prot->Py(); 
        fProtonPz = prot->Pz();
        fProtonX = trackxyz[0];
        fProtonY = trackxyz[1];
        fProtonZ = trackxyz[2]; 
        fProtonEta = prot->Eta();
        fProtonpropPx = trackProtonRot.Px(); 
        fProtonpropPy = trackProtonRot.Py(); 
        fProtonpropPz = trackProtonRot.Pz();
        fProtonDCAtoPVxy = DCAxy; 
        fProtonDCAtoPVz = DCAz; 
        fProtonPi0DCA = ProtPhotonDCA;
        fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
        fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
        fProtonNCluster = prot->GetTPCNcls();
        fProtonNITSCluster = prot->GetITSNcls();
        fProtonChi2 = prot->GetTPCchi2();  
        fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
        fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
        fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
        fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
        fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
        fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);

        fCaloPhotonELead = -999;
        AliAODCaloCells* CaloCells;
        if(fIsClusterPHOS) CaloCells = (AliAODCaloCells*)aodEvent->GetPHOSCells();
        if(fIsClusterEMCAL) CaloCells = (AliAODCaloCells*)aodEvent->GetEMCALCells();
        if(CaloCells){
          for(Int_t iCell=0; iCell<fClustNCells; iCell++){double cellamp = CaloCells->GetCellAmplitude(fCellsAbsId[iCell]); if(cellamp>fCaloPhotonELead) fCaloPhotonELead = cellamp;}
        }

        fEventhasSigmaCand = kTRUE;
        if(fSavePartCandPHOS) fSigmaPHOSCandTree->Fill();

        //Now apply some Selections to filter out potential Sigma Candidates.
        //If it passes the criteria, Sigma Proton pairs in SE and ME can be written out. 
        //Also all Protons in the Event can be written for offline Event-Mixing.
        if(!fFillPHOSPairTreeSE&&!fFillPHOSPairTreeME) continue;
        if(pi0mass<fMinCorrPi0MassPHOS) continue;
        if(pi0mass>fMaxCorrPi0MassPHOS) continue;
        if(SigmaPointingAngle>fMaxCorrSigmaPAPHOS) continue;
        if(TMath::Abs(DCAxy)<fMinCorrProtonDCAxy) continue;
        if(sigmaplusmass<fMinCorrSigmaMass) continue; 
        if(sigmaplusmass>fMaxCorrSigmaMass) continue;

        if(sigmaplusmass>1.17&&sigmaplusmass<1.21) fEventhasSigma = kTRUE;

        for(Int_t q=0; q<nProton; q++) { 
          AliAODTrack *pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
          if(!pairprot) continue;

          //If the Sigma made it here, check properties of all Protons in the Event and write them to a tree

          TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
          TVector3 deltapvec=sigmamomentum-protonmomentum;
          Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
          Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
          Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
          Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
          Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
          fSigmaProtonkstar = 1000*kstar;

          if(fSigmaProtonkstar>fMaxCorrkstar) continue;

          deltapvec=trackSigmaplusRot.Vect()-protonmomentum;
          SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentumrot.Mag()*sigmamomentumrot.Mag());
          ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
          qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
          vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
          kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
          fSigmaProtonpropkstar = 1000*kstar;

          fPairProtonIsMC = kFALSE;
          fPairProtonIsPrimary = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* PairProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(pairprot->GetLabel())));
            if(PairProtonPart){
              if(TMath::Abs(PairProtonPart->GetPdgCode())==2212){
                fPairProtonIsMC = kTRUE;
                if(PairProtonPart->IsPhysicalPrimary()) fPairProtonIsPrimary = kTRUE;
              }//MC Particle is a Proton
            }//MC Particle exists 
          }//End of isMonteCarlo

          pairprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

          fPairProtonPx = pairprot->Px();        
          fPairProtonPy = pairprot->Py();        
          fPairProtonPz = pairprot->Pz();
          fPairProtonP = pairprot->P();
          fPairProtonEta = pairprot->Eta();
          fPairProtonCharge = pairprot->Charge();        
          fPairProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kProton);        
          fPairProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kProton);
          fPairProtNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kPion);
          fPairProtNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kKaon);
          fPairProtNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(pairprot,AliPID::kElectron);
          fPairProtNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kPion);
          fPairProtNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kKaon);
          fPairProtNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(pairprot,AliPID::kElectron);
          fPairProtonChi2 = pairprot->GetTPCchi2();    
          fPairProtonCluster = pairprot->GetTPCNcls();
          fPairProtonITSCluster = pairprot->GetITSNcls();
          fPairProtonID = pairprot->GetID();
          fPairProtonStatus = pairprot->GetStatus();
          fPairProtonFilterMap = pairprot->GetFilterMap();

          //Force TOF PID for Pair Proton if requested
          if(pairprot->P()>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

          if(fPairProtonDCAtoPVxy>fMaxCorrPairProtonDCAxy) continue;
          if(fPairProtonDCAtoPVz>fMaxCorrPairProtonDCAz) continue;

          if(q!=k&&fFillPHOSPairTreeSE) fSigmaPairTreePHOSSE->Fill();
        }        

        //Continue here if no Event Mixing is requested
        if(!fFillPHOSPairTreeME) continue;

        //Get Pool from Pool Manager for given RefMult and Z Vertex values
        AliEventPool* Evpool = 0x0;
	      if(fEvPoolMgr2&&fUseAbsZCorr)  Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
	      if(fEvPoolMgr2&&!fUseAbsZCorr) Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
		    if(!Evpool){ AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ)); continue;}
        if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << " is empty!\n";*/ continue;}

        //Get Number of Events in Pool
    	  Int_t nMixEvents = Evpool->GetCurrentNEvents();
        if(fUseAbsZCorr) FillHistogram("fHistPairNMixedEventsPHOS",fRefMultComb08,TMath::Abs(primaryVtxPosZ),nMixEvents);
        else FillHistogram("fHistPairNMixedEventsPHOS",fRefMultComb08,primaryVtxPosZ,nMixEvents);

        //Now Loop over the mixed Events
			  for (Int_t iMixEvent = 0; iMixEvent < nMixEvents; iMixEvent++){
          //Retrieve Array of Protons for each mixed event 
				  TObjArray* MixedProtons = (TObjArray*)Evpool->GetEvent(iMixEvent);
      		Int_t nMixProtons = MixedProtons->GetEntriesFast();  			  
          for(Int_t iMixProton = 0; iMixProton < nMixProtons; iMixProton++){

            AliAODTrackcorrelation *mixprot = (AliAODTrackcorrelation*)MixedProtons->At(iMixProton);
            if(!mixprot) continue;

            TVector3 protonmomentum(mixprot->Px(),mixprot->Py(),mixprot->Pz());
            TVector3 deltapvec=sigmamomentum-protonmomentum;
            Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
            Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
            fSigmaProtonkstar = 1000*kstar;

            if(fSigmaProtonkstar>fMaxCorrkstar) continue;

            deltapvec=trackSigmaplusRot.Vect()-protonmomentum;
            SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentumrot.Mag()*sigmamomentumrot.Mag());
            ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
            fSigmaProtonpropkstar = 1000*kstar;

            fPairProtonIsMC = kFALSE;
            fPairProtonIsPrimary = kFALSE;

            mixprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

            fPairProtonPx = mixprot->Px();        
            fPairProtonPy = mixprot->Py();        
            fPairProtonPz = mixprot->Pz();
            fPairProtonCharge = mixprot->Charge();        
            fPairProtonNSigTPC = mixprot->NumberOfSigmasTPCProton();        
            fPairProtonNSigTOF = mixprot->NumberOfSigmasTOFProton();
            fPairProtNSigTPCKaon = mixprot->NumberOfSigmasTPCKaon();        
            fPairProtNSigTOFKaon = mixprot->NumberOfSigmasTOFKaon();
            fPairProtNSigTPCPion = mixprot->NumberOfSigmasTPCPion();        
            fPairProtNSigTOFPion = mixprot->NumberOfSigmasTOFPion();
            fPairProtonChi2 = mixprot->GetTPCchi2();    
            fPairProtonCluster = mixprot->GetTPCNcls();
            fPairProtonITSCluster = mixprot->GetITSNcls();
            fPairProtonID = mixprot->GetID();
            fPairProtonStatus = mixprot->GetStatus();
            fPairProtonFilterMap = mixprot->GetFilterMap();

            TVector3 mprt(fPairProtonPx,fPairProtonPy,fProtonPz);
            fPairProtonP = mprt.Mag();
            fProtonEta = mprt.Eta();

            //Force TOF PID for Pair Proton if requested
            if(fPairProtonP>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

            if(fPairProtonDCAtoPVxy>fMaxCorrPairProtonDCAxy) continue;
            if(fPairProtonDCAtoPVz>fMaxCorrPairProtonDCAz) continue;

            fSigmaPairTreePHOSME->Fill();

          }//End of Loop over Mixed Protons
        }//End of Loop of Mixed Events
      } //End of Calo Photon Loop

      //Do Event Mixing for the Background if requested
      if(fSavePHOSMixedBackground){

        Int_t nMixEvents = 0;

        //Get Pool from Pool Manager for given RefMult and Z Vertex values
        AliEventPool* Evpool = 0x0;
	      if(fEvPoolMgr3&&fUseAbsZ)  Evpool = fEvPoolMgr3->GetEventPool((Int_t)fRefMultComb08, (Double_t)TMath::Abs(primaryVtxPosZ));
	      if(fEvPoolMgr3&&!fUseAbsZ) Evpool = fEvPoolMgr3->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
		    if(!Evpool){AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ)); nMixEvents = -1;}
        else if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << " is empty!\n";*/ nMixEvents = -1;}

        //Get Number of Events in Pool. Number can be reduced with Setter Function to reduce Tree Size.
      	if(nMixEvents!=-1) nMixEvents = Evpool->GetCurrentNEvents();
        else nMixEvents = 0;
        if(fUseAbsZ) FillHistogram("fHistBkgNMixedEventsPHOS",fRefMultComb08,TMath::Abs(primaryVtxPosZ),nMixEvents);
        else FillHistogram("fHistBkgNMixedEventsPHOS",fRefMultComb08,primaryVtxPosZ,nMixEvents);

        //Now Loop over the mixed Events
  			for (Int_t iMixEvent = 0; iMixEvent < nMixEvents; iMixEvent++){

          //Retrieve Array of Protons for each mixed event 
  				TObjArray* MixedClusters = (TObjArray*)Evpool->GetEvent(iMixEvent);
        	Int_t nMixClusters = MixedClusters->GetEntriesFast();  			  
          for(Int_t iMixCluster = 0; iMixCluster < nMixClusters; iMixCluster++){

            AliAODClusterreduced *mixCluster = (AliAODClusterreduced*)MixedClusters->At(iMixCluster);
            if(!mixCluster) continue;

            //Now the hole Sigma reconstruction is repeated with Clusters from a different Event.

            Float_t clustE = mixCluster->E();
            Float_t clustpos[3];
            mixCluster->GetPosition(clustpos);

            //Calculating beta of cluster 
            ClusterPhoton.SetXYZ(clustpos[0],clustpos[1],clustpos[2]);
            Double_t FlightDist = ClusterPhoton.Mag();
            Double_t vTOF = 0;
            if(mixCluster->GetTOF()!=0) vTOF = FlightDist/(1e12*mixCluster->GetTOF());
            Float_t clustbeta = vTOF/c;

            //Set Calo Photon. Best assumtion: Origin = Secondary Vertex
            ClusterPhoton.SetXYZ(clustpos[0]-KFSigmaPlus.GetX(),clustpos[1]-KFSigmaPlus.GetY(),clustpos[2]-KFSigmaPlus.GetZ());
            //Calculate pion mass from available information  
            Double_t pi0mass = TMath::Sqrt(2*clustE*ConvPhoton.Mag()*(1-TMath::Cos(ConvPhoton.Angle(ClusterPhoton))));

            //Normalize Calo Photon and scale to cluster energy
            ClusterPhoton*=clustE/ClusterPhoton.Mag();
            trackSigmaplus = trackPhoton1 + trackPhoton2 + trackProton;
            Float_t sigmamassuncorr = trackSigmaplus.M();  //Save uncorrected mass

            if(pi0mass<fMinPi0MassPHOS || pi0mass>fMaxPi0MassPHOS) continue;  //Coarse mass cut to reduce combinatorics!

            //Reset Calo Photon. Best assumtion: Origin = Secondary Vertex
            ClusterPhoton.SetXYZ(clustpos[0]-KFSigmaPlus.GetX(),clustpos[1]-KFSigmaPlus.GetY(),clustpos[2]-KFSigmaPlus.GetZ());
            //Calculate cluster energy using nominal pion mass
            Double_t ClustECorr = cPi0Mass*cPi0Mass/(2*ConvPhoton.Mag()*(1-TMath::Cos(ConvPhoton.Angle(ClusterPhoton))));       
            //Normalize Calo Photon and scale to corrected energy
            ClusterPhoton*=ClustECorr/ClusterPhoton.Mag();
            //Set TLorentzvector for further calculations  
            trackPhoton2.SetXYZM(ClusterPhoton.X(),ClusterPhoton.Y(),ClusterPhoton.Z(),0);            

            trackPi0 = trackPhoton1 + trackPhoton2;  
            trackSigmaplus = trackPi0 + trackProton;

            //Calculate sigma mass
            Float_t sigmaplusmass = trackSigmaplus.M();
            if(sigmaplusmass>fMaxSigmaMassPHOS) continue;   //Limit the mass range to reduce tree size

            TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
            Float_t Rapidity = sig.Rapidity();

            TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
            TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
            Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);

            TLorentzVector trackSigmaplusRot = trackPi0 + trackProtonRot;
            TVector3 sigmamomentumrot(trackSigmaplusRot.Px(),trackSigmaplusRot.Py(),trackSigmaplusRot.Pz());
            TVector3 sigmavertexrot(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
            Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
            Lprop = sigmavertex.Mag();
            Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
            if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
            sigmavertexrot.RotateZ(Alpha);          
            Float_t AntiSigmaPointingAngle = trackProtonRot.Angle(sigmavertexrot); //Fake Pointing Angle

            Lprop = protonpath.Mag();
            Rcurve = sigmamomentumrot.Mag()*1000/(0.2998*TMath::Abs(Bz));
            Alpha = propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
            if(Lprop/(2*Rcurve)<1) Alpha = 2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
            trackSigmaplusRot.RotateZ(Alpha);
            Float_t SigmaPointingAnglerot = sigmamomentumrot.Angle(sigmavertexrot); 
            Float_t sigmaplusmassrot = trackSigmaplusRot.M();

            FillHistogram("fHistInvSigmaMassPHOSmix",sigmaplusmass);
            FillHistogram("fHistSigmaYPHOSmix",Rapidity);
            FillHistogram("fHistSigmaPAPHOSmix",SigmaPointingAngle);
            FillHistogram("fHistSigmaAntiPAPHOSmix",AntiSigmaPointingAngle);

            //Coarse topological cuts to reduce the tree size
            if(SigmaDisttoPV<fMinSigmaDCAtoPVPHOS||SigmaDisttoPV>fMaxSigmaDCAtoPVPHOS) continue;
            if(TMath::Abs(Rapidity)>fMaxSigmaYPHOS) continue;  //Cut on the rapidity
            if(SigmaPointingAngle>fMaxSigmaPAPHOS) continue;  
            if(EventTriggers&65536||EventTriggers&8){
              if(!isMonteCarlo&&SigmaPointingAngle>fMaxSigmaPAPHOSHM) continue;  //Apply stricter cut in case of high multiplicity          
            }
            if(AntiSigmaPointingAngle<fMinSigmaAntiPAPHOS) continue; 
            if(ProtPhotonDCA>fMaxProtPhotDCA) continue;

            // Fill the ME Sigma Trees
            fIsClusterEMCAL = mixCluster->IsEMCAL();
            fIsClusterPHOS = mixCluster->IsPHOS();
            fInvSigMass = sigmaplusmass;
            fInvSigpropMass = sigmaplusmassrot;
            fInvSigMassUncorr = sigmamassuncorr;
            fSigY = trackSigmaplus.Rapidity();               
            fSigYprop = trackSigmaplusRot.Rapidity();                
            fSigPA = SigmaPointingAngle; 
            fSigPAprop = SigmaPointingAnglerot; 
            fSigAntiPA = AntiSigmaPointingAngle; 
            fSigPx = trackSigmaplus.Px(); 
            fSigPy = trackSigmaplus.Py(); 
            fSigPz = trackSigmaplus.Pz(); 
            fSigPt = trackSigmaplus.Pt();
            fSigPxprop = trackSigmaplusRot.Px(); 
            fSigPyprop = trackSigmaplusRot.Py(); 
            fSigPzprop = trackSigmaplusRot.Pz();  
            fPhoton2Px = trackPhoton2.Px();
            fPhoton2Py = trackPhoton2.Py();
            fPhoton2Pz = trackPhoton2.Pz();
            fCaloPhotonX = clustpos[0];
            fCaloPhotonY = clustpos[1];
            fCaloPhotonZ = clustpos[2];
            fCaloPhotonE = clustE;
            fCaloPhotonEcorr = ClustECorr;
            fClustDispersion = mixCluster->GetDispersion();
            fClustM20 = mixCluster->GetM20();
            fClustM02 = mixCluster->GetM02();
            fClustNTracksMatched = mixCluster->GetNTracksMatched();
            fClustTrackDx = mixCluster->GetTrackDx();
            fClustTrackDz = mixCluster->GetTrackDz();
            fClustTrackD = TMath::Sqrt(fClustTrackDx*fClustTrackDx+fClustTrackDz*fClustTrackDz);
            fClustTOF = 1000000000*mixCluster->GetTOF();
            fClustBeta = clustbeta;
            fClustNCells = mixCluster->GetNCells();
            fClustDisttoBC = mixCluster->GetDistanceToBadChannel();
            for(Int_t iCell=0; iCell<20; iCell++){fCellsAbsId[iCell] = (Short_t)mixCluster->GetCellAbsId(iCell);}
            fInvPi0Mass = pi0mass; 
            fPi0Px = trackPi0.Px(); 
            fPi0Py = trackPi0.Py(); 
            fPi0Pz = trackPi0.Pz(); 
            fCaloPhotonELead = mixCluster->ELead();
            fSigmaPHOSMEBkgTree->Fill();
  
          }//End of Mixed Cluster loop
        }//End of Mixed Event loop
      }//End of Mixed Event Background
    } //End of Proton Loop
  } //End of Conversion Photon Loop

  return;  

} //End of ReconstructParticlesPHOS()

//________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesCalc() {

  KFParticle::SetField(Bz);

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TVector3 ConvPhoton, ClusterPhoton; 
  TLorentzVector trackPhoton1, trackPhoton2_1, trackPhoton2_2, trackPhoton2, trackProton, trackSigmaplus;
  KFParticle KFElectron, KFPositron, KFProton; 

  const Int_t nConvPhoton = fConvPhotonArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nConvPhoton; i++) {

    AliAODv0 *v0_1 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(i));
    if(!v0_1) continue;

    // Get daughter tracks      
    AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(0));
    AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(v0_1->GetDaughter(1));

    if(!track1 || !track2) {
      AliWarning("ERROR: Could not retrieve all AOD tracks in Pi0 reconstruction!");
      continue;
    }

    //KF calculations
    // Set up KFParticle
    trackparams[0] = v0_1->DecayVertexV0X();
    trackparams[1] = v0_1->DecayVertexV0Y();
    trackparams[2] = v0_1->DecayVertexV0Z();
    trackparams[3] = v0_1->MomPosX();
    trackparams[4] = v0_1->MomPosY();
    trackparams[5] = v0_1->MomPosZ();
    if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
    else track2->GetCovarianceXYZPxPyPz(covMatrix);
    KFPositron.Create(trackparams,covMatrix,1,cElectronMass);

    // Repeat for all other particles
    trackparams[0] = v0_1->DecayVertexV0X();
    trackparams[1] = v0_1->DecayVertexV0Y();
    trackparams[2] = v0_1->DecayVertexV0Z();
    trackparams[3] = v0_1->MomNegX();
    trackparams[4] = v0_1->MomNegY();
    trackparams[5] = v0_1->MomNegZ();
    if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
    else track2->GetCovarianceXYZPxPyPz(covMatrix);
    KFElectron.Create(trackparams,covMatrix,-1,cElectronMass);

    KFParticleCD KFPhotonCD; //Check Daughters to avoid floating point exceptions. See .h-file
    KFPhotonCD.AddDaughter(KFElectron);
    if(!KFPhotonCD.CheckDaughter(KFPositron)) continue;

    //Reconstruct the Photons with two Methods: Standard and special gamma reconstruction
    KFParticle KFPhoton(KFElectron,KFPositron);

    //Transport Photon to Conversion Point
    KFPhoton.TransportToDecayVertex();

    //Calculating DCA of Photon to PV 
    TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);                            //Prim. Vertex
    TVector3 CV1(v0_1->DecayVertexV0X(),v0_1->DecayVertexV0Y(),v0_1->DecayVertexV0Z()); //Conv. Vertices
    TVector3 p1(v0_1->Px(),v0_1->Py(),v0_1->Pz());                     //Momentum vectors of the photons
    Double_t DCAPV1 = (p1.Cross((CV1-PV))).Mag()/p1.Mag();                        //DCA to PV of Photons

    // Save Photon quality
    Int_t nMinTPCClustDaught = track1->GetTPCNcls();
    if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();

    Int_t nMinITSClustDaught = track1->GetITSNcls();
    if(track2->GetITSNcls()<nMinITSClustDaught) nMinITSClustDaught = track2->GetITSNcls();

    Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));

    Double_t MaxAlpha = TMath::Abs(v0_1->AlphaV0());
    Double_t MaxQt = TMath::Abs(v0_1->PtArmV0());
    Double_t MaxOpenAngle = v0_1->OpenAngleV0();

    TLorentzVector tl1, tl2, tlp;
    tl1.SetXYZM(v0_1->MomNegX(),v0_1->MomNegY(),v0_1->MomNegZ(),cElectronMass);
    tl2.SetXYZM(v0_1->MomPosX(),v0_1->MomPosY(),v0_1->MomPosZ(),cElectronMass);
    tlp=tl1+tl2;
    Double_t Maxphotonmass = tlp.M();

    Double_t nMaxPhotchi2 = v0_1->Chi2V0();

    Double_t MaxDaughtEta = TMath::Abs(track1->Eta());
    if(TMath::Abs(track2->Eta())>MaxDaughtEta) MaxDaughtEta = TMath::Abs(track2->Eta());

    Double_t MaxDeltaTheta = TMath::Abs(track1->Theta()-track2->Theta());

    //Set TLorentzvector for further calculations  
    trackPhoton1.SetXYZM(v0_1->Px(),v0_1->Py(),v0_1->Pz(),0);    
    //Set also TVector3 (better suited for some applications)
    ConvPhoton.SetXYZ(v0_1->Px(),v0_1->Py(),v0_1->Pz());

    for(Int_t k=0; k<nProton; k++) {

      AliAODTrack *prot;
      prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(k));
      if(!prot) continue;

      //MC treatment
      Double_t MCPi0Vtx[3];
      Int_t    SigmaMCLabel=-1;
      Bool_t   isReallySigma = kFALSE;
      Bool_t   isDalitzDecay = kFALSE;
      Bool_t   isPrimary = kFALSE;
      TLorentzVector MCSigmaMom, MCCaloMom; 

      if(isMonteCarlo){
        AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
        AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
        AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));

        Int_t Pi0ID = -1;
        if(V01Daught1&&V01Daught2){
          if(V01Daught1->GetMother()!=-1&&TMath::Abs(V01Daught1->GetMother())==TMath::Abs(V01Daught2->GetMother())){
            AliAODMCParticle* V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
            if(V0Part){
              if(V0Part->GetPdgCode()==22){
                Pi0ID = V0Part->GetMother();
              }//V0 is photon
              if(V0Part->GetPdgCode()==111){
                Pi0ID = V0Part->Label();
                isDalitzDecay = kTRUE;
              }//V0 is from dalitz decay
              Int_t SigmaID = -1;
              AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0ID)));
              if(Pi0Part){
                if(Pi0Part->GetPdgCode()==111){
                  MCPi0Vtx[0]=Pi0Part->Xv();
                  MCPi0Vtx[1]=Pi0Part->Yv();
                  MCPi0Vtx[2]=Pi0Part->Zv();
                  SigmaID = Pi0Part->GetMother();
                  MCCaloMom.SetXYZM(MCPi0Vtx[0]-V0Part->Px(),MCPi0Vtx[1]-V0Part->Py(),MCPi0Vtx[2]-V0Part->Pz(),cPi0Mass);
                }//Particle is Pi0
              }//Particle exists
              if(ProtonPart){
                if(TMath::Abs(ProtonPart->GetPdgCode())==2212&&ProtonPart->GetMother()!=-1&&ProtonPart->GetMother()==SigmaID){
                  AliAODMCParticle* SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(SigmaID)));
                  if(SigmaPart){
                    if(TMath::Abs(SigmaPart->GetPdgCode())==3222){
                      SigmaMCLabel = SigmaID;
                      isReallySigma = kTRUE;
                      if(SigmaPart->IsPhysicalPrimary()) isPrimary = kTRUE; 
                      MCSigmaMom.SetXYZM(SigmaPart->Px(),SigmaPart->Py(),SigmaPart->Pz(),cSigmaMass);
                    }//Particle is Sigma+
                  }//Sigma exists
                }//Particle is Proton and has common mother with pion 
              }//Proton exists
            }//V0 Particle exists
          }//V0 Daughters have common mother
        }//V0 Daughters exist
      }//End of MC

      FillHistogram("fHistCalcStatistics",1);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",1);

      //Reject autocorrelations
      if(fCheckProtonV0IDs&&(prot->GetID()==track1->GetID()||prot->GetID()==track2->GetID())) continue;

      FillHistogram("fHistCalcStatistics",2);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",2);

      Float_t  DCAxy = -999., DCAz = -999.;
      prot->GetImpactParameters(DCAxy,DCAz);

      if(TMath::Abs(DCAxy)<fMinProtonDCAxyPHOS) continue; //Coarse cut on DCA to PV to reduce the tree size
      if(TMath::Abs(DCAz)<fMinProtonDCAzPHOS) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxyPHOS&&DCAxy!=-999) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAzPHOS&&DCAz!=-999) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxyPHOS&&fRequireDCACutPHOS) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAzPHOS&&fRequireDCACutPHOS) continue;

      FillHistogram("fHistCalcStatistics",3);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",3);

      //Set TLorentzvector for further calculations  
      trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);

      prot->GetXYZ(trackxyz);      
      prot->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      prot->GetCovarianceXYZPxPyPz(covMatrix);
      if(prot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
      else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

      //Reconstruct pseudo sigma plus from proton and photon to access the decay vertex.
      KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
      KFSigmaPlusCD.AddDaughter(KFProton);
      if(!KFSigmaPlusCD.CheckDaughter(KFPhoton)) continue;

      FillHistogram("fHistCalcStatistics",4);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",4);

      KFParticle KFSigmaPlus(KFProton,KFPhoton);
      KFSigmaPlus.TransportToDecayVertex();

      //Calculate distance to primary vertex and reject candidates which are too close to the PV
      Double_t SigmaDisttoPV = TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY)+(KFSigmaPlus.GetZ()-primaryVtxPosZ)*(KFSigmaPlus.GetZ()-primaryVtxPosZ));

      FillHistogram("fHistSigmaDCAtoPVCalc",SigmaDisttoPV);
      if(isReallySigma) FillHistogram("fHistSigmaDCAtoPVCalcMC",SigmaDisttoPV);

      if(SigmaDisttoPV<fMinSigmaDCAtoPVPHOS||SigmaDisttoPV>fMaxSigmaDCAtoPVPHOS) continue;
      FillHistogram("fHistCalcStatistics",5);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",5);

      //Calculate distance between Photon and Proton
      Double_t ProtPhotonDCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPhoton));

      //Fill some histograms
      FillHistogram("fHistProtPhotonDCACalc",ProtPhotonDCA);
      if(isReallySigma) FillHistogram("fHistProtPhotonDCACalcMC",ProtPhotonDCA);

      if(ProtPhotonDCA>fMaxProtPhotDCA) continue;
      FillHistogram("fHistCalcStatistics",6);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",6);

      //Calculating DCA of Photon and Proton to Decay Vertex 
      Double_t DCASV1, DCATrack;
      if(fSaveAdditionalBranches){
        TVector3 SV(KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()); //Sec. Vertex
        DCASV1 = (p1.Cross((CV1-SV))).Mag()/p1.Mag();                 //DCA to SV of Photons
        Double_t Trackpnt[3];
        prot->GetXYZatR(TMath::Sqrt((KFSigmaPlus.GetX()*KFSigmaPlus.GetX())+(KFSigmaPlus.GetY()*KFSigmaPlus.GetY())),Bz,Trackpnt,0);
        TVector3 TrackpntVec(Trackpnt[0],Trackpnt[1],Trackpnt[2]);  
        DCATrack = (TrackpntVec-SV).Mag();
      }

      //Clone proton and propagate it to the secondary vertex
      TLorentzVector trackProtonRot = trackProton;
      TVector3 protonpath(KFSigmaPlus.GetX()-trackxyz[0],KFSigmaPlus.GetY()-trackxyz[1],KFSigmaPlus.GetZ()-trackxyz[2]);
      Double_t propdir = 1; 
      if(TMath::Sqrt((KFSigmaPlus.GetX()-primaryVtxPosX)*(KFSigmaPlus.GetX()-primaryVtxPosX)+(KFSigmaPlus.GetY()-primaryVtxPosY)*(KFSigmaPlus.GetY()-primaryVtxPosY))
      <TMath::Sqrt((trackxyz[0]-primaryVtxPosX)*(trackxyz[0]-primaryVtxPosX)+(trackxyz[1]-primaryVtxPosY)*(trackxyz[1]-primaryVtxPosY))) propdir = -1;
      Double_t qprot = 1; if(prot->Charge()<0) qprot=-1;
      Double_t Lprop = protonpath.Mag();
      Double_t Rcurve = trackProton.P()*1000/(0.2998*TMath::Abs(Bz));
      Double_t Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
      if(Lprop/(2*Rcurve)<1) Alpha = -2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
      trackProtonRot.RotateZ(Alpha);

      TVector3 sigmavertexrot(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
      Rcurve = trackProton.P()*1000/(0.2998*TMath::Abs(Bz));
      Lprop = sigmavertexrot.Mag();
      Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
      if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
      sigmavertexrot.RotateZ(Alpha);          

      Float_t AntiSigmaPointingAngle = trackProtonRot.Angle(sigmavertexrot); //Fake Pointing Angle
      FillHistogram("fHistSigmaAntiPACalc",AntiSigmaPointingAngle);
      if(isReallySigma) FillHistogram("fHistSigmaAntiPACalcMC",AntiSigmaPointingAngle);

      if(AntiSigmaPointingAngle<fMinSigmaAntiPACalc) continue; 
      if(AntiSigmaPointingAngle>fMaxSigmaAntiPACalc) continue; 
      FillHistogram("fHistCalcStatistics",7);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",7);

      TVector3 e1, e2, e3, ProtonOBase, Photon1OBase, Photon2OBase;
      Double_t xBase, yBase, zBase, bBase, cBase, aBase, EqParA, EqParB, EqParC, EqRes1, EqRes2, EqRes2Arg = -1;

      //Iterate to fix NaNs
      Int_t nIter = -1;
      Double_t ThetaVar=0, PhiVar=0;

      while(EqRes2Arg<0&&nIter<fMaxIter){
            //Obtain orthonormal basis
            e1 = sigmavertexrot;
            if(nIter!=-1){
              if(!fThetaFunc||!fPhiFunc){cout << "Resolution functions do not exist!\n"; break;}
              ThetaVar = fThetaFunc->GetRandom(-fThetaRange,fThetaRange);
              PhiVar = fPhiFunc->GetRandom(-fPhiRange,fPhiRange);
              e1.SetMagThetaPhi(1,sigmavertexrot.Theta()+ThetaVar,sigmavertexrot.Phi()+PhiVar);        
            }
            e1*=1./e1.Mag();
            e2 = e1.Cross(trackProtonRot.Vect()); e2*=1./e2.Mag();        
            e3 = e1.Cross(e2); e3*=1./e3.Mag();

            xBase = trackProtonRot.Px(); yBase = trackProtonRot.Py(); zBase = trackProtonRot.Pz();
            bBase = (((e1.Y()*xBase/e1.X())-yBase)/((e1.Y()*e3.X()/e1.X())-e3.Y())-((e1.Z()*xBase/e1.X())-zBase)/((e1.Z()*e3.X()/e1.X())-e3.Z()))/((e2.Z()-(e1.Z()*e2.X()/e1.X()))/((e1.Z()*e3.X()/e1.X())-e3.Z())-(e2.Y()-(e1.Y()*e2.X()/e1.X()))/((e1.Y()*e3.X()/e1.X())-e3.Y()));
            cBase = (((e1.Y()*xBase/e1.X())-yBase) + (e2.Y()-(e1.Y()*e2.X()/e1.X()))*bBase) / ((e1.Y()*e3.X()/e1.X())-e3.Y());
            aBase = (xBase-bBase*e2.X()-cBase*e3.X())/e1.X();
            ProtonOBase.SetXYZ(aBase,bBase,cBase);

            xBase = ConvPhoton.X(); yBase = ConvPhoton.Y(); zBase = ConvPhoton.Z();
            bBase = (((e1.Y()*xBase/e1.X())-yBase)/((e1.Y()*e3.X()/e1.X())-e3.Y())-((e1.Z()*xBase/e1.X())-zBase)/((e1.Z()*e3.X()/e1.X())-e3.Z()))/((e2.Z()-(e1.Z()*e2.X()/e1.X()))/((e1.Z()*e3.X()/e1.X())-e3.Z())-(e2.Y()-(e1.Y()*e2.X()/e1.X()))/((e1.Y()*e3.X()/e1.X())-e3.Y()));
            cBase = (((e1.Y()*xBase/e1.X())-yBase) + (e2.Y()-(e1.Y()*e2.X()/e1.X()))*bBase) / ((e1.Y()*e3.X()/e1.X())-e3.Y());
            aBase = (xBase-bBase*e2.X()-cBase*e3.X())/e1.X();
            Photon1OBase.SetXYZ(aBase,bBase,cBase);

            Photon2OBase.SetY(-Photon1OBase.Y());
            Photon2OBase.SetZ(-Photon1OBase.Z()-ProtonOBase.Z());

            //Calculate missing momentum component
            EqParA = ((cPi0Mass*cPi0Mass)/2+Photon1OBase.Y()*Photon2OBase.Y()+Photon1OBase.Z()*Photon2OBase.Z())/Photon1OBase.Mag();
            EqParB = Photon1OBase.X()/Photon1OBase.Mag();
            EqParC = Photon2OBase.Y()*Photon2OBase.Y()+Photon2OBase.Z()*Photon2OBase.Z();
            EqRes1 = EqParA*EqParB/(1-(EqParB*EqParB));
            EqRes2Arg = (EqParA*EqParA*EqParB*EqParB-(EqParB*EqParB-1)*(EqParA*EqParA-EqParC));
            if(EqRes2Arg>=0) EqRes2 = TMath::Sqrt(EqRes2Arg)/(EqParB*EqParB-1);
            else EqRes2 = -999;
            nIter++;
      }      

      if(EqRes2Arg>=0) FillHistogram("fHistNIterCalc",nIter);
      else FillHistogram("fHistNIterCalc",nIter+1);
      if(isReallySigma){
          if(EqRes2Arg>=0) FillHistogram("fHistNIterCalcMC",nIter);
          else FillHistogram("fHistNIterCalcMC",nIter+1);
      }

      //Check for NaNs
      if(EqRes2Arg<0) continue;

      FillHistogram("fHistCalcStatistics",8);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",8);

      FillHistogram("fHistEqPar1Calc",EqRes1);
      if(isReallySigma){
          FillHistogram("fHistEqPar1CalcMC",EqRes1);
      }

      //Reject unphysical phase space
      if(EqRes1<0&&fRejectNegPar1) continue;

      FillHistogram("fHistCalcStatistics",9);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",9);

      Photon2OBase.SetX(EqRes1+EqRes2);
      ClusterPhoton = Photon2OBase.X()*e1 + Photon2OBase.Y()*e2 + Photon2OBase.Z()*e3;
      trackPhoton2_1.SetXYZM(ClusterPhoton.X(),ClusterPhoton.Y(),ClusterPhoton.Z(),0);

      Photon2OBase.SetX(EqRes1-EqRes2);
      ClusterPhoton = Photon2OBase.X()*e1 + Photon2OBase.Y()*e2 + Photon2OBase.Z()*e3;
      trackPhoton2_2.SetXYZM(ClusterPhoton.X(),ClusterPhoton.Y(),ClusterPhoton.Z(),0);

      //Check which sign is the physical one
      if(TMath::Abs((trackPhoton1+trackPhoton2_1+trackProtonRot).M()-cSigmaMass)<TMath::Abs((trackPhoton1+trackPhoton2_2+trackProtonRot).M()-cSigmaMass)) trackPhoton2 = trackPhoton2_1;
      else trackPhoton2 = trackPhoton2_2;

      trackSigmaplus = trackPhoton1 + trackPhoton2 + trackProtonRot;
      Float_t pi0mass = (trackPhoton1+trackPhoton2).M();
      Float_t sigmaplusmass = trackSigmaplus.M();

      TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
      Float_t Rapidity = sig.Rapidity();

      sigmavertexrot.SetXYZ(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
      Rcurve = trackSigmaplus.P()*1000/(0.2998*TMath::Abs(Bz));
      Lprop = sigmavertexrot.Mag();
      Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi()/2;
      if(Lprop/(2*Rcurve)<1) Alpha = -1*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
      sigmavertexrot.RotateZ(Alpha);          

      TLorentzVector trackSigmaplusRot = trackSigmaplus;
      Alpha = propdir*qprot*TMath::Sign(1,Bz)*TMath::Pi();
      if(Lprop/(2*Rcurve)<1) Alpha = 2*propdir*qprot*TMath::Sign(1,Bz)*TMath::ASin(Lprop/(2*Rcurve));
      trackSigmaplusRot.RotateZ(Alpha);

      //Fill some more histograms
      FillHistogram("fHistDeltaThetaCalc",ThetaVar);
      FillHistogram("fHistDeltaPhiCalc",PhiVar);
      FillHistogram("fHistInvSigmaMassCalc",sigmaplusmass);
      FillHistogram("fHistSigmaYCalc",Rapidity);
      if(isReallySigma){
          FillHistogram("fHistDeltaThetaCalcMC",ThetaVar);
          FillHistogram("fHistDeltaPhiCalcMC",PhiVar);
          FillHistogram("fHistInvSigmaMassCalcMC",sigmaplusmass);
          FillHistogram("fHistSigmaYCalcMC",Rapidity);
      }

      //Coarse topological cuts to reduce the tree size
      if(sigmaplusmass>fMaxSigmaMassCalc) continue; //Limit the mass range to reduce tree size

      FillHistogram("fHistCalcStatistics",10);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",10);

      if(TMath::Abs(Rapidity)>fMaxSigmaYCalc) continue; //Cut on the rapidity

      FillHistogram("fHistCalcStatistics",11);
      if(isReallySigma) FillHistogram("fHistCalcStatisticsMC",11);

      //MC for Correlations
      if(isMonteCarlo&&fSaveAddMCBranches){
        fTrackLabel = -999;
        fTrackPDGCode = -999;
        fTrackMotherID = -999;
        fTrackMotherPDGCode = -999;
        fTrackMotherMCPx = -999;
        fTrackMotherMCPy = -999;
        fTrackMotherMCPz = -999;
        fTrackMCPx = -999;
        fTrackMCPy = -999;
        fTrackMCPz = -999;
        fPhoton1Label = -1;
        fPhoton1PDGCode = -999;
        fPhoton1MotherID = -999;
        fPhoton1MotherPDGCode = -999;
        fPhoton1GMotherID = -999;
        fPhoton1GMotherPDGCode = -999;
        fPhoton1MCPx = -999;
        fPhoton1MCPy = -999;
        fPhoton1MCPz = -999;
        AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
        if(ProtonPart){
          fTrackLabel = ProtonPart->GetLabel();
          fTrackPDGCode = ProtonPart->GetPdgCode();
          fTrackMotherID = ProtonPart->GetMother();
          if(ProtonPart->IsPhysicalPrimary()) fTrackMotherID = -1;
          fTrackMCPx = ProtonPart->Px();
          fTrackMCPy = ProtonPart->Py();
          fTrackMCPz = ProtonPart->Pz();
          AliAODMCParticle* ProtonPartM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ProtonPart->GetMother())));
          if(ProtonPartM){
            fTrackMotherPDGCode = ProtonPartM->GetPdgCode();
            fTrackMotherMCPx = ProtonPartM->Px();
            fTrackMotherMCPy = ProtonPartM->Py();
            fTrackMotherMCPz = ProtonPartM->Pz();
          }//Proton Mother exists
        }//Track exists
        AliAODMCParticle* V01Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
        AliAODMCParticle* V01Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
        if(V01Daught1&&V01Daught2){
          if(V01Daught1->GetMother()!=-1&&V01Daught1->GetMother()==V01Daught2->GetMother()){
            AliAODMCParticle* V0Part1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V01Daught1->GetMother())));
            if(V0Part1){
              fPhoton1Label = V0Part1->GetLabel();
              fPhoton1PDGCode = V0Part1->GetPdgCode();
              fPhoton1MotherID = V0Part1->GetMother();
              if(V0Part1->IsPhysicalPrimary()) fPhoton1MotherID = -1;
              fPhoton1MCPx = V0Part1->Px();
              fPhoton1MCPy = V0Part1->Py();
              fPhoton1MCPz = V0Part1->Pz();
              AliAODMCParticle* V0Part1M = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
              if(V0Part1M){
                fPhoton1MotherPDGCode = V0Part1M->GetPdgCode();
                fPhoton1GMotherID = V0Part1M->GetMother();
                if(V0Part1M->IsPhysicalPrimary()) fPhoton1GMotherID = -1;
                AliAODMCParticle* V0Part1GM = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1M->GetMother())));
                if(V0Part1GM){
                  fPhoton1GMotherPDGCode = V0Part1GM->GetPdgCode();        
                }//V0 G Mother exists
              }//V0 Mother exits
            }//MC V0 exists
          }//Real V0
        }//V01s Daughters exist
      }//IsMonteCarlo

      // Fill the Sigma Candidate Trees
      fSigEventID = fGlobalEventID;
      fSigBField = Bz;
      fSigCentrality = Centrality;
      fSigRefMultComb08 = fRefMultComb08;
      fSigRunnumber = aodEvent->GetRunNumber();
      fSigTriggerMask = EventTriggers;
      fPrimVertX = primaryVtxPosX;
      fPrimVertY = primaryVtxPosY;
      fPrimVertZ = primaryVtxPosZ;
      fSigCharge = prot->Charge();
      fSigProtonID = prot->GetID();
      fSigProtonStatus = prot->GetStatus();
      fSigProtonFilterMap = prot->GetFilterMap();
      fProtonEta = prot->Eta();
      fProtonDCAtoPVxy = DCAxy;
      fProtonDCAtoPVz = DCAz;
      fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
      fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
      fProtonNCluster = prot->GetTPCNcls();
      fProtonNITSCluster = prot->GetITSNcls();
      fProtonChi2 = prot->GetTPCchi2();  
      fProtonPx = prot->Px();
      fProtonPy = prot->Py();
      fProtonPz = prot->Pz();
      fProtonpropPx = trackProtonRot.Px();
      fProtonpropPy = trackProtonRot.Py();
      fProtonpropPz = trackProtonRot.Pz();
      fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
      fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
      fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
      fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
      fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
      fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);
      fIsV01Onthefly = v0_1->GetOnFlyStatus();
      fPhotonsMaxalpha = MaxAlpha;
      fPhotonsMaxqt = MaxQt;
      fPhotonsMaxOpenAngle = MaxOpenAngle;
      fPhotonsMaxinvmass = Maxphotonmass;
      fPhotonsMaxNSigTPC = nMaxNsigTPCDaught;
      fPhotonsMaxChi2 = nMaxPhotchi2;
      fPhotonDaughtMaxEta = MaxDaughtEta;
      fPhotonsMaxDeltaTheta = MaxDeltaTheta;
      fPhoton1Radius = TMath::Sqrt(v0_1->DecayVertexV0X()*v0_1->DecayVertexV0X()+v0_1->DecayVertexV0Y()*v0_1->DecayVertexV0Y());
      fPhotonsMinCluster = nMinTPCClustDaught;
      fPhotonsMinITSCluster = nMinITSClustDaught;
      fPhoton1Px = v0_1->Px();
      fPhoton1Py = v0_1->Py();
      fPhoton1Pz = v0_1->Pz();
      fProtonPi0DCA = ProtPhotonDCA;
      fSigPxprop = trackSigmaplusRot.Px();
      fSigPyprop = trackSigmaplusRot.Py();
      fSigPzprop = trackSigmaplusRot.Pz();
      fSigPt = trackSigmaplusRot.Pt();
      fInvSigpropMass = sigmaplusmass;
      fSigYprop = trackSigmaplusRot.Rapidity();                
      fSigAntiPA = AntiSigmaPointingAngle;
      fSigFlightDist = SigmaDisttoPV;
      fSigDecayVertX = KFSigmaPlus.GetX();
      fSigDecayVertY = KFSigmaPlus.GetY();
      fSigDecayVertZ = KFSigmaPlus.GetZ();
      fPhoton2Px = trackPhoton2.Px();
      fPhoton2Py = trackPhoton2.Py();
      fPhoton2Pz = trackPhoton2.Pz();
      fNIter = nIter;
      fDeltaPhi = PhiVar;
      fDeltaTheta = ThetaVar;
      fIsMCSigma = isReallySigma;
      fIsMCDalitz = isDalitzDecay;
      fIsMCPrimary = (isReallySigma&&isPrimary);
      fSigMCLabel = SigmaMCLabel;
      fPrimVertXMC = primaryVtxPosXMC;
      fPrimVertYMC = primaryVtxPosYMC;
      fPrimVertZMC = primaryVtxPosZMC;
      fSigDecayVertXMC = MCPi0Vtx[0];
      fSigDecayVertYMC = MCPi0Vtx[1];
      fSigDecayVertZMC = MCPi0Vtx[2];
      fSigPxMC = MCSigmaMom.Px();
      fSigPyMC = MCSigmaMom.Py();
      fSigPzMC = MCSigmaMom.Pz();
      fCaloPhotonPxMC = MCCaloMom.Px();
      fCaloPhotonPyMC = MCCaloMom.Py();
      fCaloPhotonPzMC = MCCaloMom.Pz();
      fSigRefMultComb05 = fRefMultComb05;
      fSigRefMultComb10 = fRefMultComb10;
      fConvPhotonX = v0_1->DecayVertexV0X();
      fConvPhotonY = v0_1->DecayVertexV0Y();
      fConvPhotonZ = v0_1->DecayVertexV0Z();
      fPhoton1DCAPV = DCAPV1;
      fPhoton1DCASV = DCASV1;
      fTrackDCASV = DCATrack;
      fKFChi2 = KFSigmaPlus.GetChi2();
      fPhoton1CPA = v0_1->CosPointingAngle(primaryVtxPos);
      fProtonX = trackxyz[0];
      fProtonY = trackxyz[1];
      fProtonZ = trackxyz[2];

      if(fSavePartCandCalc) fSigmaCalcCandTree->Fill();
      fEventhasSigmaCand = kTRUE;

    } //End of Proton Loop
  } //End of Conversion Photon Loop
  return;  

} //End of ReconstructParticlesCalc()

//_____________________________________________________________________________

//////////Functions for Filling Histograms in TList "fOutputList"///////////////
//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillHistogram(const char * key, Double_t x) const{

  TH1F* Histo = dynamic_cast<TH1F*>(fOutputList->FindObject(key));
  
  if(!Histo){
  AliWarning(Form("Error: Histogram '%s' does not exist in TList 'fOutputList'!",key));      
  return;
  }    
  else Histo->Fill(x);
  return;
}
  
//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillHistogram(const char * key, Double_t x, Double_t y) const{

  TH2F* Histo = dynamic_cast<TH2F*>(fOutputList->FindObject(key));
  
  if(!Histo){
  AliWarning(Form("Error: Histogram '%s' does not exist in TList 'fOutputList'!",key));      
  return;
  }  
  else Histo->Fill(x,y);
  return;
}
  
//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillHistogram(const char * key, Double_t x, Double_t y, Double_t z) const{

  TH3F* Histo = dynamic_cast<TH3F*>(fOutputList->FindObject(key));
  
  if(!Histo){
  AliWarning(Form("Error: Histogram '%s' does not exist in TList 'fOutputList'!",key));      
  return;
  } 
  else Histo->Fill(x,y,z);
  return;
}  

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::Terminate(Option_t *)
{
    // terminate. Called at the END of the analysis (when all events are processed).
}

//_____________________________________________________________________________
