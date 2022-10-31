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

ClassImp(AliAnalysisTaskSigmaPlus) // classimp: necessary for root
ClassImp(AliAODTrackreduced)
ClassImp(AliAODTrackcorrelation)

AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus() : AliAnalysisTaskSE(), 
fOutputList(0), aodEvent(0x0), mcEvent(0x0), AODMCTrackArray(0x0), 
fPIDResponse(0), fEvPoolMgr(0x0), fEvPoolMgr2(0x0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaCandTreeExtra(0x0), fSigmaPairTree(0x0), fProtonTree(0x0), 
fSigmaPairTreeSE(0x0), fSigmaPairTreeME(0x0), fSigmaMEBackgroundTree(0x0), fSigmaCandTreerot(0x0),
cElectronMass(0), cProtonMass(0), cSigmaMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(-999), primaryVtxPosY(-999), primaryVtxPosZ(-999),
primaryVtxPosXMC(-999), primaryVtxPosYMC(-999), primaryVtxPosZMC(-999), 
nTracks(0), Centrality(0), EventTriggers(0),
fRefMultComb05(0), fRefMultComb08(0), fRefMultComb10(0), fGlobalEventID(0), fMeanProtonpt(-1),
fEventhasSigma(kFALSE), fEventhasProton(kFALSE), fDebug(kFALSE), fRemoveGenPileup(kFALSE),

fProcessProtons(kTRUE), 
fProcessMCParticles(kTRUE), 
fProcessV0s(kTRUE), 
fProcessAddPhoton(kFALSE), 
fProcessElectrons(kFALSE), 
fProcessReco(kTRUE), 
fProcessRecoOff(kFALSE), 
fProcessRecoOff2(kFALSE), 
fProcessOneGamma(kFALSE), 
fSavePartCand(kTRUE),
fSavePairs(kFALSE),
fSaveAllProtons(kFALSE), 
fFillredPairTreeSE(kTRUE),
fFillredPairTreeME(kTRUE),
fSaveMixedBackground(kFALSE),
fSaveRotateBackground(kFALSE),

fMaxVertexZ(10),

fEvPoolSize(10),
fEvTrackSize(10000),
fCentralityBins(29),
fMinCentBin(-5),
fMaxCentBin(140),
fZvtxBins(20),
fMinZBin(-10),
fMaxZBin(10),

fEvPoolSize2(20),
fEvTrackSize2(10000),
fCentralityBins2(29),
fMinCentBin2(-5),
fMaxCentBin2(140),
fZvtxBins2(20),
fMinZBin2(-10),
fMaxZBin2(10),
fPtBins2(1),
fMinPtBin2(0),
fMaxPtBin2(10),
fRequireSigma(kTRUE),
fRequireProton(kTRUE),

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

fMinDCADaughtPV(0.1),
fMaxDCADaught(100),

fMaxElecEta(0.9),    
fMinTPCClustElec(40),
fMaxNsigElecTPC(3),
fMinNsigHadronTPC(0),
fMaxNsigElecTOF(3),  
fMaxElecpt(5), 

fCleanAutoCorr(kTRUE),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.06),
fMaxSigmaY(0.9),  
fMaxSigmaMass(1.4),
fMinProtonDCAxy(0.005),
fMinProtonDCAz(-1),
fMaxProtonDCAxy(5),
fMaxProtonDCAz(9999),
fRequireDCACut(kFALSE),
flowkstar(0.3),
fverylowkstar(0.15),
fveryverylowkstar(0.05),

fMinCorrPi0Mass(0.1),    
fMaxCorrPi0Mass(0.16),    
fMaxCorrSigmaPA(0.06),     
fMinCorrSigmaMass(1.13),   
fMaxCorrSigmaMass(1.25),   
fMinCorrProtonDCAxy(0.005),
fMaxCorrPairProtonDCAxy(0.3),
fMaxCorrPairProtonDCAz(0.6),
fMaxCorrkstar(0.6),

fIsMCSigma(kFALSE),
fIsMCPrimary(kFALSE),
fIsGoodCandidate(kFALSE),
fIsV01fromFinder(kFALSE),
fIsV02fromFinder(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fHas4DiffIDs(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(0),
fSigMCLabel(0),
fSigProtonID(-999),
fSigProtonStatus(0),
fSigEventID(0),
fSigCentrality(-999),
fSigRefMultComb05(-999),
fSigRefMultComb08(-999),
fSigRefMultComb10(-999),
fSigBField(-999),
fInvSigMass(-999),
fInvSigpropMass(-999),
fSigPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fSigpropPx(-999),
fSigpropPy(-999),
fSigpropPz(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
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
fProtonPxatDCA(-999),
fProtonPyatDCA(-999),
fProtonPzatDCA(-999),
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
fnPair(-999),
fnPairlowkstar(-999),
fnPairverylowkstar(-999),
fnPairveryverylowkstar(-999),

fIsV0fromFinder(kTRUE),
fIsV0Onthefly(kTRUE),
fPhotonPx(999),        
fPhotonPy(999),        
fPhotonPz(999),
fPhotonRadius(-999),
fPhotonDCAPV(-999),
fExtPhotProtDCA(999),

fPairProtonIsMC(kFALSE),
fPairProtonIsPrimary(kFALSE),
fPairProtonPx(-999),
fPairProtonPy(-999),
fPairProtonPz(-999),
fPairProtonPxatDCA(-999),
fPairProtonPyatDCA(-999),
fPairProtonPzatDCA(-999),
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
fSigmaProtonkstar(-999),
fSigmaProtonpropkstar(-999)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus(const char* name) : AliAnalysisTaskSE(name),
fOutputList(0), aodEvent(0x0), mcEvent(0x0), AODMCTrackArray(0x0), 
fPIDResponse(0), fEvPoolMgr(0x0), fEvPoolMgr2(0x0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaCandTreeExtra(0x0), fSigmaPairTree(0x0), fProtonTree(0x0), 
fSigmaPairTreeSE(0x0), fSigmaPairTreeME(0x0), fSigmaMEBackgroundTree(0x0), fSigmaCandTreerot(0x0),
cElectronMass(0), cProtonMass(0), cSigmaMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(-999), primaryVtxPosY(-999), primaryVtxPosZ(-999),
primaryVtxPosXMC(-999), primaryVtxPosYMC(-999), primaryVtxPosZMC(-999), 
nTracks(0), Centrality(0), EventTriggers(0),
fRefMultComb05(0), fRefMultComb08(0), fRefMultComb10(0), fGlobalEventID(0), fMeanProtonpt(-1),
fEventhasSigma(kFALSE), fEventhasProton(kFALSE), fDebug(kFALSE), fRemoveGenPileup(kFALSE),

fProcessProtons(kTRUE), 
fProcessMCParticles(kTRUE), 
fProcessV0s(kTRUE), 
fProcessAddPhoton(kFALSE), 
fProcessElectrons(kFALSE), 
fProcessReco(kTRUE), 
fProcessRecoOff(kFALSE), 
fProcessRecoOff2(kFALSE), 
fProcessOneGamma(kFALSE),
fSavePartCand(kTRUE),
fSavePairs(kFALSE),
fSaveAllProtons(kFALSE), 
fFillredPairTreeSE(kTRUE),
fFillredPairTreeME(kTRUE),
fSaveMixedBackground(kFALSE),
fSaveRotateBackground(kFALSE),

fMaxVertexZ(10),

fEvPoolSize(10),
fEvTrackSize(10000),
fCentralityBins(29),
fMinCentBin(-5),
fMaxCentBin(140),
fZvtxBins(20),
fMinZBin(-10),
fMaxZBin(10),

fEvPoolSize2(20),
fEvTrackSize2(10000),
fCentralityBins2(29),
fMinCentBin2(-5),
fMaxCentBin2(140),
fZvtxBins2(20),
fMinZBin2(-10),
fMaxZBin2(10),
fPtBins2(1),
fMinPtBin2(0),
fMaxPtBin2(10),
fRequireSigma(kTRUE),
fRequireProton(kTRUE),

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

fMinDCADaughtPV(0.1),
fMaxDCADaught(100),

fMaxElecEta(0.9),    
fMinTPCClustElec(40),
fMaxNsigElecTPC(3),
fMinNsigHadronTPC(0),
fMaxNsigElecTOF(3),  
fMaxElecpt(5), 

fCleanAutoCorr(kTRUE),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.06),
fMaxSigmaY(0.9),  
fMaxSigmaMass(1.4),
fMinProtonDCAxy(0.005),
fMinProtonDCAz(-1),
fMaxProtonDCAxy(5),
fMaxProtonDCAz(9999),
fRequireDCACut(kFALSE),
flowkstar(0.3),
fverylowkstar(0.15),
fveryverylowkstar(0.05),

fMinCorrPi0Mass(0.1),    
fMaxCorrPi0Mass(0.16),    
fMaxCorrSigmaPA(0.06),     
fMinCorrSigmaMass(1.13),   
fMaxCorrSigmaMass(1.25),   
fMinCorrProtonDCAxy(0.005),
fMaxCorrPairProtonDCAxy(0.3),
fMaxCorrPairProtonDCAz(0.6),
fMaxCorrkstar(0.6),

fIsMCSigma(kFALSE),
fIsMCPrimary(kFALSE),
fIsGoodCandidate(kFALSE),
fIsV01fromFinder(kFALSE),
fIsV02fromFinder(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fHas4DiffIDs(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(0),
fSigMCLabel(0),
fSigProtonID(-999),
fSigProtonStatus(0),
fSigEventID(0),
fSigCentrality(-999),
fSigRefMultComb05(-999),
fSigRefMultComb08(-999),
fSigRefMultComb10(-999),
fSigBField(-999),
fInvSigMass(-999),
fInvSigpropMass(-999),
fSigPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fSigpropPx(-999),
fSigpropPy(-999),
fSigpropPz(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
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
fProtonPxatDCA(-999),
fProtonPyatDCA(-999),
fProtonPzatDCA(-999),
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
fnPair(-999),
fnPairlowkstar(-999),
fnPairverylowkstar(-999),
fnPairveryverylowkstar(-999),

fIsV0fromFinder(kTRUE),
fIsV0Onthefly(kTRUE),
fPhotonPx(999),        
fPhotonPy(999),        
fPhotonPz(999),
fPhotonRadius(-999),
fPhotonDCAPV(-999),
fExtPhotProtDCA(999),

fPairProtonIsMC(kFALSE),
fPairProtonIsPrimary(kFALSE),
fPairProtonPx(-999),
fPairProtonPy(-999),
fPairProtonPz(-999),
fPairProtonPxatDCA(-999),
fPairProtonPyatDCA(-999),
fPairProtonPzatDCA(-999),
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

    DefineOutput(2, TTree::Class());    //Additional Output: Sigma Candidate Tree
    DefineOutput(3, TTree::Class());    //Additional Output: Extra Sigma Candidate Tree
    DefineOutput(4, TTree::Class());    //Additional Output: Sigma Proton Pair Tree
    DefineOutput(5, TTree::Class());    //Additional Output: Proton Tree
    DefineOutput(6, TTree::Class());    //Additional Output: reduced Sigma Proton Pair Tree SE
    DefineOutput(7, TTree::Class());    //Additional Output: reduced Sigma Proton Pair Tree ME
    DefineOutput(8, TTree::Class());    //Additional Output: Tree of Sigma Mixed Event Background
    DefineOutput(9, TTree::Class());    //Additional Output: Tree of Sigma Mixed Event Background
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::~AliAnalysisTaskSigmaPlus()
{
    // destructor
    // delete objects from memory at the end of the task
    if(fOutputList) delete fOutputList;     

    if(fSigmaCandTree) delete fSigmaCandTree;

    if(fSigmaCandTreeExtra) delete fSigmaCandTreeExtra;

    if(fSigmaPairTree) delete fSigmaPairTree;

    if(fProtonTree) delete fProtonTree;

    if(fSigmaPairTreeSE) delete fSigmaPairTreeSE;

    if(fSigmaPairTreeME) delete fSigmaPairTreeME;

    if(fSigmaMEBackgroundTree) delete fSigmaMEBackgroundTree;

    if(fSigmaCandTreerot) delete fSigmaCandTreerot;

    if(fEvPoolMgr) delete fEvPoolMgr;

    if(fEvPoolMgr2) delete fEvPoolMgr2;
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlus::UserCreateOutputObjects()
{
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

    const Int_t nCentralityBins2 = fCentralityBins2;
    const Int_t nZvtxBins2 = fZvtxBins2;
    const Int_t nPtBins2 = fPtBins2;

    if(fFillredPairTreeME){
      //Create Pool Manager for Event Mixing
      Int_t poolSize2 = fEvPoolSize2;
      Int_t trackDepth2 = fEvTrackSize2;
      Double_t centBinz2[nCentralityBins2+1];
      for(Int_t b=0;b<=nCentralityBins2;b++){centBinz2[b]=fMinCentBin2+b*(fMaxCentBin2-fMinCentBin2)/nCentralityBins2;}
      Double_t vertexBinz2[nZvtxBins2+1];
      for(Int_t b=0;b<=nZvtxBins2;b++){vertexBinz2[b]=fMinZBin2+b*(fMaxZBin2-fMinZBin2)/nZvtxBins2;}
      Double_t ptBinz2[nPtBins2+1];
      for(Int_t b=0;b<=nPtBins2;b++){ptBinz2[b]=fMinPtBin2+b*(fMaxPtBin2-fMinPtBin2)/nPtBins2;}
      Double_t* centBins2 = centBinz2; 
      Double_t* vertexBins2 = vertexBinz2;
      Double_t* ptBins2 = vertexBinz2;
      fEvPoolMgr2 = new AliEventPoolManager(poolSize2, trackDepth2, nCentralityBins2, centBins2, nZvtxBins2, vertexBins2, nPtBins2, ptBins2);
    }

    // Create TTree of Sigma Candidates
    fSigmaCandTree = new TTree("fSigmaCandTree","Tree of Sigma Candidates");
    fSigmaCandTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaCandTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaCandTree->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaCandTree->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaCandTree->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaCandTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaCandTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaCandTree->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaCandTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaCandTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaCandTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaCandTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaCandTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaCandTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCandTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaCandTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaCandTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaCandTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCandTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaCandTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaCandTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaCandTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaCandTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaCandTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaCandTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaCandTree->Branch("fSigpropPx",&fSigpropPx,"fSigpropPx/F");
    fSigmaCandTree->Branch("fSigpropPy",&fSigpropPy,"fSigpropPy/F");
    fSigmaCandTree->Branch("fSigpropPz",&fSigpropPz,"fSigpropPz/F");
    fSigmaCandTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaCandTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaCandTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCandTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
    fSigmaCandTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
    fSigmaCandTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
    fSigmaCandTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaCandTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaCandTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaCandTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaCandTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaCandTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaCandTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaCandTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaCandTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaCandTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaCandTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaCandTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaCandTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaCandTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaCandTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaCandTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaCandTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaCandTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaCandTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaCandTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaCandTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
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
    fSigmaCandTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaCandTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaCandTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaCandTree->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaCandTree->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaCandTree->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaCandTree->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaCandTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCandTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCandTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaCandTree->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaCandTree->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaCandTree->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaCandTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaCandTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaCandTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaCandTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaCandTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaCandTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaCandTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaCandTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaCandTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaCandTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaCandTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaCandTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaCandTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaCandTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaCandTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");

    // Create TTree of Extra Sigma Candidates
    fSigmaCandTreeExtra = new TTree("fSigmaCandTreeExtra","Tree of Extra Sigma Candidates");
    fSigmaCandTreeExtra->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaCandTreeExtra->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaCandTreeExtra->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");
    fSigmaCandTreeExtra->Branch("fIsV0fromFinder",&fIsV0fromFinder,"fIsV0fromFinder/O");
    fSigmaCandTreeExtra->Branch("fIsV0Onthefly",&fIsV0Onthefly,"fIsV0Onthefly/O");
    fSigmaCandTreeExtra->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTreeExtra->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaCandTreeExtra->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaCandTreeExtra->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaCandTreeExtra->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaCandTreeExtra->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaCandTreeExtra->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCandTreeExtra->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaCandTreeExtra->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaCandTreeExtra->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaCandTreeExtra->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCandTreeExtra->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaCandTreeExtra->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaCandTreeExtra->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaCandTreeExtra->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaCandTreeExtra->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaCandTreeExtra->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaCandTreeExtra->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaCandTreeExtra->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaCandTreeExtra->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaCandTreeExtra->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaCandTreeExtra->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaCandTreeExtra->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaCandTreeExtra->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaCandTreeExtra->Branch("fPhotonPx",&fPhotonPx,"fPhotonPx/F");
    fSigmaCandTreeExtra->Branch("fPhotonPy",&fPhotonPy,"fPhotonPy/F");
    fSigmaCandTreeExtra->Branch("fPhotonPz",&fPhotonPz,"fPhotonPz/F");
    fSigmaCandTreeExtra->Branch("fPhotonRadius",&fPhotonRadius,"fPhotonRadius/F");
    fSigmaCandTreeExtra->Branch("fPhotonDCAPV",&fPhotonDCAPV,"fPhotonDCAPV/F");
    fSigmaCandTreeExtra->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");                           
    fSigmaCandTreeExtra->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");                           
    fSigmaCandTreeExtra->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");                         
    fSigmaCandTreeExtra->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");                      
    fSigmaCandTreeExtra->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");                             
    fSigmaCandTreeExtra->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F"); 
    fSigmaCandTreeExtra->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaCandTreeExtra->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");     
    fSigmaCandTreeExtra->Branch("fExtPhotProtDCA",&fExtPhotProtDCA,"fExtPhotProtDCA/F");
    fSigmaCandTreeExtra->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCandTreeExtra->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCandTreeExtra->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaCandTreeExtra->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaCandTreeExtra->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaCandTreeExtra->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaCandTreeExtra->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaCandTreeExtra->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaCandTreeExtra->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaCandTreeExtra->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTreeExtra->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaCandTreeExtra->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaCandTreeExtra->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaCandTreeExtra->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaCandTreeExtra->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");

    // Create TTree of Sigma Proton Pairs
    fSigmaPairTree = new TTree("fSigmaPairTree","Tree of Sigma Proton Pairs");
    fSigmaPairTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaPairTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaPairTree->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaPairTree->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaPairTree->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaPairTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaPairTree->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaPairTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTree->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaPairTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaPairTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPairTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaPairTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaPairTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaPairTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaPairTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaPairTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaPairTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaPairTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaPairTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaPairTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPairTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaPairTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaPairTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaPairTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaPairTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaPairTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaPairTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaPairTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaPairTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaPairTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaPairTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaPairTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaPairTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaPairTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaPairTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
    fSigmaPairTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaPairTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaPairTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaPairTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaPairTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaPairTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");
    fSigmaPairTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaPairTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaPairTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaPairTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaPairTree->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaPairTree->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaPairTree->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaPairTree->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaPairTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaPairTree->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaPairTree->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaPairTree->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaPairTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaPairTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaPairTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaPairTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaPairTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaPairTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaPairTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaPairTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaPairTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaPairTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaPairTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaPairTree->Branch("fnPair",&fnPair,"fnPair/S");
    fSigmaPairTree->Branch("fnPairlowkstar",&fnPairlowkstar,"fnPairlowkstar/S");
    fSigmaPairTree->Branch("fnPairverylowkstar",&fnPairverylowkstar,"fnPairverylowkstar/S");
    fSigmaPairTree->Branch("fnPairveryverylowkstar",&fnPairveryverylowkstar,"fnPairveryverylowkstar/S");
    fSigmaPairTree->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
    fSigmaPairTree->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    fSigmaPairTree->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTree->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTree->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fSigmaPairTree->Branch("fPairProtonPxatDCA",&fPairProtonPxatDCA,"fPairProtonPxatDCA/F");
    fSigmaPairTree->Branch("fPairProtonPyatDCA",&fPairProtonPyatDCA,"fPairProtonPyatDCA/F");
    fSigmaPairTree->Branch("fPairProtonPzatDCA",&fPairProtonPzatDCA,"fPairProtonPzatDCA/F");
    fSigmaPairTree->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/F");
    fSigmaPairTree->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTree->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTree->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTree->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTree->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTree->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTree->Branch("fPairProtNSigTPCElec",&fPairProtNSigTPCElec,"fPairProtNSigTPCElec/F");
    fSigmaPairTree->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTree->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTree->Branch("fPairProtNSigTOFElec",&fPairProtNSigTOFElec,"fPairProtNSigTOFElec/F");
    fSigmaPairTree->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    fSigmaPairTree->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/I");
    fSigmaPairTree->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/I");
    fSigmaPairTree->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTree->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");

    // Create TTree of Protons
    fProtonTree = new TTree("fProtonTree","Tree of Protons in Sigma Events");
    fProtonTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fProtonTree->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fProtonTree->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fProtonTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fProtonTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fProtonTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fProtonTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fProtonTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fProtonTree->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
    fProtonTree->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    fProtonTree->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fProtonTree->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fProtonTree->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
    fProtonTree->Branch("fPairProtonPxatDCA",&fPairProtonPxatDCA,"fPairProtonPxatDCA/F");
    fProtonTree->Branch("fPairProtonPyatDCA",&fPairProtonPyatDCA,"fPairProtonPyatDCA/F");
    fProtonTree->Branch("fPairProtonPzatDCA",&fPairProtonPzatDCA,"fPairProtonPzatDCA/F");
    fProtonTree->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/F");    
    fProtonTree->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fProtonTree->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fProtonTree->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fProtonTree->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fProtonTree->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fProtonTree->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fProtonTree->Branch("fPairProtNSigTPCElec",&fPairProtNSigTPCElec,"fPairProtNSigTPCElec/F");
    fProtonTree->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fProtonTree->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fProtonTree->Branch("fPairProtNSigTOFElec",&fPairProtNSigTOFElec,"fPairProtNSigTOFElec/F");
    fProtonTree->Branch("fPairProtonChi2",&fPairProtonChi2,"fPairProtonChi2/F");
    fProtonTree->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/I");
    fProtonTree->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/I");

    // Create reduced TTree of Sigma Proton Pairs in Same Event
    fSigmaPairTreeSE = new TTree("fSigmaPairTreeSE","Tree of Sigma Proton Pairs in Same Event");
    fSigmaPairTreeSE->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaPairTreeSE->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaPairTreeSE->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaPairTreeSE->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaPairTreeSE->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaPairTreeSE->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreeSE->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaPairTreeSE->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaPairTreeSE->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreeSE->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreeSE->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaPairTreeSE->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreeSE->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreeSE->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreeSE->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreeSE->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaPairTreeSE->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPairTreeSE->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaPairTreeSE->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreeSE->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaPairTreeSE->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreeSE->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaPairTreeSE->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaPairTreeSE->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaPairTreeSE->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaPairTreeSE->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaPairTreeSE->Branch("fSigpropPx",&fSigpropPx,"fSigpropPx/F");
    fSigmaPairTreeSE->Branch("fSigpropPy",&fSigpropPy,"fSigpropPy/F");
    fSigmaPairTreeSE->Branch("fSigpropPz",&fSigpropPz,"fSigpropPz/F");
    fSigmaPairTreeSE->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaPairTreeSE->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaPairTreeSE->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPairTreeSE->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
    fSigmaPairTreeSE->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
    fSigmaPairTreeSE->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaPairTreeSE->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaPairTreeSE->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaPairTreeSE->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaPairTreeSE->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaPairTreeSE->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaPairTreeSE->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaPairTreeSE->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaPairTreeSE->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaPairTreeSE->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaPairTreeSE->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaPairTreeSE->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaPairTreeSE->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaPairTreeSE->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaPairTreeSE->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaPairTreeSE->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaPairTreeSE->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
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
    fSigmaPairTreeSE->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaPairTreeSE->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaPairTreeSE->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaPairTreeSE->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaPairTreeSE->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaPairTreeSE->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaPairTreeSE->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaPairTreeSE->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreeSE->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreeSE->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaPairTreeSE->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaPairTreeSE->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaPairTreeSE->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaPairTreeSE->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaPairTreeSE->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaPairTreeSE->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaPairTreeSE->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreeSE->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreeSE->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreeSE->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaPairTreeSE->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaPairTreeSE->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaPairTreeSE->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaPairTreeSE->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    fSigmaPairTreeSE->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreeSE->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
    fSigmaPairTreeSE->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    fSigmaPairTreeSE->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreeSE->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreeSE->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/F");
    fSigmaPairTreeSE->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreeSE->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreeSE->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
//    fSigmaPairTreeSE->Branch("fPairProtonPxatDCA",&fPairProtonPxatDCA,"fPairProtonPxatDCA/F");
//    fSigmaPairTreeSE->Branch("fPairProtonPyatDCA",&fPairProtonPyatDCA,"fPairProtonPyatDCA/F");
//    fSigmaPairTreeSE->Branch("fPairProtonPzatDCA",&fPairProtonPzatDCA,"fPairProtonPzatDCA/F");
    fSigmaPairTreeSE->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreeSE->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreeSE->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreeSE->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreeSE->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreeSE->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/I");
    fSigmaPairTreeSE->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/I");

    // Create reduced TTree of Sigma Proton Pairs in Mixed Event
    fSigmaPairTreeME = new TTree("fSigmaPairTreeME","Tree of Sigma Proton Pairs in Mixed Event");
    fSigmaPairTreeME->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaPairTreeME->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaPairTreeME->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaPairTreeME->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaPairTreeME->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaPairTreeME->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaPairTreeME->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaPairTreeME->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaPairTreeME->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaPairTreeME->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaPairTreeME->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaPairTreeME->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaPairTreeME->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaPairTreeME->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaPairTreeME->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaPairTreeME->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaPairTreeME->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaPairTreeME->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaPairTreeME->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaPairTreeME->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaPairTreeME->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaPairTreeME->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaPairTreeME->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaPairTreeME->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaPairTreeME->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaPairTreeME->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaPairTreeME->Branch("fSigpropPx",&fSigpropPx,"fSigpropPx/F");
    fSigmaPairTreeME->Branch("fSigpropPy",&fSigpropPy,"fSigpropPy/F");
    fSigmaPairTreeME->Branch("fSigpropPz",&fSigpropPz,"fSigpropPz/F");
    fSigmaPairTreeME->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaPairTreeME->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaPairTreeME->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaPairTreeME->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
    fSigmaPairTreeME->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
    fSigmaPairTreeME->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
    fSigmaPairTreeME->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaPairTreeME->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaPairTreeME->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaPairTreeME->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaPairTreeME->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaPairTreeME->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaPairTreeME->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaPairTreeME->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaPairTreeME->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaPairTreeME->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaPairTreeME->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaPairTreeME->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaPairTreeME->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaPairTreeME->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaPairTreeME->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaPairTreeME->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaPairTreeME->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaPairTreeME->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaPairTreeME->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaPairTreeME->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaPairTreeME->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
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
    fSigmaPairTreeME->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaPairTreeME->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaPairTreeME->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaPairTreeME->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaPairTreeME->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaPairTreeME->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaPairTreeME->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaPairTreeME->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaPairTreeME->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaPairTreeME->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaPairTreeME->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaPairTreeME->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaPairTreeME->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaPairTreeME->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaPairTreeME->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaPairTreeME->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaPairTreeME->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaPairTreeME->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaPairTreeME->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaPairTreeME->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaPairTreeME->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaPairTreeME->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaPairTreeME->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaPairTreeME->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");
    fSigmaPairTreeME->Branch("fSigmaProtonkstar",&fSigmaProtonkstar,"fSigmaProtonkstar/F");
    fSigmaPairTreeME->Branch("fSigmaProtonpropkstar",&fSigmaProtonpropkstar,"fSigmaProtonpropkstar/F");
    fSigmaPairTreeME->Branch("fPairProtonIsMC",&fPairProtonIsMC,"fPairProtonIsMC/O");
    fSigmaPairTreeME->Branch("fPairProtonIsPrimary",&fPairProtonIsPrimary,"fPairProtonIsPrimary/O");
    fSigmaPairTreeME->Branch("fPairProtonID",&fPairProtonID,"fPairProtonID/I");
    fSigmaPairTreeME->Branch("fPairProtonStatus",&fPairProtonStatus,"fPairProtonStatus/l");
    fSigmaPairTreeME->Branch("fPairProtonCharge",&fPairProtonCharge,"fPairProtonCharge/F");
    fSigmaPairTreeME->Branch("fPairProtonPx",&fPairProtonPx,"fPairProtonPx/F");
    fSigmaPairTreeME->Branch("fPairProtonPy",&fPairProtonPy,"fPairProtonPy/F");
    fSigmaPairTreeME->Branch("fPairProtonPz",&fPairProtonPz,"fPairProtonPz/F");
//    fSigmaPairTreeME->Branch("fPairProtonPxatDCA",&fPairProtonPxatDCA,"fPairProtonPxatDCA/F");
//    fSigmaPairTreeME->Branch("fPairProtonPyatDCA",&fPairProtonPyatDCA,"fPairProtonPyatDCA/F");
//    fSigmaPairTreeME->Branch("fPairProtonPzatDCA",&fPairProtonPzatDCA,"fPairProtonPzatDCA/F");
    fSigmaPairTreeME->Branch("fPairProtonDCAtoPVxy",&fPairProtonDCAtoPVxy,"fPairProtonDCAtoPVxy/F");
    fSigmaPairTreeME->Branch("fPairProtonDCAtoPVz",&fPairProtonDCAtoPVz,"fPairProtonDCAtoPVz/F");
    fSigmaPairTreeME->Branch("fPairProtonNSigTPC",&fPairProtonNSigTPC,"fPairProtonNSigTPC/F");
    fSigmaPairTreeME->Branch("fPairProtonNSigTOF",&fPairProtonNSigTOF,"fPairProtonNSigTOF/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTPCKaon",&fPairProtNSigTPCKaon,"fPairProtNSigTPCKaon/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTOFKaon",&fPairProtNSigTOFKaon,"fPairProtNSigTOFKaon/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTPCPion",&fPairProtNSigTPCPion,"fPairProtNSigTPCPion/F");
    fSigmaPairTreeME->Branch("fPairProtNSigTOFPion",&fPairProtNSigTOFPion,"fPairProtNSigTOFPion/F");
    fSigmaPairTreeME->Branch("fPairProtonCluster",&fPairProtonCluster,"fPairProtonCluster/I");
    fSigmaPairTreeME->Branch("fPairProtonITSCluster",&fPairProtonITSCluster,"fPairProtonITSCluster/I");

    // Create TTree of Sigma Candidate Mixed Event Background
    fSigmaMEBackgroundTree = new TTree("fSigmaMEBackgroundTree","Tree of Sigma Mixed Event Background");
    fSigmaMEBackgroundTree->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaMEBackgroundTree->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaMEBackgroundTree->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaMEBackgroundTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaMEBackgroundTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaMEBackgroundTree->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaMEBackgroundTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaMEBackgroundTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaMEBackgroundTree->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaMEBackgroundTree->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaMEBackgroundTree->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaMEBackgroundTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaMEBackgroundTree->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaMEBackgroundTree->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaMEBackgroundTree->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaMEBackgroundTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaMEBackgroundTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaMEBackgroundTree->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaMEBackgroundTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaMEBackgroundTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaMEBackgroundTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaMEBackgroundTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaMEBackgroundTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaMEBackgroundTree->Branch("fSigpropPx",&fSigpropPx,"fSigpropPx/F");
    fSigmaMEBackgroundTree->Branch("fSigpropPy",&fSigpropPy,"fSigpropPy/F");
    fSigmaMEBackgroundTree->Branch("fSigpropPz",&fSigpropPz,"fSigpropPz/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaMEBackgroundTree->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaMEBackgroundTree->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaMEBackgroundTree->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaMEBackgroundTree->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
    fSigmaMEBackgroundTree->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaMEBackgroundTree->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaMEBackgroundTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaMEBackgroundTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaMEBackgroundTree->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaMEBackgroundTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaMEBackgroundTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaMEBackgroundTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaMEBackgroundTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaMEBackgroundTree->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaMEBackgroundTree->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaMEBackgroundTree->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaMEBackgroundTree->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaMEBackgroundTree->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaMEBackgroundTree->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaMEBackgroundTree->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaMEBackgroundTree->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaMEBackgroundTree->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaMEBackgroundTree->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaMEBackgroundTree->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaMEBackgroundTree->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaMEBackgroundTree->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaMEBackgroundTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaMEBackgroundTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaMEBackgroundTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaMEBackgroundTree->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaMEBackgroundTree->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaMEBackgroundTree->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaMEBackgroundTree->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");

    // Create TTree of rotated Sigma Candidates
    fSigmaCandTreerot = new TTree("fSigmaCandTreerot","Tree of rotated Sigma Candidates");
    fSigmaCandTreerot->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaCandTreerot->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaCandTreerot->Branch("fIsGoodCandidate",&fIsGoodCandidate,"fIsGoodCandidate/O");    
    fSigmaCandTreerot->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaCandTreerot->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaCandTreerot->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaCandTreerot->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaCandTreerot->Branch("fHas4DiffIDs",&fHas4DiffIDs,"fHas4DiffIDs/O");
    fSigmaCandTreerot->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTreerot->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/i");
    fSigmaCandTreerot->Branch("fSigMCLabel",&fSigMCLabel,"fSigMCLabel/I");
    fSigmaCandTreerot->Branch("fSigProtonID",&fSigProtonID,"fSigProtonID/I");
    fSigmaCandTreerot->Branch("fSigProtonStatus",&fSigProtonStatus,"fSigProtonStatus/l");
    fSigmaCandTreerot->Branch("fSigEventID",&fSigEventID,"fSigEventID/l");
    fSigmaCandTreerot->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCandTreerot->Branch("fSigRefMultComb05",&fSigRefMultComb05,"fSigRefMultComb05/S");
    fSigmaCandTreerot->Branch("fSigRefMultComb08",&fSigRefMultComb08,"fSigRefMultComb08/S");
    fSigmaCandTreerot->Branch("fSigRefMultComb10",&fSigRefMultComb10,"fSigRefMultComb10/S");
    fSigmaCandTreerot->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCandTreerot->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaCandTreerot->Branch("fInvSigpropMass",&fInvSigpropMass,"fInvSigpropMass/F");
    fSigmaCandTreerot->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaCandTreerot->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaCandTreerot->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaCandTreerot->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaCandTreerot->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaCandTreerot->Branch("fSigpropPx",&fSigpropPx,"fSigpropPx/F");
    fSigmaCandTreerot->Branch("fSigpropPy",&fSigpropPy,"fSigpropPy/F");
    fSigmaCandTreerot->Branch("fSigpropPz",&fSigpropPz,"fSigpropPz/F");
    fSigmaCandTreerot->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaCandTreerot->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaCandTreerot->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCandTreerot->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaCandTreerot->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaCandTreerot->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaCandTreerot->Branch("fSigDecayVertXMC",&fSigDecayVertXMC,"fSigDecayVertXMC/F");
    fSigmaCandTreerot->Branch("fSigDecayVertYMC",&fSigDecayVertYMC,"fSigDecayVertYMC/F");
    fSigmaCandTreerot->Branch("fSigDecayVertZMC",&fSigDecayVertZMC,"fSigDecayVertZMC/F");
    fSigmaCandTreerot->Branch("fSigPxMC",&fSigPxMC,"fSigPxMC/F");
    fSigmaCandTreerot->Branch("fSigPyMC",&fSigPyMC,"fSigPyMC/F");
    fSigmaCandTreerot->Branch("fSigPzMC",&fSigPzMC,"fSigPzMC/F");
    fSigmaCandTreerot->Branch("fPrimVertXMC",&fPrimVertXMC,"fPrimVertXMC/F");
    fSigmaCandTreerot->Branch("fPrimVertYMC",&fPrimVertYMC,"fPrimVertYMC/F");
    fSigmaCandTreerot->Branch("fPrimVertZMC",&fPrimVertZMC,"fPrimVertZMC/F");
    fSigmaCandTreerot->Branch("fPhoton1Px",&fPhoton1Px,"fPhoton1Px/F");
    fSigmaCandTreerot->Branch("fPhoton1Py",&fPhoton1Py,"fPhoton1Py/F");
    fSigmaCandTreerot->Branch("fPhoton1Pz",&fPhoton1Pz,"fPhoton1Pz/F");
    fSigmaCandTreerot->Branch("fPhoton2Px",&fPhoton2Px,"fPhoton2Px/F");
    fSigmaCandTreerot->Branch("fPhoton2Py",&fPhoton2Py,"fPhoton2Py/F");
    fSigmaCandTreerot->Branch("fPhoton2Pz",&fPhoton2Pz,"fPhoton2Pz/F");
    fSigmaCandTreerot->Branch("fPhotonDaughtMaxEta",&fPhotonDaughtMaxEta,"fPhotonDaughtMaxEta/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxDeltaTheta",&fPhotonsMaxDeltaTheta,"fPhotonsMaxDeltaTheta/F");
    fSigmaCandTreerot->Branch("fPhoton1CPA",&fPhoton1CPA,"fPhoton1CPA/F");
    fSigmaCandTreerot->Branch("fPhoton2CPA",&fPhoton2CPA,"fPhoton2CPA/F");
    fSigmaCandTreerot->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaCandTreerot->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaCandTreerot->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaCandTreerot->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
    fSigmaCandTreerot->Branch("fPhotonsMinCluster",&fPhotonsMinCluster,"fPhotonsMinCluster/F");
    fSigmaCandTreerot->Branch("fPhotonsMinITSCluster",&fPhotonsMinITSCluster,"fPhotonsMinITSCluster/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxalpha",&fPhotonsMaxalpha,"fPhotonsMaxalpha/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxqt",&fPhotonsMaxqt,"fPhotonsMaxqt/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxOpenAngle",&fPhotonsMaxOpenAngle,"fPhotonsMaxOpenAngle/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxinvmass",&fPhotonsMaxinvmass,"fPhotonsMaxinvmass/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxNSigTPC",&fPhotonsMaxNSigTPC,"fPhotonsMaxNSigTPC/F");
    fSigmaCandTreerot->Branch("fPhotonsMaxChi2",&fPhotonsMaxChi2,"fPhotonsMaxChi2/F");    
    fSigmaCandTreerot->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaCandTreerot->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaCandTreerot->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaCandTreerot->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaCandTreerot->Branch("fPi0DecayVertX",&fPi0DecayVertX,"fPi0DecayVertX/F");
    fSigmaCandTreerot->Branch("fPi0DecayVertY",&fPi0DecayVertY,"fPi0DecayVertY/F");
    fSigmaCandTreerot->Branch("fPi0DecayVertZ",&fPi0DecayVertZ,"fPi0DecayVertZ/F");
    fSigmaCandTreerot->Branch("fPi0PhotPhotDCA",&fPi0PhotPhotDCA,"fPi0PhotPhotDCA/F");
    fSigmaCandTreerot->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCandTreerot->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCandTreerot->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
//    fSigmaCandTreerot->Branch("fProtonPxatDCA",&fProtonPxatDCA,"fProtonPxatDCA/F");
//    fSigmaCandTreerot->Branch("fProtonPyatDCA",&fProtonPyatDCA,"fProtonPyatDCA/F");
//    fSigmaCandTreerot->Branch("fProtonPzatDCA",&fProtonPzatDCA,"fProtonPzatDCA/F");
    fSigmaCandTreerot->Branch("fProtonpropPx",&fProtonpropPx,"fProtonpropPx/F");
    fSigmaCandTreerot->Branch("fProtonpropPy",&fProtonpropPy,"fProtonpropPy/F");
    fSigmaCandTreerot->Branch("fProtonpropPz",&fProtonpropPz,"fProtonpropPz/F");
    fSigmaCandTreerot->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTreerot->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTreerot->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaCandTreerot->Branch("fProtonNSigTPC",&fProtonNSigTPC,"fProtonNSigTPC/F");
    fSigmaCandTreerot->Branch("fProtonNSigTOF",&fProtonNSigTOF,"fProtonNSigTOF/F");
    fSigmaCandTreerot->Branch("fProtonNCluster",&fProtonNCluster,"fProtonNCluster/I");
    fSigmaCandTreerot->Branch("fProtonNITSCluster",&fProtonNITSCluster,"fProtonNITSCluster/I");
    fSigmaCandTreerot->Branch("fProtonChi2",&fProtonChi2,"fProtonChi2/F");
    fSigmaCandTreerot->Branch("fProtonNSigTPCPion",&fProtonNSigTPCPion,"fProtonNSigTPCPion/F");
    fSigmaCandTreerot->Branch("fProtonNSigTPCKaon",&fProtonNSigTPCKaon,"fProtonNSigTPCKaon/F");
    fSigmaCandTreerot->Branch("fProtonNSigTPCElec",&fProtonNSigTPCElec,"fProtonNSigTPCElec/F");
    fSigmaCandTreerot->Branch("fProtonNSigTOFPion",&fProtonNSigTOFPion,"fProtonNSigTOFPion/F");
    fSigmaCandTreerot->Branch("fProtonNSigTOFKaon",&fProtonNSigTOFKaon,"fProtonNSigTOFKaon/F");
    fSigmaCandTreerot->Branch("fProtonNSigTOFElec",&fProtonNSigTOFElec,"fProtonNSigTOFElec/F");

    //Save Particle Masses and other constants for later use
    cElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();      
    cProtonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    
    cSigmaMass    = TDatabasePDG::Instance()->GetParticle(3222)->Mass();    
    cPi0Mass      = TDatabasePDG::Instance()->GetParticle(111)->Mass();     
    c = 2.99792457999999984e-02; // [cm/ps]

/**************************Histograms********************************/

    //Book Keeper for used Cuts 
    TH1D* fHistCutBookKeeper           = new TH1D("fHistCutBookKeeper", "Book Keeper for used Cuts", 84, 0.5, 84.5);

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
    TH1F* fHistProtonStatistics        = new TH1F("fHistProtonStatistics", "Proton Counter;Stage;Counts",7,0.5,7.5);
    TH1F* fHistProtonStatisticsMC      = new TH1F("fHistProtonStatisticsMC", "Proton Counter MC;Stage;Counts",7,0.5,7.5);
    TH1F* fHistProtonStatisticsSigmaMC = new TH1F("fHistProtonStatisticsSigmaMC", "Proton Counter Sigma MC;Stage;Counts",7,0.5,7.5);
    TH1F* fHistAddV0Statistics         = new TH1F("fHistAddV0Statistics", "Additional V0 Counter;Stage;Counts",11,0.5,11.5);
    TH1F* fHistAddV0StatisticsMC       = new TH1F("fHistAddV0StatisticsMC", "Additional V0 Counter MC;Stage;Counts",11,0.5,11.5);
    TH1F* fHistAddV0StatisticsSigmaMC  = new TH1F("fHistAddV0StatisticsSigmaMC", "Additional V0 Counter Sigma MC;Stage;Counts",11,0.5,11.5);
    TH1F* fHistGammaPairStats          = new TH1F("fHistGammaPairStats", "Gamma Pair Counter;;Counts",6,0.5,6.5);
    TH1F* fHistGammaPairStatsOneadd    = new TH1F("fHistGammaPairStatsOneadd", "Gamma Pair Counter, One Additional;;Counts",6,0.5,6.5);
    TH1F* fHistGammaPairStatsOnlyadd   = new TH1F("fHistGammaPairStatsOnlyadd", "Gamma Pair Counter, One Additional;;Counts",6,0.5,6.5);
    TH1F* fHistSigmaCounter            = new TH1F("fHistSigmaCounter", "#Sigma Counter", 6, 0.5, 6.5);
                                       
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
    TH1F* fHistMeanProtonpt            = new TH1F("fHistMeanProtonpt","Mean transverse momentum of selected Protons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",5002,-0.002,5);   

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

    //Additional Photon QA                    
    TH1F* fHistPairDCAtoPV             = new TH1F("fHistPairDCAtoPV","DCA to PV of Pair Daughters;DCA (cm);Counts/(0.1 mm)",2000,0,20); 
    TH1F* fHistPairDCA                 = new TH1F("fHistPairDCA","DCA between Pair Daughters;DCA (cm);Counts/(0.5mm)",200,0,10); 
    TH1F* fHistPairDCAMC               = new TH1F("fHistPairDCAMC","DCA between Pair Daughters MC;DCA (cm);Counts/(0.5mm)",200,0,10);
    TH1F* fHistPairRadius              = new TH1F("fHistPairRadius","Radius of Pairs;r_{xy} (cm);Counts/mm",3000,0,300);      
    TH1F* fHistPairRadiusMC            = new TH1F("fHistPairRadiusMC","Radius of Pairs MC;r_{xy} (cm);Counts/mm",3000,0,300);   
    TH1F* fHistPairPhotonDCAPV         = new TH1F("fHistPairPhotonDCAPV","DCA to PV of Photons;DCA (cm);Counts/mm)",100,0,10);
    TH1F* fHistPairPhotonDCAPVMC       = new TH1F("fHistPairPhotonDCAPVMC","DCA to PV of MC Photons;DCA (cm);Counts/mm)",100,0,10);
    TH1F* fHistPairPhotonDCAPVMCSigma  = new TH1F("fHistPairPhotonDCAPVMCSigma","DCA to PV of MC Photons from #Sigma^{+}/#bar#Sigma^{-};DCA (cm);Counts/mm)",100,0,10);
    TH1F* fHistPairCPA                 = new TH1F("fHistPairCPA","Cosine of Pointing Angle of Pairs;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistPairCPAMC               = new TH1F("fHistPairCPAMC","Cosine of Pointing Angle of Pairs MC;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistPairCPAMCSigma          = new TH1F("fHistPairCPAMCSigma","Cosine of Pointing Angle of Pairs MC (#gamma from #Sigma);cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH2F* fHistPairArmPod              = new TH2F("fHistPairArmPod", "Armenteros-Podolanski Plot of Pairs;#alpha;q_{t} (GeV/#it{c})",300,-1.5,1.5,200,-0.05,0.4);       
    TH2F* fHistPairArmPodMC            = new TH2F("fHistPairArmPodMC", "Armenteros-Podolanski Plot of Pairs MC;#alpha;q_{t} (GeV/#it{c})",300,-1.5,1.5,200,-0.05,0.4);
    TH1F* fHistPairOpenAngle           = new TH1F("fHistPairOpenAngle","Total Opening Angle of Pairs;#xi (rad);Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistPairOpenAngleMC         = new TH1F("fHistPairOpenAngleMC","Total Opening Angle of Pairs MC;#xi (rad);Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistPairDeltaTheta          = new TH1F("fHistPairDeltaTheta","#Delta#Theta of Pair Daughters;#Delta#Theta (rad);Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistPairDeltaThetaMC        = new TH1F("fHistPairDeltaThetaMC","#Delta#Theta of Pair Daughters MC;#Delta#Theta (rad);Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistPairInvMass             = new TH1F("fHistPairInvMass", "Invariant mass of Pairs;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistPairInvMassMC           = new TH1F("fHistPairInvMassMC", "Invariant mass of Pairs MC;m_{inv} (GeV/#it{c}^{2});Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistPairptMC                = new TH1F("fHistPairptMC","Transverse momentum of MC Pair Photons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistPairptwCutsMC           = new TH1F("fHistPairptwCutsMC","Transverse momentum of MC Pair Photons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistPairptwCuts             = new TH1F("fHistPairptwCuts","Transverse momentum of Pair Photons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);   

    //Electron QA                      
    TH2F* fHistElectronEtaPhiMC        = new TH2F("fHistElectronEtaPhiMC","#eta vs. #phi of Electrons MC;#eta;#phi (rad)",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistElectronChi2MC          = new TH1F("fHistElectronChi2MC","#chi^{2}/NDF of Electrons MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistElectronTPCClusterMC    = new TH1F("fHistElectronTPCClusterMC","Number of TPC Clusters of Electrons MC",161,-0.5,160.5);
    TH1F* fHistElectronITSClusterMC    = new TH1F("fHistElectronITSClusterMC","Number of ITS Clusters of Electrons MC",13,-0.5,12.5);
    TH2F* fHistElectronpvsNSigmaTPC    = new TH2F("fHistElectronpvsNSigmaTPC","Momentum vs N sigma TPC Electrons;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigTPCwCuts = new TH2F("fHistElectronpvsNSigTPCwCuts","Momentum vs N sigma TPC Electrons with Cuts;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTPCMC  = new TH2F("fHistElectronpvsNSigmaTPCMC","Momentum vs N sigma TPC Electrons MC;p (GeV/#it{c});n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTOF    = new TH2F("fHistElectronpvsNSigmaTOF","Momentum vs N sigma TOF Electrons;p (GeV/#it{c});n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTOFMC  = new TH2F("fHistElectronpvsNSigmaTOFMC","Momentum vs N sigma TOF Electrons MC;p (GeV/#it{c});n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH1F* fHistElectronptMC            = new TH1F("fHistElectronptMC","Transverse momentum of MC Electrons;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistElectronptwCutsMC       = new TH1F("fHistElectronptwCutsMC","Transverse momentum of MC Electrons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistElectronptwCuts         = new TH1F("fHistElectronptwCuts","Transverse momentum of Electrons with Cuts;#it{p}_{T} (GeV/#it{c});Counts/(MeV/#it{c})",12000,0,12);   

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
    TH1F* fHistInvSigmaMassmix         = new TH1F("fHistInvSigmaMassmix","Invariant mass of mixed #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistInvSigmaMassrot         = new TH1F("fHistInvSigmaMassrot","Invariant mass of rotated #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} (GeV/#it{c}^{2});Counts/(10 MeV/#it{c}^{2})", 500, 0.5, 5.5);
    TH1F* fHistMCOneGammaSigmaPA       = new TH1F("fHistMCOneGammaSigmaPA","#Sigma^{+} PA, One #gamma MC;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistMCPrimOneGammaSigmaPA   = new TH1F("fHistMCPrimOneGammaSigmaPA","#Sigma^{+} PA, One #gamma primary MC;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistOneGammaSigmaPA         = new TH1F("fHistOneGammaSigmaPA","#Sigma^{+} PA, One #gamma;PA (rad);Counts/(0.005)",600,0,3);
    TH1F* fHistPi0VertexMC             = new TH1F("fHistPi0VertexMC","#pi^{0} MC Decay Radius;r (cm);Counts/(mm)", 200, 0, 20);
    TH1F* fHistMCSigmaProtonkstar      = new TH1F("fHistMCSigmaProtonkstar","MC #Sigma^{+}/#bar#Sigma^{-}-p k*;k* (GeV/#it{c});Counts/(5 MeV/#it{c})",1000,0,5);
    TH1F* fHistSigmaProtonkstar        = new TH1F("fHistSigmaProtonkstar","#Sigma^{+}-p k*;k* (GeV/#it{c});Counts/(5 MeV/#it{c})",1000,0,5);
    TH1F* fHistAntiSigmaProtonkstar    = new TH1F("fHistAntiSigmaProtonkstar","#bar{#Sigma^{-}}-#bar{p} k*;k* (GeV/#it{c});Counts/(5 MeV/#it{c})",1000,0,5);

    //Sigma Momentum Resolution
    TH1F* fHistSigmaPxResnoprop        = new TH1F("fHistSigmaPxResnoprop","#Sigma p_{x} - MC #Sigma p_{x}, no propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResnoprop        = new TH1F("fHistSigmaPyResnoprop","#Sigma p_{y} - MC #Sigma p_{y}, no propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResnoprop        = new TH1F("fHistSigmaPzResnoprop","#Sigma p_{z} - MC #Sigma p_{z}, no propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxResprop          = new TH1F("fHistSigmaPxResprop","#Sigma p_{x} - MC #Sigma p_{x}, ETP propagation;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResprop          = new TH1F("fHistSigmaPyResprop","#Sigma p_{y} - MC #Sigma p_{y}, ETP propagation;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResprop          = new TH1F("fHistSigmaPzResprop","#Sigma p_{z} - MC #Sigma p_{z}, ETP propagation;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxResdEdxcorr      = new TH1F("fHistSigmaPxResdEdxcorr","#Sigma p_{x} - MC #Sigma p_{x}, propagation with dEdx correction;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyResdEdxcorr      = new TH1F("fHistSigmaPyResdEdxcorr","#Sigma p_{y} - MC #Sigma p_{y}, propagation with dEdx correction;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzResdEdxcorr      = new TH1F("fHistSigmaPzResdEdxcorr","#Sigma p_{z} - MC #Sigma p_{z}, propagation with dEdx correction;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPxRespropKF        = new TH1F("fHistSigmaPxRespropKF","#Sigma p_{x} - MC #Sigma p_{x}, propagation with KFParticle;#Deltap_{x} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPyRespropKF        = new TH1F("fHistSigmaPyRespropKF","#Sigma p_{y} - MC #Sigma p_{y}, propagation with KFParticle;#Deltap_{y} (MeV/c);Counts/(MeV/c)", 400, -100, 100);
    TH1F* fHistSigmaPzRespropKF        = new TH1F("fHistSigmaPzRespropKF","#Sigma p_{z} - MC #Sigma p_{z}, propagation with KFParticle;#Deltap_{z} (MeV/c);Counts/(MeV/c)", 400, -100, 100);

    //Event Mixing
    TH3F* fHistPairNMixedEvents1       = new TH3F("fHistPairNMixedEvents1","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins2, fMinCentBin2, fMaxCentBin2, nZvtxBins2, fMinZBin2, fMaxZBin2, fEvPoolSize2+1, -0.5, fEvPoolSize2+0.5);
    TH3F* fHistPairNMixedEvents2       = new TH3F("fHistPairNMixedEvents2","Number of Mixed Events;Mean Proton #it{p}_{T} (GeV/#it{c});fRefMultComb08;N Mixed Events", nPtBins2, fMinPtBin2, fMaxPtBin2, nCentralityBins2, fMinCentBin2, fMaxCentBin2, fEvPoolSize2+1, -0.5, fEvPoolSize2+0.5);
    TH3F* fHistPairNMixedEvents3       = new TH3F("fHistPairNMixedEvents3","Number of Mixed Events;Mean Proton #it{p}_{T} (GeV/#it{c});Vertex Z (cm);N Mixed Events", nPtBins2, fMinPtBin2, fMaxPtBin2, nZvtxBins2, fMinZBin2, fMaxZBin2, fEvPoolSize2+1, -0.5, fEvPoolSize2+0.5);
    TH3F* fHistBkgNMixedEvents         = new TH3F("fHistBkgNMixedEvents","Number of Mixed Events;fRefMultComb08;Vertex Z (cm);N Mixed Events", nCentralityBins, fMinCentBin, fMaxCentBin, nZvtxBins, fMinZBin, fMaxZBin, fEvPoolSize+1, -0.5, fEvPoolSize+0.5);

    //KF Checks
    TH1F* fHistPhotonKFCheck           = new TH1F("fHistPhotonKFCheck","Check Photon Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistPi0KFCheck              = new TH1F("fHistPi0KFCheck","Check Pion Daughters to avoid floating point exceptions", 4, 0.5, 4.5);
    TH1F* fHistSigmaKFCheck            = new TH1F("fHistSigmaKFCheck","Check Sigma Daughters to avoid floating point exceptions", 4, 0.5, 4.5);

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
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(19,"fPtBins2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(20,"fMinPtBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(21,"fMaxPtBin2");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(22,"fRequireSigma");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(23,"fRequireProton");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(24,"fMaxProtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(25,"fMinTPCClustProt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(26,"fMaxNsigProtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(27,"fRequireProtonTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(28,"fRequireProtonTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(29,"fRequireProtonTOFforPairs");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(30,"fMaxNsigProtTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(31,"fMaxpOnlyTPCPID");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(32,"fMinProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(33,"fMaxProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(34,"fStrictMaxProtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(35,"fStrictMinTPCClustProt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(36,"fStrictMaxNsigProtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(37,"fStrictMaxNsigProtTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(38,"fStrictMaxpOnlyTPCPID");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(39,"fStrictMinProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(40,"fStrictMaxProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(41,"fMaxMCEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(42,"fMaxDaughtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(43,"fMinTPCClustDaught");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(44,"fMaxNsigDaughtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(45,"fMaxalpha");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(46,"fMaxqt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(47,"fMaxopenangle");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(48,"fMaxdeltatheta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(49,"fMinV0CPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(50,"fMinV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(51,"fMaxV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(52,"fMaxphotonmass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(53,"fMinDCADaughtPV");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(54,"fMaxDCADaught");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(55,"fMaxElecEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(56,"fMinTPCClustElec");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(57,"fMaxNsigElecTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(58,"fMinNsigHadronTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(59,"fMaxNsigElecTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(60,"fMaxElecpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(61,"fCleanAutoCorr");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(62,"fMinPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(63,"fMaxPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(64,"fMaxSigmaPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(65,"fMaxSigmaY");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(66,"fMaxSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(67,"fMinProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(68,"fMinProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(69,"fMaxProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(70,"fMaxProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(71,"fRequireDCACut");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(72,"flowkstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(73,"fverylowkstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(74,"fveryverylowkstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(75,"fMinCorrPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(76,"fMaxCorrPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(77,"fMaxCorrSigmaPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(78,"fMinCorrSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(79,"fMaxCorrSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(80,"fMinCorrProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(81,"fMaxCorrPairProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(82,"fMaxCorrPairProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(83,"fMaxCorrkstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(84,"Number of Fills");

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
    fHistCutBookKeeper->SetBinContent(19,fPtBins2);
    fHistCutBookKeeper->SetBinContent(20,fMinPtBin2);
    fHistCutBookKeeper->SetBinContent(21,fMaxPtBin2);
    fHistCutBookKeeper->SetBinContent(22,fRequireSigma);
    fHistCutBookKeeper->SetBinContent(23,fRequireProton);
    fHistCutBookKeeper->SetBinContent(24,fMaxProtEta);
    fHistCutBookKeeper->SetBinContent(25,fMinTPCClustProt);
    fHistCutBookKeeper->SetBinContent(26,fMaxNsigProtTPC);
    fHistCutBookKeeper->SetBinContent(27,fRequireProtonTPC);
    fHistCutBookKeeper->SetBinContent(28,fRequireProtonTOF);
    fHistCutBookKeeper->SetBinContent(29,fRequireProtonTOFforPairs);
    fHistCutBookKeeper->SetBinContent(30,fMaxNsigProtTOF);
    fHistCutBookKeeper->SetBinContent(31,fMaxpOnlyTPCPID);
    fHistCutBookKeeper->SetBinContent(32,fMinProtpt);
    fHistCutBookKeeper->SetBinContent(33,fMaxProtpt);
    fHistCutBookKeeper->SetBinContent(34,fStrictMaxProtEta);
    fHistCutBookKeeper->SetBinContent(35,fStrictMinTPCClustProt);
    fHistCutBookKeeper->SetBinContent(36,fStrictMaxNsigProtTPC);
    fHistCutBookKeeper->SetBinContent(37,fStrictMaxNsigProtTOF);
    fHistCutBookKeeper->SetBinContent(38,fStrictMaxpOnlyTPCPID);
    fHistCutBookKeeper->SetBinContent(39,fStrictMinProtpt);
    fHistCutBookKeeper->SetBinContent(40,fStrictMaxProtpt);
    fHistCutBookKeeper->SetBinContent(41,fMaxMCEta);
    fHistCutBookKeeper->SetBinContent(42,fMaxDaughtEta);
    fHistCutBookKeeper->SetBinContent(43,fMinTPCClustDaught);
    fHistCutBookKeeper->SetBinContent(44,fMaxNsigDaughtTPC);
    fHistCutBookKeeper->SetBinContent(45,fMaxalpha);
    fHistCutBookKeeper->SetBinContent(46,fMaxqt);
    fHistCutBookKeeper->SetBinContent(47,fMaxopenangle);
    fHistCutBookKeeper->SetBinContent(48,fMaxdeltatheta);
    fHistCutBookKeeper->SetBinContent(49,fMinV0CPA);
    fHistCutBookKeeper->SetBinContent(50,fMinV0Radius);
    fHistCutBookKeeper->SetBinContent(51,fMaxV0Radius);
    fHistCutBookKeeper->SetBinContent(52,fMaxphotonmass);
    fHistCutBookKeeper->SetBinContent(53,fMinDCADaughtPV);
    fHistCutBookKeeper->SetBinContent(54,fMaxDCADaught);
    fHistCutBookKeeper->SetBinContent(55,fMaxElecEta);
    fHistCutBookKeeper->SetBinContent(56,fMinTPCClustElec);
    fHistCutBookKeeper->SetBinContent(57,fMaxNsigElecTPC);
    fHistCutBookKeeper->SetBinContent(58,fMinNsigHadronTPC);
    fHistCutBookKeeper->SetBinContent(59,fMaxNsigElecTOF);
    fHistCutBookKeeper->SetBinContent(60,fMaxElecpt);
    fHistCutBookKeeper->SetBinContent(61,fCleanAutoCorr); 
    fHistCutBookKeeper->SetBinContent(62,fMinPi0Mass);
    fHistCutBookKeeper->SetBinContent(63,fMaxPi0Mass);
    fHistCutBookKeeper->SetBinContent(64,fMaxSigmaPA);
    fHistCutBookKeeper->SetBinContent(65,fMaxSigmaY);
    fHistCutBookKeeper->SetBinContent(66,fMaxSigmaMass);
    fHistCutBookKeeper->SetBinContent(67,fMinProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(68,fMinProtonDCAz);
    fHistCutBookKeeper->SetBinContent(69,fMaxProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(70,fMaxProtonDCAz);
    fHistCutBookKeeper->SetBinContent(71,fRequireDCACut);
    fHistCutBookKeeper->SetBinContent(72,flowkstar);
    fHistCutBookKeeper->SetBinContent(73,fverylowkstar);
    fHistCutBookKeeper->SetBinContent(74,fveryverylowkstar);
    fHistCutBookKeeper->SetBinContent(75,fMinCorrPi0Mass);
    fHistCutBookKeeper->SetBinContent(76,fMaxCorrPi0Mass);
    fHistCutBookKeeper->SetBinContent(77,fMaxCorrSigmaPA);
    fHistCutBookKeeper->SetBinContent(78,fMinCorrSigmaMass);
    fHistCutBookKeeper->SetBinContent(79,fMaxCorrSigmaMass);
    fHistCutBookKeeper->SetBinContent(80,fMinCorrProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(81,fMaxCorrPairProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(82,fMaxCorrPairProtonDCAz);
    fHistCutBookKeeper->SetBinContent(83,fMaxCorrkstar);
    fHistCutBookKeeper->SetBinContent(84,1);

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
    fHistProtonStatistics->GetXaxis()->SetBinLabel(2, "No doublecount");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(3, "Passed eta cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(4, "Passed cluster cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(5, "Passed TPC cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(6, "Passed TOF cut");
    fHistProtonStatistics->GetXaxis()->SetBinLabel(7, "Passed pt cut");

    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(1, "Tracks");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(2, "No doublecount");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(3, "Passed eta cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(4, "Passed cluster cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(5, "Passed TPC cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(6, "Passed TOF cut");
    fHistProtonStatisticsMC->GetXaxis()->SetBinLabel(7, "Passed pt cut");

    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(1, "Tracks");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(2, "No doublecount");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(3, "Passed eta cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(4, "Passed cluster cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(5, "Passed TPC cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(6, "Passed TOF cut");
    fHistProtonStatisticsSigmaMC->GetXaxis()->SetBinLabel(7, "Passed pt cut");

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
    fOutputList->Add(fHistMeanProtonpt);
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
    //Additional Photon QA
    fOutputList->Add(fHistPairDCAtoPV);
    fOutputList->Add(fHistPairDCA);
    fOutputList->Add(fHistPairDCAMC);
    fOutputList->Add(fHistPairRadius);
    fOutputList->Add(fHistPairRadiusMC);
    fOutputList->Add(fHistPairPhotonDCAPV);
    fOutputList->Add(fHistPairPhotonDCAPVMC);
    fOutputList->Add(fHistPairPhotonDCAPVMCSigma);
    fOutputList->Add(fHistPairCPA);
    fOutputList->Add(fHistPairCPAMC);
    fOutputList->Add(fHistPairCPAMCSigma);
    fOutputList->Add(fHistPairArmPod);
    fOutputList->Add(fHistPairArmPodMC);
    fOutputList->Add(fHistPairOpenAngle);
    fOutputList->Add(fHistPairOpenAngleMC);
    fOutputList->Add(fHistPairDeltaTheta);
    fOutputList->Add(fHistPairDeltaThetaMC);
    fOutputList->Add(fHistPairInvMass);
    fOutputList->Add(fHistPairInvMassMC);
    fOutputList->Add(fHistPairptMC);
    fOutputList->Add(fHistPairptwCutsMC);
    fOutputList->Add(fHistPairptwCuts);
    //Electron QA
    fOutputList->Add(fHistElectronEtaPhiMC);
    fOutputList->Add(fHistElectronChi2MC);
    fOutputList->Add(fHistElectronTPCClusterMC);
    fOutputList->Add(fHistElectronITSClusterMC);
    fOutputList->Add(fHistElectronpvsNSigmaTPC);
    fOutputList->Add(fHistElectronpvsNSigTPCwCuts);
    fOutputList->Add(fHistElectronpvsNSigmaTPCMC);
    fOutputList->Add(fHistElectronpvsNSigmaTOF);
    fOutputList->Add(fHistElectronpvsNSigmaTOFMC);
    fOutputList->Add(fHistElectronptMC);
    fOutputList->Add(fHistElectronptwCutsMC);
    fOutputList->Add(fHistElectronptwCuts);
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
    fOutputList->Add(fHistInvSigmaMassmix);         
    fOutputList->Add(fHistInvSigmaMassrot);         
    fOutputList->Add(fHistMCOneGammaSigmaPA);
    fOutputList->Add(fHistMCPrimOneGammaSigmaPA);
    fOutputList->Add(fHistOneGammaSigmaPA);
    fOutputList->Add(fHistPi0VertexMC);
    fOutputList->Add(fHistMCSigmaProtonkstar);
    fOutputList->Add(fHistSigmaProtonkstar);
    fOutputList->Add(fHistAntiSigmaProtonkstar);

    //Sigma Momentum Resolution 
    fOutputList->Add(fHistSigmaPxResnoprop);
    fOutputList->Add(fHistSigmaPyResnoprop);
    fOutputList->Add(fHistSigmaPzResnoprop);
    fOutputList->Add(fHistSigmaPxResprop);
    fOutputList->Add(fHistSigmaPyResprop);
    fOutputList->Add(fHistSigmaPzResprop);
    fOutputList->Add(fHistSigmaPxResdEdxcorr);
    fOutputList->Add(fHistSigmaPyResdEdxcorr);
    fOutputList->Add(fHistSigmaPzResdEdxcorr);
    fOutputList->Add(fHistSigmaPxRespropKF);
    fOutputList->Add(fHistSigmaPyRespropKF);
    fOutputList->Add(fHistSigmaPzRespropKF);

    //Event Mixing
    fOutputList->Add(fHistPairNMixedEvents1);
    fOutputList->Add(fHistPairNMixedEvents2);
    fOutputList->Add(fHistPairNMixedEvents3);
    fOutputList->Add(fHistBkgNMixedEvents);

    //KF Checks
    fOutputList->Add(fHistPhotonKFCheck);
    fOutputList->Add(fHistPi0KFCheck);
    fOutputList->Add(fHistSigmaKFCheck);

/**************************************************************************/
    PostData(1, fOutputList);         
    PostData(2, fSigmaCandTree);         
    PostData(3, fSigmaCandTreeExtra);         
    PostData(4, fSigmaPairTree);         
    PostData(5, fProtonTree);         
    PostData(6, fSigmaPairTreeSE);         
    PostData(7, fSigmaPairTreeME);         
    PostData(8, fSigmaMEBackgroundTree);         
    PostData(9, fSigmaCandTreerot);         

}//end of UserCreateOutputObjects()

//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlus::UserExec(Option_t *)
{  

  // Main loop. Called once for each event
  if(fDebug) cout << "Now in UserExec(). Start of Event\n";

  // Reset bools
  fEventhasSigma = kFALSE;
  fEventhasProton = kFALSE;

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

  Bool_t isSimPileup = kFALSE;
  AliAODMCHeader* aodMCHeader = NULL;
  aodMCHeader=(AliAODMCHeader*)mcEvent->GetHeader();
  if(aodMCHeader){
    TString genname = (TString)aodMCHeader->GetGeneratorName();
    isSimPileup = AliAnalysisUtils::IsPileupInGeneratedEvent(aodMCHeader,genname);
  }
  FillHistogram("fHistMCGenPileup",isSimPileup);
  if(fRemoveGenPileup&&isSimPileup) return;

  if(isMonteCarlo) {
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
    // Check MC primary vertex position
    Double_t primaryVtxPosMC[3] = {-999,-999,-999};
    const AliAODVertex *aodVtxMC = aodEvent->GetPrimaryVertex();
    if(!aodVtxMC) {AliWarning("No primary vertex in MC!");}
    else aodVtxMC->GetXYZ(primaryVtxPosMC);
    primaryVtxPosXMC=primaryVtxPosMC[0];
    primaryVtxPosYMC=primaryVtxPosMC[1];
    primaryVtxPosZMC=primaryVtxPosMC[2];
    FillHistogram("fHistVertexZMC",primaryVtxPosZMC);
  }
  
  //Magnetic Field
  Bz = aodEvent->GetMagneticField();    
  // Set Magnetic field for ALL KFParticles
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

  //Process Protons
  if(fDebug&&fProcessProtons) cout << "Processing Protons\n";
  if(fProcessProtons) FillProtonArray();

  //Process MC Particles
  if(fDebug&&isMonteCarlo&&fProcessMCParticles) cout << "Processing MC Particles\n";
  if(isMonteCarlo && fProcessMCParticles) ProcessMCParticles();

  //Process V0s
  if(fDebug&&fProcessV0s) cout << "Processing V0s\n";
  if(fProcessV0s) FillV0PhotonArray();

  //Find Additional Photons
  if(fDebug&&fProcessAddPhoton) cout << "Searching for additional Photons\n";
  if(fProcessAddPhoton) FindAddPhotons();

  //Process Electrons
  if(fDebug&&fProcessElectrons) cout << "Processing Electrons\n";
  if(fProcessElectrons) FillElectronArray(); 

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessReco) cout << "Reconstructing Pi0s and Sigmas\n";
  if(fProcessReco) ReconstructParticles();

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessRecoOff) cout << "Reconstructing Pi0s and Sigmas. One offline Photon\n";
  if(fProcessRecoOff) ReconstructParticlesOff();

  //Reconstruct Pi0 and Sigma+
  if(fDebug&&fProcessRecoOff2) cout << "Reconstructing Pi0s and Sigmas. Two offline Photons\n";
  if(fProcessRecoOff2) ReconstructParticlesOff2();

  //Reconstruct Sigma+
  if(fDebug&&fProcessOneGamma) cout << "Reconstructing Pi0s and Sigmas. Exotic decay channel\n";
  if(fProcessOneGamma) ReconstructParticlesOneGamma();

  //Update Event Pool at the End of the Event if Event Mixing is enabled
  if(fSaveMixedBackground){
  if(fDebug) cout << "Updating Event Pool\n";
   
    //Get Pool from Pool Manager for given RefMult and Z Vertex values
    AliEventPool* Evpool = 0x0;
	  if(fEvPoolMgr) Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
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
        redprot->InitfromTrack(prot, fPIDResponse);
        ProtonObjArray->Add(redprot);
      }

      //Clone it and update the Pool    
      TObjArray* ProtonCloneArray = (TObjArray*)ProtonObjArray->Clone();    
      ProtonCloneArray->SetOwner(kTRUE);
      Evpool->UpdatePool(ProtonCloneArray);
    }
    else{AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ));}
  } //End of Pool updating

  //Update Event Pool at the End of the Event if Event Mixing is enabled
  if(fFillredPairTreeME){  
  if(!fRequireProton||fEventhasProton){
  if(!fRequireSigma||fEventhasSigma){
  
    if(fDebug) cout << "Updating Event Pool\n";
   
    //Get Pool from Pool Manager for given RefMult and Z Vertex values
    AliEventPool* Evpool = 0x0;
	  if(fEvPoolMgr2) Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ, (Double_t)fMeanProtonpt);
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
    else if(fMeanProtonpt>0){AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f, fMeanProtonpt = %f", fRefMultComb08, primaryVtxPosZ, fMeanProtonpt));}
  } //Require Sigma
  } //Require Proton
  } //End of Pool updating

  if(fDebug) cout << "End of Event processing. Calling PostData()\n";

/************************End of Event Processing**************************************/

    PostData(1, fOutputList); // stream the analysis results of the current event to output manager        
    PostData(2, fSigmaCandTree);
    PostData(3, fSigmaCandTreeExtra);
    PostData(4, fSigmaPairTree);
    PostData(5, fProtonTree);
    PostData(6, fSigmaPairTreeSE);         
    PostData(7, fSigmaPairTreeME); 
    PostData(8, fSigmaMEBackgroundTree);                                   
    PostData(9, fSigmaCandTreerot);         

    if(fDebug) cout << "Returning from UserExec()\n";
    return;

}//end of UserExec()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillProtonArray() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fProtonArray.clear();
  fProtonArray2.clear();
  Int_t countProton = 0;
  Double_t Protonptsum = 0;
  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

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
          if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()) isPrimaryProton = kTRUE;
          if(mcPart->IsSecondaryFromWeakDecay()) isProtonfromWeakDecay = kTRUE;
          if(mcPart->IsSecondaryFromMaterial()) isProtonfromMaterial = kTRUE;
          if(mcPart->GetMother()!=-1){
            AliAODMCParticle* ProtonMotherPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mcPart->GetMother()));
            if(ProtonMotherPart){
              if(ProtonMotherPart->GetPdgCode()==3222 || ProtonMotherPart->GetPdgCode()==-3222){ 
                isProtonfromSigma = kTRUE;
                if(ProtonMotherPart->IsPrimary()||ProtonMotherPart->IsPhysicalPrimary()) isPrimarySigma = kTRUE;
              }//Is really Sigma
              if(ProtonMotherPart->GetPdgCode()==2114) isProtonfromLambda = kTRUE;
            }//MC Mother exists
          }//Proton has Mother
        }//is Proton or Anti-Proton
      }//MC Particle exists
    }//MC treatment

    FillHistogram("fHistProtonStatistics",1);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",1);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",1);

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

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

    const Double_t len = aodTrack->GetIntegratedLength();
    const Double_t tim = aodTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(aodTrack->GetTPCmomentum());
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

    Protonptsum += pt;
    countProton++;

  } // End of proton track loop

  if(countProton!=0) fMeanProtonpt = Protonptsum/countProton;
  else fMeanProtonpt = -0.0015;
  FillHistogram("fHistMeanProtonpt",fMeanProtonpt);

  if(countProton!=0) fEventhasProton = kTRUE;

return;

}//End of FillProtonArray  

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

    if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){
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

    if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){
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
      if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimProtonPt",mcPart->Pt());}
    }
    if(MCPartPDGCode == -2212) {FillHistogram("fHistMCCounter",4); 
      FillHistogram("fHistMCAntiProtonPt",mcPart->Pt()); 
      if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimAntiProtonPt",mcPart->Pt());}
    }

    if(MCPartPDGCode == 2214)  {FillHistogram("fHistMCCounter",5); FillHistogram("fHistMCDeltaPt",mcPart->Pt());}
    if(MCPartPDGCode == -2214) {FillHistogram("fHistMCCounter",6); FillHistogram("fHistMCAntiDeltaPt",mcPart->Pt());}
    
    if(MCPartPDGCode == 3222) {FillHistogram("fHistMCCounter",7); 
      FillHistogram("fHistSigmaMotherPart",1);
      FillHistogram("fHistMCSigmaPt",mcPart->Pt()); 
      FillHistogram("fHistMCSigmaOrigin",TMath::Sqrt(mcPart->Xv()*mcPart->Xv()+mcPart->Yv()*mcPart->Yv())); 
      if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimSigmaPt",mcPart->Pt()); FillHistogram("fHistSigmaMotherPart",2); continue;}
    
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
      if(mcPart->IsPrimary()||mcPart->IsPhysicalPrimary()){FillHistogram("fHistMCPrimAntiSigmaPt",mcPart->Pt()); FillHistogram("fHistAntiSigmaMotherPart",2); continue;}

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
    Double_t DCAPV = (p.Cross((CV-PV))).Mag()/p.Mag();                        //DCA to PV of Photons

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

void AliAnalysisTaskSigmaPlus::FindAddPhotons() {
        
  const AliAODVertex *vtxT3D=aodEvent->GetPrimaryVertex();  

  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  TArrayI neg(nTracks);
  TArrayI pos(nTracks);
    
  Long_t nneg=0, npos=0, nvtx=0;

  for(Int_t i=0; i<nTracks; i++) {
    
    AliAODTrack *aodTrack=dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));

    // Acceptance Cut
    if(TMath::Abs(aodTrack->Eta()) > fMaxDaughtEta) continue; 

    //Check if Track is already used by one of the "official" Finders
    Int_t nFound = fFinderVector.size();
    Bool_t isused = kFALSE;
    for(Int_t iID = 0; iID<nFound; iID++){
      if(aodTrack->GetID()==fFinderVector[iID]) isused = kTRUE; 
    }
    if(isused) continue;

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

    // Cluster-based rejection
    if(aodTrack->GetTPCNcls() < fMinTPCClustDaught) continue;
                
    // Daughter track PID using TPC (only for PCM!)
    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kElectron);
    if(TMath::Abs(nSigmaTPCelectron) > fMaxNsigDaughtTPC) continue;    

    Float_t  d = -999., dz = -999.;
    aodTrack->GetImpactParameters(d,dz);
    Float_t  dcapv = TMath::Sqrt(d*d+dz*dz);

    FillHistogram("fHistPairDCAtoPV",dcapv);

    //Select on single-track to PV DCA here, do not call that O(N^2)
    if(TMath::Abs(d)<fMinDCADaughtPV) continue;

    if(aodTrack->Charge() < 0.) neg[nneg++]=i;
    if(aodTrack->Charge() > 0.) pos[npos++]=i;
  }    

  PairIndexArray.clear();

  for (Int_t i=0; i<nneg; i++) {
    Long_t nidx=neg[i];
    AliAODTrack *ntrk=dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(nidx));
    if(!ntrk) continue;
        
    for (Int_t k=0; k<npos; k++) {
      Int_t pidx=pos[k];
      AliAODTrack *ptrk=dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(pidx));
      if(!ptrk) continue;
            
      FillHistogram("fHistAddV0Statistics",1); //number of considered pairs
      FillHistogram("fHistAddV0StatisticsMC",1);
      FillHistogram("fHistAddV0StatisticsSigmaMC",1);
            
      Bool_t   isReallyPhoton = kFALSE;
      Bool_t   isPhotonfromSigma = kFALSE;

      // AOD MC treatment
      if(isMonteCarlo){     
        AliAODMCParticle* V0Part = NULL;
        AliAODMCParticle* NPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ntrk->GetLabel())));
        AliAODMCParticle* PPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(ptrk->GetLabel())));
        if(NPart&&PPart){
          if(NPart->GetMother()==PPart->GetMother()&&NPart->GetMother()!=-1){
            V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NPart->GetMother()));    
            if(V0Part){if(V0Part->GetPdgCode()==22){ 
              isReallyPhoton = kTRUE;
                    
              AliAODMCParticle* V0Mother = NULL;
              if(V0Part->GetMother()!=-1) V0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part->GetMother())));
              if(V0Mother){
                if(TMath::Abs(V0Mother->GetPdgCode())==3222) FillHistogram("fHistSigmaCounter",6);
                if(V0Mother->GetPdgCode()==111){
                        
                AliAODMCParticle* Pi0Mother = NULL;
                if(V0Mother->GetMother()!=-1) Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Mother->GetMother())));
                if(TMath::Abs(Pi0Mother?Pi0Mother->GetPdgCode():0)==3222) isPhotonfromSigma = kTRUE;
              
              }}//Mother of Photon exists and is a Pi0
            }}//Mother exists and is a Photon
          }//Both Tracks have a common Mother
        }//Both Tracks have matched MC Particle
      }//End of isMonteCarlo

      //OTF not available for this pair
      FillHistogram("fHistAddV0Statistics",2); //number of pairs not found by the Finders
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",2);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",2);

      //Create AliExternalTrackParams 
      AliExternalTrackParam nt, pt;
      nt.CopyFromVTrack(ntrk); //Copy properties from tracks
      pt.CopyFromVTrack(ptrk);
      AliExternalTrackParam *ntp=&nt, *ptp=&pt;
      Double_t xn, xp, dca;
                        
      //Re-propagate to closest position to the primary vertex
      Double_t dztemp[2], covartemp[3];
      //Safety margin: 250 cm 
      ntp->PropagateToDCA(vtxT3D,Bz,250,dztemp,covartemp);
      ptp->PropagateToDCA(vtxT3D,Bz,250,dztemp,covartemp);

      dca=nt.GetDCA(&pt,Bz,xn,xp);

      FillHistogram("fHistPairDCA",dca);
      if(isReallyPhoton) FillHistogram("fHistPairDCAMC",dca);

      if (dca > fMaxDCADaught) continue;
            
      FillHistogram("fHistAddV0Statistics",3); //passed dca cut
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",3);      
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",3);      

      nt.PropagateTo(xn,Bz);
      pt.PropagateTo(xp,Bz);
            
      AliESDv0 vertex(nt,nidx,pt,pidx);
            
      //Experimental: refit V0
      vertex.Refit();
                        
      Double_t x=vertex.Xv(), y=vertex.Yv();
      Double_t r2=x*x + y*y;
      Double_t r=TMath::Sqrt(r2);

      FillHistogram("fHistPairRadius",r);
      if(isReallyPhoton) FillHistogram("fHistPairRadiusMC",r);
      
      if (r < fMinV0Radius) continue;
      if (r > fMaxV0Radius) continue;
            
      FillHistogram("fHistAddV0Statistics",4); //passed radius cut
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",4);      
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",4);      

      Float_t cpa=vertex.GetV0CosineOfPointingAngle(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);

      FillHistogram("fHistPairCPA",cpa);
      if(isReallyPhoton) FillHistogram("fHistPairCPAMC",cpa);
      if(isPhotonfromSigma) FillHistogram("fHistPairCPAMCSigma",cpa);  

      //Simple cosine cut (no pt dependence for now)
      if (cpa < fMinV0CPA) continue;
            
      FillHistogram("fHistAddV0Statistics",5); //passed cosPA
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",5);

      //Now: Photon selection. Use same cuts as for V0s from Finders

      //Get reconstructed cartesian momentum
      Double_t Momvec[3];
      TVector3 vecN, vecP;
      vertex.GetNPxPyPz(Momvec[0],Momvec[1],Momvec[2]);
      vecN.SetXYZ(Momvec[0],Momvec[1],Momvec[2]); //negative daughter
      vertex.GetPPxPyPz(Momvec[0],Momvec[1],Momvec[2]);
      vecP.SetXYZ(Momvec[0],Momvec[1],Momvec[2]); //positive daughter 

      // Reconstruct photon with TLorentzVector
      TLorentzVector electron, positron, photon;
      electron.SetXYZM(vecN(0),vecN(1),vecN(2),cElectronMass);
      positron.SetXYZM(vecP(0),vecP(1),vecP(2),cElectronMass);
      photon=electron+positron;

      //Calculating DCA of Photons to PV 
      TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);          //Prim. Vertex
      TVector3 CV(vertex.Xv(),vertex.Yv(),vertex.Zv());                   //Conv. Vertex
      TVector3 p(photon.Px(),photon.Py(),photon.Pz()); //Momentum vectors of the photons
      Double_t DCAPV = (p.Cross((CV-PV))).Mag()/p.Mag();      //DCA to PV of Photons

      // Calculate photon invariant mass with TL
      Double_t photonmass = photon.M();

      Double_t thetaPos = ptrk->Theta();
      Double_t thetaNeg = ntrk->Theta();
      Double_t alpha    = vertex.AlphaV0();
      Double_t qt       = vertex.PtArmV0();
      Double_t totangle = TMath::Abs(vecP.Angle(vecN));
  
      // Angle calculation
      Double_t deltatheta = thetaPos - thetaNeg;

      FillHistogram("fHistPairArmPod",alpha,qt);
      FillHistogram("fHistPairOpenAngle",totangle);
      FillHistogram("fHistPairDeltaTheta",deltatheta);
      FillHistogram("fHistPairInvMass",photonmass);
      FillHistogram("fHistPairPhotonDCAPV",DCAPV);

      if(isReallyPhoton){
        FillHistogram("fHistPairptMC",photon.Perp());
        FillHistogram("fHistPairArmPodMC",alpha,qt);
        FillHistogram("fHistPairOpenAngleMC",totangle);
        FillHistogram("fHistPairDeltaThetaMC",deltatheta);
        FillHistogram("fHistPairInvMassMC",photonmass);
        FillHistogram("fHistPairPhotonDCAPVMC",DCAPV);
      }

      if(isPhotonfromSigma) FillHistogram("fHistPairPhotonDCAPVMCSigma",DCAPV);

      // Armenteros-Podolanski Cuts
      if(TMath::Abs(alpha) > fMaxalpha) continue;
      FillHistogram("fHistAddV0Statistics",6); //passed alpha
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",6);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",6);

      if(TMath::Abs(qt) > fMaxqt) continue;
      FillHistogram("fHistAddV0Statistics",7); //passed qt
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",7);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",7);

      // Angle Cut
      if(TMath::Abs(totangle) > fMaxopenangle) continue;
      FillHistogram("fHistAddV0Statistics",8); //passed openangle
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",8);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",8);

      if(TMath::Abs(deltatheta) > fMaxdeltatheta) continue;
      FillHistogram("fHistAddV0Statistics",9); //passed deltatheta
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",9);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",9);

      // Inv. Mass Cut
      if(photonmass > fMaxphotonmass) continue;
      FillHistogram("fHistAddV0Statistics",10); //passed inv mass 
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",10);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",10);

      if(isReallyPhoton) FillHistogram("fHistPairptwCutsMC",photon.Perp());
      FillHistogram("fHistPairptwCuts",photon.Perp());

      //Store ESD IDs of used electrons
      fV0ParticleIDArray.push_back(ntrk->GetID());
      fV0ParticleIDArray.push_back(ptrk->GetID());

      //Store the AOD! Track indices as pairs
      PairIndexArray.push_back(std::make_pair(nidx,pidx));  

    }//End of Positive Particle Track Loop
  }//End of Negative Particle Track Loop
  
  return;

}//End of FindAddPhotons() const

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillElectronArray() {
  
  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Electron and Positron Arrays and reset counters
  fElectronArray.clear();
  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  //Loop for Electron Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {

    //Initialisation of local variables
    Bool_t   isTPCElectron = kFALSE;
    Bool_t   isTOFElectron = kFALSE;
    Bool_t   isTPCHadron = kFALSE;

    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack) {
      AliWarning("No AOD Track!");
      continue;
    }

    //Check MC Truth
    Bool_t isReallyElectron = kFALSE;
    Bool_t isElectronfromSigma = kFALSE;
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        if(TMath::Abs(mcPart->GetPdgCode())==11){
          isReallyElectron = kTRUE;
          if(mcPart->GetMother()!=-1){
            AliAODMCParticle* photon = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
            if(photon){
              if(photon->GetPdgCode()==22&&photon->GetMother()!=-1){
                AliAODMCParticle* pion = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(photon->GetMother())));
                if(pion){
                  if(pion->GetPdgCode()==111&&pion->GetMother()!=-1){
                    AliAODMCParticle* sigma = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(pion->GetMother())));
                    if(sigma){
                      if(TMath::Abs(sigma->GetPdgCode())==3222) isElectronfromSigma = kTRUE;
                    }
                  }
                }
              }              
            }
          }
        }
      }
    }

    ////Check if electron is used in V0
    //Int_t nV0IDs = fV0ParticleIDArray.size();
    //Bool_t isusedV0 = kFALSE;
    //for(Int_t iID = 0; iID<nV0IDs; iID++){
    //  if(aodTrack->GetID()==fV0ParticleIDArray[iID]) isusedV0 = kTRUE; 
    //}
    //if(isusedV0) continue;

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

    Double_t eta           = aodTrack->Eta();
    Double_t phi           = aodTrack->Phi();
    Double_t p             = aodTrack->GetTPCmomentum();
    Double_t pt            = aodTrack->Pt();
    Int_t    nTPCClustElec = aodTrack->GetTPCNcls();
    Int_t    nITSClustElec = aodTrack->GetITSNcls();
    Double_t chi2          = aodTrack->GetTPCchi2();

    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kElectron);
    Double_t nSigmaTOFelectron = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kElectron);

    if(TMath::Abs(nSigmaTPCelectron)<fMaxNsigElecTPC) isTPCElectron = kTRUE;
    if(TMath::Abs(nSigmaTOFelectron)<fMaxNsigElecTOF) isTOFElectron = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kMuon))<fMinNsigHadronTPC)   isTPCHadron = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion))<fMinNsigHadronTPC)   isTPCHadron = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon))<fMinNsigHadronTPC)   isTPCHadron = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton))<fMinNsigHadronTPC) isTPCHadron = kTRUE;

    //QA Histograms
    if(isReallyElectron) FillHistogram("fHistElectronEtaPhiMC",eta,phi);

    //Acceptance Cut
    if(TMath::Abs(eta) > fMaxElecEta) continue; 

    //MC Histograms
    if(isReallyElectron){
      FillHistogram("fHistElectronChi2MC",chi2);
      FillHistogram("fHistElectronTPCClusterMC",nTPCClustElec);
      FillHistogram("fHistElectronITSClusterMC",nITSClustElec);
      FillHistogram("fHistElectronpvsNSigmaTPCMC",p,nSigmaTPCelectron);
      FillHistogram("fHistElectronpvsNSigmaTOFMC",p,nSigmaTOFelectron);
      FillHistogram("fHistElectronptMC",pt);
    }

    // Cluster Cut
    if(nTPCClustElec<fMinTPCClustElec) continue;

    //PID Histograms
    FillHistogram("fHistElectronpvsNSigmaTPC",p,nSigmaTPCelectron);
    FillHistogram("fHistElectronpvsNSigmaTOF",p,nSigmaTOFelectron);

    // TPC PID Cut
    if(!isTPCElectron) continue;
    // TOF PID/Hadron rejection
    if(!isTOFElectron&&isTPCHadron) continue;

    FillHistogram("fHistElectronpvsNSigTPCwCuts",p,nSigmaTPCelectron);

    //Histograms for Purity and Efficiency Calcultation
    if(isReallyElectron) FillHistogram("fHistElectronptwCutsMC",pt);
    FillHistogram("fHistElectronptwCuts",pt);

  	if(pt>fMaxElecpt) continue;

    // Store (Anti-)Proton candidates after selection
    fElectronArray.push_back(iTrack);

  } // End of electron track loop

return;

}//End of FillElectronArray  

//_____________________________________________________________________________

AliESDv0 AliAnalysisTaskSigmaPlus::Tracks2V0vertex(Int_t v0index) const {

  Int_t nidx = PairIndexArray[v0index].first;
  Int_t pidx = PairIndexArray[v0index].second;
  AliAODTrack *ntrk=dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(nidx));
  AliAODTrack *ptrk=dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(pidx));

  if(!ntrk||!ptrk){
    cout << "WARNING: Tracks for ESD Vertex do not exist!\nTrack Pair Array Size is: " << PairIndexArray.size() << ", Index is: " << v0index << "\n";
  }

  //Create AliExternalTrackParams 
  AliExternalTrackParam nt, pt;
  nt.CopyFromVTrack(ntrk); //Copy properties from tracks
  pt.CopyFromVTrack(ptrk);
  AliExternalTrackParam *ntp=&nt, *ptp=&pt;
  Double_t xn, xp, dca;

  const AliAODVertex *vtxT3D=aodEvent->GetPrimaryVertex();  

  //Re-propagate to closest position to the primary vertex if asked to do so
  Double_t dztemp[2], covartemp[3];
  //Safety margin: 250 -> exceedingly large... not sure this makes sense, but ok
  ntp->PropagateToDCA(vtxT3D,Bz,250,dztemp,covartemp);
  ptp->PropagateToDCA(vtxT3D,Bz,250,dztemp,covartemp);

  dca=nt.GetDCA(&pt,Bz,xn,xp);

  nt.PropagateTo(xn,Bz);
  pt.PropagateTo(xp,Bz);

  //Initizialize ESDv0 to use the refit function
  AliESDv0 vertex(nt,nidx,pt,pidx);
  vertex.Refit();

  Float_t cpa=vertex.GetV0CosineOfPointingAngle(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);
  vertex.SetDcaV0Daughters(dca);
  vertex.SetV0CosineOfPointingAngle(cpa);

  return vertex;

}

//________________________________________________________________________

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
                        if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE; 
                        MCSigmaMom.SetXYZM(Pi0Mother->Px(),Pi0Mother->Py(),Pi0Mother->Pz(),cSigmaMass);
                      }
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==2214) isReallyPi0fromDelta = kTRUE;
                    }

                  }}//Photon Mother exists and is Pi0
                }//Photon and Single electron have common mother
              }}//MC Photon exists
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
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          if(isReallySigma&&isPrimary) FillHistogram("fHistMCPrimSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          //Propagate the Proton to the proper Vertex. 3 methods implemented: AliExternalTrackParam, AliTrackerBase (with dEdx Corr.!), and KFParticle
          AliExternalTrackParam ETPProton1; 
          AliExternalTrackParam* ETPProton2 = new AliExternalTrackParam();
          ETPProton1.CopyFromVTrack(prot); //Copy properties from proton track
          ETPProton2->CopyFromVTrack(prot); 
          
          //Create VVertex for the decay vertex from KF    
          AliAODVertex *Sigmavtx3D = new AliAODVertex();  
          Sigmavtx3D->SetX(KFSigmaPlus.GetX());
          Sigmavtx3D->SetY(KFSigmaPlus.GetY());
          Sigmavtx3D->SetZ(KFSigmaPlus.GetZ());
          Double_t Sigmacovmatrix[6];
          Sigmacovmatrix[0] = KFSigmaPlus.GetCovariance(0);
          Sigmacovmatrix[1] = KFSigmaPlus.GetCovariance(1);
          Sigmacovmatrix[2] = KFSigmaPlus.GetCovariance(2);
          Sigmacovmatrix[3] = KFSigmaPlus.GetCovariance(3);
          Sigmacovmatrix[4] = KFSigmaPlus.GetCovariance(4);
          Sigmacovmatrix[5] = KFSigmaPlus.GetCovariance(5);          
          Sigmavtx3D->SetCovMatrix(Sigmacovmatrix);
          Double_t dztemp[2], covartemp[3];
          //Propagate the Proton to the decay vertex. Safety margin: 250. Might be exceedingly large, but ok
          ETPProton1.PropagateToDCA(Sigmavtx3D,Bz,250,dztemp,covartemp);
    
          AliTrackerBase* ATP = new AliTrackerBase();
          ATP->PropagateTrackTo(ETPProton2, KFSigmaPlus.GetX(), cSigmaMass, 0.01, kTRUE, 0, prot->Charge(), kTRUE, kFALSE);
          float sigposxyz[3] = {KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()};
          KFProton.TransportToPoint(sigposxyz);

          //Calculate momentum of the Sigma with the propagated momentum of the proton
          TLorentzVector trackProtonm1, trackProtonm2, trackProtonm3; 
          TLorentzVector trackSigmaplusm1, trackSigmaplusm2, trackSigmaplusm3; 
          trackProtonm1.SetXYZM(ETPProton1.Px(),ETPProton1.Py(),ETPProton1.Pz(),cProtonMass);
          trackProtonm2.SetXYZM(ETPProton2->Px(),ETPProton2->Py(),ETPProton2->Pz(),cProtonMass);
          trackProtonm3.SetXYZM(KFProton.Px(),KFProton.Py(),KFProton.Pz(),cProtonMass);
          trackSigmaplusm1 = trackPi0 + trackProtonm1;
          trackSigmaplusm2 = trackPi0 + trackProtonm2;
          trackSigmaplusm3 = trackPi0 + trackProtonm3;

          //Check the deviation from MC
          if(isReallySigma){
            FillHistogram("fHistSigmaPxResnoprop",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResnoprop",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResnoprop",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResprop",1000*(trackSigmaplusm1.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResprop",1000*(trackSigmaplusm1.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResprop",1000*(trackSigmaplusm1.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResdEdxcorr",1000*(trackSigmaplusm2.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResdEdxcorr",1000*(trackSigmaplusm2.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResdEdxcorr",1000*(trackSigmaplusm2.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxRespropKF",1000*(trackSigmaplusm3.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyRespropKF",1000*(trackSigmaplusm3.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzRespropKF",1000*(trackSigmaplusm3.Pz()-MCSigmaMom.Pz()));              
          }

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {              //Check if there are Proton Sigma Pairs
            if(q==k) continue;                          //Skip if its the same Proton
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;  //Use only Particle-Particle/Antiparticle-Antiparticle Pairs
            pairs++;

            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());   //Now calculate kstar. Quite complicated. Definition e.g. in ArXiv:2012.09806v1
            TVector3 deltapvec=sigmamomentum-protonmomentum;
            Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
            Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));
            
            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar);   //Do some preselection and save as QA
            if(prot->Charge()>0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistSigmaProtonkstar",kstar);                                    
            if(prot->Charge()<0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistAntiSigmaProtonkstar",kstar);                                    
            if(kstar<flowkstar) pairslowkstar++;
            if(kstar<fverylowkstar) pairsverylowkstar++;
            if(kstar<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          if(pi0mass>fMinCorrPi0Mass&&pi0mass<fMaxCorrPi0Mass&&SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
          else fIsGoodCandidate = kFALSE;
          fIsV01fromFinder = kTRUE;
          fIsV02fromFinder = kTRUE;
          fIsV01Onthefly = v0_1->GetOnFlyStatus();
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fHas4DiffIDs = hasdiffindices;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = EventTriggers;
          fSigMCLabel = Pi0MotherLabel;
          fSigProtonID = prot->GetID();
          fSigProtonStatus = prot->GetStatus();
          fSigEventID = fGlobalEventID;
          fSigCentrality = Centrality;
          fSigRefMultComb05 = fRefMultComb05;
          fSigRefMultComb08 = fRefMultComb08;
          fSigRefMultComb10 = fRefMultComb10;
          fSigBField = Bz;
          fInvSigMass = sigmaplusmass; 
          fSigPA = SigmaPointingAngle; 
          fSigCharge = prot->Charge(); 
          fSigPx = trackSigmaplus.Px(); 
          fSigPy = trackSigmaplus.Py(); 
          fSigPz = trackSigmaplus.Pz(); 
          fPrimVertX = primaryVtxPosX; 
          fPrimVertY = primaryVtxPosY; 
          fPrimVertZ = primaryVtxPosZ; 
          fSigDecayVertX = KFSigmaPlus.GetX(); 
          fSigDecayVertY = KFSigmaPlus.GetY(); 
          fSigDecayVertZ = KFSigmaPlus.GetZ();
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
          fProtonPxatDCA = prot->PxAtDCA(); 
          fProtonPyatDCA = prot->PyAtDCA(); 
          fProtonPzatDCA = prot->PzAtDCA(); 
          fProtonpropPx = ETPProton1.Px(); 
          fProtonpropPy = ETPProton1.Py(); 
          fProtonpropPz = ETPProton1.Pz();
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
          fnPair = pairs;
          fnPairlowkstar = pairslowkstar; 
          fnPairverylowkstar = pairsverylowkstar; 
          fnPairveryverylowkstar = pairsveryverylowkstar; 
          if(fSavePartCand) fSigmaCandTree->Fill();

          //Now apply some Selections to filter out potential Sigma Candidates.
          //If it passes the criteria, Sigma Proton pairs in SE and ME can be written out. 
          //Also all Protons in the Event can be written for offline Event-Mixing.
          if(!fSavePairs&&!fSaveAllProtons&&!fFillredPairTreeSE&&!fFillredPairTreeME) continue;
          if(pi0mass<fMinCorrPi0Mass) continue;
          if(pi0mass>fMaxCorrPi0Mass) continue;
          if(SigmaPointingAngle>fMaxCorrSigmaPA) continue;
          if(TMath::Abs(DCAxy)<fMinCorrProtonDCAxy) continue;
          if(sigmaplusmass<fMinCorrSigmaMass) continue; 
          if(sigmaplusmass>fMaxCorrSigmaMass) continue;

          if(sigmaplusmass>1.17&&sigmaplusmass<1.21) fEventhasSigma = kTRUE;

          for(Int_t q=0; q<nProton; q++) { 
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
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

            fPairProtonIsMC = kFALSE;
            fPairProtonIsPrimary = kFALSE;
            if(isMonteCarlo){
              AliAODMCParticle* PairProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(pairprot->GetLabel())));
              if(PairProtonPart){
                if(TMath::Abs(PairProtonPart->GetPdgCode())==2212){
                  fPairProtonIsMC = kTRUE;
                  if(PairProtonPart->IsPrimary()||PairProtonPart->IsPhysicalPrimary()) fPairProtonIsPrimary = kTRUE;
                }//MC Particle is a Proton
              }//MC Particle exists 
            }//End of isMonteCarlo

            pairprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

            fPairProtonPx = pairprot->Px();        
            fPairProtonPy = pairprot->Py();        
            fPairProtonPz = pairprot->Pz();
            //fPairProtonPxatDCA = pairprot->PxAtDCA();        
            //fPairProtonPyatDCA = pairprot->PyAtDCA();        
            //fPairProtonPzatDCA = pairprot->PzAtDCA();
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

            //Force TOF PID for Pair Proton if requested
            if(pairprot->P()>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

            if(fSaveAllProtons) fProtonTree->Fill();

            if(kstar>fMaxCorrkstar) continue;
            if(fPairProtonDCAtoPVxy>fMaxCorrPairProtonDCAxy) continue;
            if(fPairProtonDCAtoPVz>fMaxCorrPairProtonDCAz) continue;

            if(q!=k&&fSavePairs) fSigmaPairTree->Fill();
            if(q!=k&&fFillredPairTreeSE) fSigmaPairTreeSE->Fill();
          }        

          //Continue here if no Event Mixing is requested
          if(!fFillredPairTreeME) continue;

          //Get Pool from Pool Manager for given RefMult and Z Vertex values
        	AliEventPool* Evpool = 0x0;
	        if(fEvPoolMgr2) Evpool = fEvPoolMgr2->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ, (Double_t)fMeanProtonpt);
		      if(!Evpool){  if(fMeanProtonpt>0) AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f, fMeanProtonpt = %f", fRefMultComb08, primaryVtxPosZ, fMeanProtonpt)); continue;}
          if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << ", fMeanProtonpt = " << fMeanProtonpt << " is empty!\n";*/ continue;}

          //Get Number of Events in Pool
    			Int_t nMixEvents = Evpool->GetCurrentNEvents();
          FillHistogram("fHistPairNMixedEvents1",fRefMultComb08,primaryVtxPosZ,nMixEvents);
          FillHistogram("fHistPairNMixedEvents2",fMeanProtonpt,fRefMultComb08,nMixEvents);
          FillHistogram("fHistPairNMixedEvents3",fMeanProtonpt,primaryVtxPosZ,nMixEvents);

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

              fPairProtonIsMC = kFALSE;
              fPairProtonIsPrimary = kFALSE;

              mixprot->GetImpactParameters(fPairProtonDCAtoPVxy,fPairProtonDCAtoPVz);

              fPairProtonPx = mixprot->Px();        
              fPairProtonPy = mixprot->Py();        
              fPairProtonPz = mixprot->Pz();
              fPairProtonPxatDCA = -999;//mixprot->PxAtDCA();        
              fPairProtonPyatDCA = -999;//mixprot->PyAtDCA();        
              fPairProtonPzatDCA = -999;//mixprot->PzAtDCA();
              fPairProtonCharge = mixprot->Charge();        
              fPairProtonNSigTPC = mixprot->NumberOfSigmasTPCProton();        
              fPairProtonNSigTOF = mixprot->NumberOfSigmasTOFProton();
              fPairProtNSigTPCKaon = mixprot->NumberOfSigmasTPCKaon();        
              fPairProtNSigTOFKaon = mixprot->NumberOfSigmasTOFKaon();
              fPairProtNSigTPCPion = mixprot->NumberOfSigmasTPCPion();        
              fPairProtNSigTOFPion = mixprot->NumberOfSigmasTOFPion();
              fPairProtonCluster = mixprot->GetTPCNcls();
              fPairProtonITSCluster = mixprot->GetITSNcls();
              fPairProtonID = mixprot->GetID();
              fPairProtonStatus = mixprot->GetStatus();

              Double_t MixProtonMom = TMath::Sqrt(fPairProtonPx*fPairProtonPx+fPairProtonPy*fPairProtonPy+fPairProtonPz*fPairProtonPz);

              //Force TOF PID for Pair Proton if requested
              if(MixProtonMom>fMaxpOnlyTPCPID&&fRequireProtonTOFforPairs&&TMath::Abs(fPairProtonNSigTOF)>fMaxNsigProtTOF) continue;

              if(kstar>fMaxCorrkstar) continue;
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
  	      if(fEvPoolMgr) Evpool = fEvPoolMgr->GetEventPool((Int_t)fRefMultComb08, (Double_t)primaryVtxPosZ);
  		    if(!Evpool){AliWarning(Form("No pool found for fRefMultComb08 = %hd, primaryVtxPosZ = %f", fRefMultComb08, primaryVtxPosZ)); nMixEvents = -1;}
          else if(Evpool->GetCurrentNEvents()==0) {/*cout << "Pool for fRefMultComb08 = "<< fRefMultComb08 << ", primaryVtxPosZ = " << primaryVtxPosZ << " is empty!\n";*/ nMixEvents = -1;}

          //Get Number of Events in Pool. Number can be reduced with Setter Function to reduce Tree Size.
      		if(nMixEvents!=-1) nMixEvents = Evpool->GetCurrentNEvents();
          else nMixEvents = 0;
          FillHistogram("fHistBkgNMixedEvents",fRefMultComb08,primaryVtxPosZ,nMixEvents);

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

            // Fill the ME Sigma Trees
            if(pi0mass>fMinCorrPi0Mass&&pi0mass<fMaxCorrPi0Mass&&SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
            else fIsGoodCandidate = kFALSE;
            fIsV01fromFinder = kTRUE;
            fIsV02fromFinder = kTRUE;
            fIsV01Onthefly = v0_1->GetOnFlyStatus();
            fIsV02Onthefly = v0_2->GetOnFlyStatus();
            fHas4DiffIDs = hasdiffindices;
            fSigRunnumber = aodEvent->GetRunNumber();
            fSigTriggerMask = EventTriggers;
            fSigProtonID = mixprot->GetID();
            fSigProtonStatus = mixprot->GetStatus();
            fSigEventID = fGlobalEventID;
            fSigCentrality = Centrality;
            fSigRefMultComb05 = fRefMultComb05;
            fSigRefMultComb08 = fRefMultComb08;
            fSigRefMultComb10 = fRefMultComb10;
            fSigBField = Bz;
            fInvSigMass = sigmaplusmass; 
            fSigPA = SigmaPointingAngle; 
            fSigCharge = mixprot->Charge(); 
            fSigPx = trackSigmaplus.Px(); 
            fSigPy = trackSigmaplus.Py(); 
            fSigPz = trackSigmaplus.Pz(); 
            fPrimVertX = primaryVtxPosX; 
            fPrimVertY = primaryVtxPosY; 
            fPrimVertZ = primaryVtxPosZ; 
            fSigDecayVertX = KFSigmaPlus.GetX(); 
            fSigDecayVertY = KFSigmaPlus.GetY(); 
            fSigDecayVertZ = KFSigmaPlus.GetZ();
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
            fProtonPxatDCA = mixprot->PxAtDCA(); 
            fProtonPyatDCA = mixprot->PyAtDCA(); 
            fProtonPzatDCA = mixprot->PzAtDCA();  
            fProtonpropPx = -999; 
            fProtonpropPy = -999; 
            fProtonpropPz = -999;
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

        //Do Rotation Method for the Background if requested
        if(fSaveRotateBackground){

        //Loop over Protons again and rotate them by 180°
        for(Int_t k=0; k<nProton; k++) {

          AliAODTrack *prot;
          prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(k));
          if(!prot) continue;

          trackProton.SetXYZM(-(prot->Px()),-(prot->Py()),-(prot->Pz()),cProtonMass); //Invert momentum components
          trackSigmaplus = trackPi0 + trackProton;

          Float_t sigmaplusmass = trackSigmaplus.M();
          FillHistogram("fHistInvSigmaMassrot",sigmaplusmass);
          if(sigmaplusmass>fMaxSigmaMass) continue;   //Limit the mass range to reduce tree size

          prot->GetXYZ(trackxyz);      
          prot->GetPxPyPz(trackpxpypz);
          //for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
          trackparams[0] = (2*primaryVtxPosX)-trackxyz[0];            //Point reflection of space coordinates
          trackparams[1] = (2*primaryVtxPosY)-trackxyz[1];
          trackparams[2] = (2*primaryVtxPosZ)-trackxyz[2];
          for (Int_t q = 0; q<3;q++) {trackparams[q+3] = -trackpxpypz[q];} //Invert momentum components
          prot->GetCovarianceXYZPxPyPz(covMatrix);
          if(prot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
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

          Bool_t isReallySigma = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){if(TMath::Abs(ProtonPart->GetPdgCode())==2212){
              if(TMath::Abs(ProtonPart->GetMother())==TMath::Abs(Pi0MotherLabel)){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE;}
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);
          DCAxy = -DCAxy; DCAz = -DCAz; //Invert DCA

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

          TLorentzVector sig; sig.SetXYZM(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz(),cSigmaMass);
          Float_t Rapidity = sig.Rapidity();
          if(isReallySigma) FillHistogram("fHistMCSigmaYrot",Rapidity);
          FillHistogram("fHistSigmaYrot",Rapidity);
          if(TMath::Abs(Rapidity)>fMaxSigmaY) continue;  //Cut on the rapidity

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPArot",SigmaPointingAngle);
          if(isReallySigma&&isPrimary) FillHistogram("fHistMCPrimSigmaPArot",SigmaPointingAngle);
          FillHistogram("fHistSigmaPArot",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          if(pi0mass>fMinCorrPi0Mass&&pi0mass<fMaxCorrPi0Mass&&SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
          else fIsGoodCandidate = kFALSE;
          fIsV01fromFinder = kTRUE;
          fIsV02fromFinder = kTRUE;
          fIsV01Onthefly = v0_1->GetOnFlyStatus();
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fHas4DiffIDs = hasdiffindices;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = EventTriggers;
          fSigMCLabel = Pi0MotherLabel;
          fSigProtonID = prot->GetID();
          fSigProtonStatus = prot->GetStatus();
          fSigEventID = fGlobalEventID;
          fSigCentrality = Centrality;
          fSigRefMultComb05 = fRefMultComb05;
          fSigRefMultComb08 = fRefMultComb08;
          fSigRefMultComb10 = fRefMultComb10;
          fSigBField = Bz;
          fInvSigMass = sigmaplusmass; 
          fSigPA = SigmaPointingAngle; 
          fSigCharge = prot->Charge(); 
          fSigPx = trackSigmaplus.Px(); 
          fSigPy = trackSigmaplus.Py(); 
          fSigPz = trackSigmaplus.Pz(); 
          fPrimVertX = primaryVtxPosX; 
          fPrimVertY = primaryVtxPosY; 
          fPrimVertZ = primaryVtxPosZ; 
          fSigDecayVertX = KFSigmaPlus.GetX(); 
          fSigDecayVertY = KFSigmaPlus.GetY(); 
          fSigDecayVertZ = KFSigmaPlus.GetZ();
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
          fProtonPx = -(prot->Px()); 
          fProtonPy = -(prot->Py()); 
          fProtonPz = -(prot->Pz());
          fProtonPxatDCA = -(prot->PxAtDCA()); 
          fProtonPyatDCA = -(prot->PyAtDCA()); 
          fProtonPzatDCA = -(prot->PzAtDCA());   
          fProtonpropPx = -999; 
          fProtonpropPy = -999; 
          fProtonpropPz = -999;
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
          fnPair = -999;
          fnPairlowkstar = -999; 
          fnPairverylowkstar = -999; 
          fnPairveryverylowkstar = -999; 
          fSigmaCandTreerot->Fill();

        }//End of rotated Proton Loop

        }//End of Rotation

      /************************End of Sigma+ reconstruction*****************************/        

    }//End of Photon 2 loop
  }//End of Photon 1 loop
  
return;

} //End of ReconstructParticles()

//________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesOff() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TLorentzVector trackPhoton1, trackPhoton2, trackPi0, trackProton, trackSigmaplus;
  KFParticle KFElectron1, KFElectron2, KFPositron1, KFPositron2, KFProton; 

  const Int_t nConvPhoton = fConvPhotonArray.size();
  const Int_t nOffPhoton = PairIndexArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nOffPhoton; i++) {

    AliESDv0 v0_1 = (AliESDv0)Tracks2V0vertex(i);

    for(Int_t j=i; j<nConvPhoton; j++) {

      AliAODv0 *v0_2 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(j));
      if(!v0_2) continue;

      // Get daughter tracks      
      AliAODTrack* track1 = (AliAODTrack*)aodEvent->GetTrack(v0_1.GetNindex());
      AliAODTrack* track2 = (AliAODTrack*)aodEvent->GetTrack(v0_1.GetPindex());
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
                        if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE; 
                        MCSigmaMom.SetXYZM(Pi0Mother->Px(),Pi0Mother->Py(),Pi0Mother->Pz(),cSigmaMass);
                      }
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==2214) isReallyPi0fromDelta = kTRUE;
                    }

                  }}//Photon Mother exists and is Pi0
                }//Photon and Single electron have common mother
              }}//MC Photon exists
            }//Electrons have common mother
          }//Daughters are all electrons
        }//V0 Daughters exist
      }//End of isMonteCarlo

      trackPhoton1.SetXYZM(v0_1.Px(),v0_1.Py(),v0_1.Pz(),0);
      trackPhoton2.SetXYZM(v0_2->Px(),v0_2->Py(),v0_2->Pz(),0);
      trackPi0 = trackPhoton1 + trackPhoton2;
      Double_t pi0mass = trackPi0.M();

      //KF Pi0 calculations
      // Set up KFParticle
      trackparams[0] = v0_1.Xv();
      trackparams[1] = v0_1.Yv();
      trackparams[2] = v0_1.Zv();
      v0_1.GetPPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
      if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron1.Create(trackparams,covMatrix,1,cElectronMass);

      // Repeat for all other particles
      trackparams[0] = v0_1.Xv();
      trackparams[1] = v0_1.Yv();
      trackparams[2] = v0_1.Zv();
      v0_1.GetNPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
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

      //Reconstruct the Photons with two Methods: Standard and special gamma reconstruction
      KFParticle KFPhoton1(KFElectron1,KFPositron1);
      KFParticle KFPhoton2(KFElectron2,KFPositron2);

      //Transport Photons to Conversion Points
      KFPhoton1.TransportToDecayVertex();
      KFPhoton2.TransportToDecayVertex();

      Double_t PhotPhotDCA = TMath::Abs(KFPhoton1.GetDistanceFromParticle(KFPhoton2));

      //Calculating DCA of Photons to PV 
      TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);              //Prim. Vertex
      TVector3 CV1(v0_1.Xv(),v0_1.Yv(),v0_1.Zv());                          //Conv. Vertices
      TVector3 CV2(v0_2->DecayVertexV0X(),v0_2->DecayVertexV0Y(),v0_2->DecayVertexV0Z());
      TVector3 p1(v0_1.Px(),v0_1.Py(),v0_1.Pz());          //Momentum vectors of the photons
      TVector3 p2(v0_2->Px(),v0_2->Py(),v0_2->Pz());   
      Double_t DCAPV1 = (p1.Cross((CV1-PV))).Mag()/p1.Mag(); //DCA to PV of Photons
      Double_t DCAPV2 = (p2.Cross((CV2-PV))).Mag()/p2.Mag(); //using line-point distance equation

      Bool_t hasdiffindices = kFALSE;
      if(track1->GetID()!=track3->GetID()&&track1->GetID()!=track4->GetID()&&track2->GetID()!=track3->GetID()&&track2->GetID()!=track4->GetID()) hasdiffindices = kTRUE;

      FillHistogram("fHistGammaPairStatsOneadd",1);
      if(isReallyPi0) FillHistogram("fHistGammaPairStatsOneadd",2);
      if(isSameGamma) FillHistogram("fHistGammaPairStatsOneadd",3);
      if(hasdiffindices) FillHistogram("fHistGammaPairStatsOneadd",4);
      if(isReallyPi0&&hasdiffindices) FillHistogram("fHistGammaPairStatsOneadd",5);
      if(isSameGamma&&hasdiffindices) FillHistogram("fHistGammaPairStatsOneadd",6);

  	  FillHistogram("fHistGammaPairInvMassOneAdd",pi0mass);
  	  if(hasdiffindices) FillHistogram("fHistGammaPairInvMassOneAdd2",pi0mass);
  	  else FillHistogram("fHistGammaPairInvMassOneAdd3",pi0mass);
      if(isReallyPi0){FillHistogram("fHistGammaPairInvMassOneAddMC",pi0mass);}

      if(!hasdiffindices&&fCleanAutoCorr) continue;  //Discard Photon Pairs that share at least one Track. Is most probably auto-correlation!
      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

      // Save Photon quality
      Int_t nMinTPCClustDaught = track1->GetTPCNcls();
      if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();
      if(track3->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track3->GetTPCNcls();
      if(track4->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track4->GetTPCNcls();

      Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron));

      Double_t MaxAlpha = TMath::Abs(v0_1.AlphaV0());
      if(TMath::Abs(v0_2->AlphaV0())>MaxAlpha) MaxAlpha = TMath::Abs(v0_2->AlphaV0());

      Double_t MaxQt = TMath::Abs(v0_1.PtArmV0());
      if(TMath::Abs(v0_2->PtArmV0())>MaxQt) MaxQt = TMath::Abs(v0_2->PtArmV0());

      //Get reconstructed cartesian momentum
      Double_t Momvector[3];
      TVector3 vecN1, vecP1;
      v0_1.GetNPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecN1.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //negative daughter
      v0_1.GetPPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecP1.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //positive daughter 

      Double_t MaxOpenAngle = v0_2->OpenAngleV0();
      if(vecP1.Angle(vecN1)>MaxOpenAngle) MaxOpenAngle = vecP1.Angle(vecN1);

      TLorentzVector tl1, tl2, tlp;
      tl1.SetXYZM(v0_2->MomNegX(),v0_2->MomNegY(),v0_2->MomNegZ(),cElectronMass);
      tl2.SetXYZM(v0_2->MomPosX(),v0_2->MomPosY(),v0_2->MomPosZ(),cElectronMass);
      tlp=tl1+tl2;

      Double_t Maxphotonmass = tlp.M();

      tl1.SetXYZM(vecN1[0],vecN1[1],vecN1[2],cElectronMass);
      tl2.SetXYZM(vecP1[0],vecP1[1],vecP1[2],cElectronMass);
      tlp=tl1+tl2;

      if(tlp.M()>Maxphotonmass) Maxphotonmass = tlp.M();

      Double_t nMaxPhotchi2 = v0_2->Chi2V0();

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
          if(!KFPi0CD.CheckDaughter(KFPhoton2)) continue;

          KFParticle KFPi0(KFPhoton1,KFPhoton2);
          KFPi0.TransportToDecayVertex();

          Double_t KFPi0DCAPV = TMath::Sqrt(KFPi0.GetX()*KFPi0.GetX()+KFPi0.GetY()*KFPi0.GetY());  
          if(isReallyPi0){ 
            FillHistogram("fHistPi0VertexMC",MCPi0DCAPV);
            FillHistogram("fHistPi0VertexvsMC",KFPi0DCAPV-MCPi0DCAPV);
          }

          KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
          KFSigmaPlusCD.AddDaughter(KFProton);
          if(!KFSigmaPlusCD.CheckDaughter(KFPi0)) continue;
    
          KFParticle KFSigmaPlus(KFProton,KFPi0);
          KFSigmaPlus.TransportToDecayVertex();
          Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));

          Bool_t isReallySigma = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){if(TMath::Abs(ProtonPart->GetPdgCode())==2212){
              if(TMath::Abs(ProtonPart->GetMother())==TMath::Abs(Pi0MotherLabel)){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",2);}
                else if (isReallyPi0fromDelta) FillHistogram("fHistMCdeltamass",sigmaplusmass); 
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          if(isReallySigma&&isPrimary) FillHistogram("fHistMCPrimSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          //Propagate the Proton to the proper Vertex. 3 methods implemented: AliExternalTrackParam, AliTrackerBase (with dEdx Corr.!), and KFParticle
          AliExternalTrackParam ETPProton1; 
          AliExternalTrackParam* ETPProton2 = new AliExternalTrackParam();
          ETPProton1.CopyFromVTrack(prot); //Copy properties from proton track
          ETPProton2->CopyFromVTrack(prot); 

          //Create VVertex for the decay vertex from KF    
          AliAODVertex *Sigmavtx3D = new AliAODVertex();  
          Sigmavtx3D->SetX(KFSigmaPlus.GetX());
          Sigmavtx3D->SetY(KFSigmaPlus.GetY());
          Sigmavtx3D->SetZ(KFSigmaPlus.GetZ());
          Double_t Sigmacovmatrix[6];
          Sigmacovmatrix[0] = KFSigmaPlus.GetCovariance(0);
          Sigmacovmatrix[1] = KFSigmaPlus.GetCovariance(1);
          Sigmacovmatrix[2] = KFSigmaPlus.GetCovariance(2);
          Sigmacovmatrix[3] = KFSigmaPlus.GetCovariance(3);
          Sigmacovmatrix[4] = KFSigmaPlus.GetCovariance(4);
          Sigmacovmatrix[5] = KFSigmaPlus.GetCovariance(5);          
          Sigmavtx3D->SetCovMatrix(Sigmacovmatrix);
          Double_t dztemp[2], covartemp[3];
          //Propagate the Proton to the decay vertex. Safety margin: 250. Might be exceedingly large, but ok
          ETPProton1.PropagateToDCA(Sigmavtx3D,Bz,250,dztemp,covartemp);

          AliTrackerBase* ATP = new AliTrackerBase();
          ATP->PropagateTrackTo(ETPProton2, KFSigmaPlus.GetX(), cSigmaMass, 0.01, kTRUE, 0, prot->Charge(), kTRUE, kFALSE);
          float sigposxyz[3] = {KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()};
          KFProton.TransportToPoint(sigposxyz);

          //Calculate momentum of the Sigma with the propagated momentum of the proton
          TLorentzVector trackProtonm1, trackProtonm2, trackProtonm3; 
          TLorentzVector trackSigmaplusm1, trackSigmaplusm2, trackSigmaplusm3; 
          trackProtonm1.SetXYZM(ETPProton1.Px(),ETPProton1.Py(),ETPProton1.Pz(),cProtonMass);
          trackProtonm2.SetXYZM(ETPProton2->Px(),ETPProton2->Py(),ETPProton2->Pz(),cProtonMass);
          trackProtonm3.SetXYZM(KFProton.Px(),KFProton.Py(),KFProton.Pz(),cProtonMass);
          trackSigmaplusm1 = trackPi0 + trackProtonm1;
          trackSigmaplusm2 = trackPi0 + trackProtonm2;
          trackSigmaplusm3 = trackPi0 + trackProtonm3;

          //Check the deviation from MC
          if(isReallySigma){
            FillHistogram("fHistSigmaPxResnoprop",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResnoprop",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResnoprop",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResprop",1000*(trackSigmaplusm1.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResprop",1000*(trackSigmaplusm1.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResprop",1000*(trackSigmaplusm1.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResdEdxcorr",1000*(trackSigmaplusm2.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResdEdxcorr",1000*(trackSigmaplusm2.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResdEdxcorr",1000*(trackSigmaplusm2.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxRespropKF",1000*(trackSigmaplusm3.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyRespropKF",1000*(trackSigmaplusm3.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzRespropKF",1000*(trackSigmaplusm3.Pz()-MCSigmaMom.Pz()));              
          }

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {
            if(q==k) continue;
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;
            pairs++;

            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 deltapvec=sigmamomentum-protonmomentum;
            Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
            Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));

            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar);
            if(prot->Charge()>0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistSigmaProtonkstar",kstar);                                    
            if(prot->Charge()<0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistAntiSigmaProtonkstar",kstar);                                    
            if(kstar<flowkstar) pairslowkstar++;
            if(kstar<fverylowkstar) pairsverylowkstar++;
            if(kstar<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          if(pi0mass>fMinCorrPi0Mass&&pi0mass<fMaxCorrPi0Mass&&SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
          else fIsGoodCandidate = kFALSE;
          fIsV01fromFinder = kFALSE;
          fIsV02fromFinder = kTRUE;
          fIsV01Onthefly = kFALSE;
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fHas4DiffIDs = hasdiffindices;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = (Int_t)aodEvent->GetTriggerMask();
          fSigMCLabel = Pi0MotherLabel;
          fSigProtonID = prot->GetID();
          fSigProtonStatus = prot->GetStatus();
          fSigEventID = fGlobalEventID;
          fSigCentrality = Centrality;
          fSigRefMultComb05 = fRefMultComb05;
          fSigRefMultComb08 = fRefMultComb08;
          fSigRefMultComb10 = fRefMultComb10;
          fSigBField = Bz;
          fInvSigMass = sigmaplusmass; 
          fSigPA = SigmaPointingAngle; 
          fSigCharge = prot->Charge(); 
          fSigPx = trackSigmaplus.Px(); 
          fSigPy = trackSigmaplus.Py(); 
          fSigPz = trackSigmaplus.Pz(); 
          fPrimVertX = primaryVtxPosX; 
          fPrimVertY = primaryVtxPosY; 
          fPrimVertZ = primaryVtxPosZ; 
          fSigDecayVertX = KFSigmaPlus.GetX(); 
          fSigDecayVertY = KFSigmaPlus.GetY(); 
          fSigDecayVertZ = KFSigmaPlus.GetZ();
          fSigDecayVertXMC = MCPi0decayX;
          fSigDecayVertYMC = MCPi0decayY;
          fSigDecayVertZMC = MCPi0decayZ;
          fSigPxMC = MCSigmaMom.Px();        
          fSigPyMC = MCSigmaMom.Py();        
          fSigPzMC = MCSigmaMom.Pz();         
          fPhoton1Radius = TMath::Sqrt(v0_1.Xv()*v0_1.Xv()+v0_1.Yv()*v0_1.Yv());
          fPhoton2Radius = TMath::Sqrt(v0_2->DecayVertexV0X()*v0_2->DecayVertexV0X()+v0_2->DecayVertexV0Y()*v0_2->DecayVertexV0Y());
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2;
          fPhotonsMinCluster   = nMinTPCClustDaught;
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
          fProtonpropPx = ETPProton1.Px(); 
          fProtonpropPy = ETPProton1.Py(); 
          fProtonpropPz = ETPProton1.Pz();
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
          fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
          fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
          fProtonNCluster = prot->GetTPCNcls();
          fProtonChi2 = prot->GetTPCchi2();    
          fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
          fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
          fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
          fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
          fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
          fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);
          fnPair = pairs;
          fnPairlowkstar = pairslowkstar; 
          fnPairverylowkstar = pairsverylowkstar; 
          fnPairveryverylowkstar = pairsveryverylowkstar; 
          fSigmaCandTree->Fill();

        }//End of Proton loop

      /************************End of Sigma+ reconstruction*****************************/        
        
      }//End of Photon 2 loop
    }//End of Photon 1 loop

return;

} //End of ReconstructParticlesOff()

//________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesOff2() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TLorentzVector trackPhoton1, trackPhoton2, trackPi0, trackProton, trackSigmaplus;
  KFParticle KFElectron1, KFElectron2, KFPositron1, KFPositron2, KFProton; 

  const Int_t nOffPhoton = PairIndexArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nOffPhoton-1; i++) {

    AliESDv0 v0_1 = (AliESDv0)Tracks2V0vertex(i);

    for(Int_t j=i+1; j<nOffPhoton; j++) {

      AliESDv0 v0_2 = (AliESDv0)Tracks2V0vertex(j);

      // Get daughter tracks      
      AliAODTrack* track1 = (AliAODTrack*)aodEvent->GetTrack(v0_1.GetNindex());
      AliAODTrack* track2 = (AliAODTrack*)aodEvent->GetTrack(v0_1.GetPindex());
      AliAODTrack* track3 = (AliAODTrack*)aodEvent->GetTrack(v0_2.GetNindex());
      AliAODTrack* track4 = (AliAODTrack*)aodEvent->GetTrack(v0_2.GetPindex());

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
                        if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE;
                        MCSigmaMom.SetXYZM(Pi0Mother->Px(),Pi0Mother->Py(),Pi0Mother->Pz(),cSigmaMass);  
                      }
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==2214) isReallyPi0fromDelta = kTRUE;
                    }

                  }}//Photon Mother exists and is Pi0
                }//Photon and Single electron have common mother
              }}//MC Photon exists
            }//Electrons have common mother
          }//Daughters are all electrons
        }//V0 Daughters exist
      }//End of isMonteCarlo

      trackPhoton1.SetXYZM(v0_1.Px(),v0_1.Py(),v0_1.Pz(),0);
      trackPhoton2.SetXYZM(v0_2.Px(),v0_2.Py(),v0_2.Pz(),0);
      trackPi0 = trackPhoton1 + trackPhoton2;
      Double_t pi0mass = trackPi0.M();

      //KF Pi0 calculations
      // Set up KFParticle
      trackparams[0] = v0_1.Xv();
      trackparams[1] = v0_1.Yv();
      trackparams[2] = v0_1.Zv();
      v0_1.GetPPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
      if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron1.Create(trackparams,covMatrix,1,cElectronMass);

      // Repeat for all other particles
      trackparams[0] = v0_1.Xv();
      trackparams[1] = v0_1.Yv();
      trackparams[2] = v0_1.Zv();
      v0_1.GetNPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
      if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron1.Create(trackparams,covMatrix,-1,cElectronMass);

      trackparams[0] = v0_2.Xv();
      trackparams[1] = v0_2.Yv();
      trackparams[2] = v0_2.Zv();
      v0_2.GetPPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
      if(track3->Charge()>0) track3->GetCovarianceXYZPxPyPz(covMatrix);
      else track4->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron2.Create(trackparams,covMatrix,1,cElectronMass);

      trackparams[0] = v0_2.Xv();
      trackparams[1] = v0_2.Yv();
      trackparams[2] = v0_2.Zv();
      v0_2.GetNPxPyPz(trackparams[3],trackparams[4],trackparams[5]);
      if(track3->Charge()<0) track3->GetCovarianceXYZPxPyPz(covMatrix);
      else track4->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron2.Create(trackparams,covMatrix,-1,cElectronMass);

      //Reconstruct the Photons with two Methods: Standard and special gamma reconstruction
      KFParticle KFPhoton1(KFElectron1,KFPositron1);
      KFParticle KFPhoton2(KFElectron2,KFPositron2);

      //Transport Photons to Conversion Points
      KFPhoton1.TransportToDecayVertex();
      KFPhoton2.TransportToDecayVertex();

      Double_t PhotPhotDCA = TMath::Abs(KFPhoton1.GetDistanceFromParticle(KFPhoton2));

      //Calculating DCA of Photons to PV 
      TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);              //Prim. Vertex
      TVector3 CV1(v0_1.Xv(),v0_1.Yv(),v0_1.Zv());                          //Conv. Vertices
      TVector3 CV2(v0_2.Xv(),v0_2.Yv(),v0_2.Zv());
      TVector3 p1(v0_1.Px(),v0_1.Py(),v0_1.Pz());          //Momentum vectors of the photons
      TVector3 p2(v0_2.Px(),v0_2.Py(),v0_2.Pz());   
      Double_t DCAPV1 = (p1.Cross((CV1-PV))).Mag()/p1.Mag(); //DCA to PV of Photons
      Double_t DCAPV2 = (p2.Cross((CV2-PV))).Mag()/p2.Mag(); //using line-point distance equation

      Bool_t hasdiffindices = kFALSE;
      if(track1->GetID()!=track3->GetID()&&track1->GetID()!=track4->GetID()&&track2->GetID()!=track3->GetID()&&track2->GetID()!=track4->GetID()) hasdiffindices = kTRUE;

      FillHistogram("fHistGammaPairStatsOnlyadd",1);
      if(isReallyPi0) FillHistogram("fHistGammaPairStatsOnlyadd",2);
      if(isSameGamma) FillHistogram("fHistGammaPairStatsOnlyadd",3);
      if(hasdiffindices) FillHistogram("fHistGammaPairStatsOnlyadd",4);
      if(isReallyPi0&&hasdiffindices) FillHistogram("fHistGammaPairStatsOnlyadd",5);
      if(isSameGamma&&hasdiffindices) FillHistogram("fHistGammaPairStatsOnlyadd",6);

  	  FillHistogram("fHistGammaPairInvMassOnlyAdd",pi0mass);
  	  if(hasdiffindices) FillHistogram("fHistGammaPairInvMassOnlyAdd2",pi0mass);
  	  else FillHistogram("fHistGammaPairInvMassOnlyAdd3",pi0mass);
      if(isReallyPi0){FillHistogram("fHistGammaPairInvMassOnlyAddMC",pi0mass);}

      if(!hasdiffindices&&fCleanAutoCorr) continue;  //Discard Photon Pairs that share at least one Track. Is most probably auto-correlation!
      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

      // Save Photon quality
      Int_t nMinTPCClustDaught = track1->GetTPCNcls();
      if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();
      if(track3->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track3->GetTPCNcls();
      if(track4->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track4->GetTPCNcls();

      Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track3,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track4,AliPID::kElectron));

      Double_t MaxAlpha = TMath::Abs(v0_1.AlphaV0());
      if(TMath::Abs(v0_2.AlphaV0())>MaxAlpha) MaxAlpha = TMath::Abs(v0_2.AlphaV0());

      Double_t MaxQt = TMath::Abs(v0_1.PtArmV0());
      if(TMath::Abs(v0_2.PtArmV0())>MaxQt) MaxQt = TMath::Abs(v0_2.PtArmV0());

      //Get reconstructed cartesian momentum
      Double_t Momvector[3];
      TVector3 vecN1, vecP1, vecN2, vecP2;
      v0_1.GetNPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecN1.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //negative daughter
      v0_1.GetPPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecP1.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //positive daughter 
      v0_2.GetNPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecN2.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //negative daughter
      v0_2.GetPPxPyPz(Momvector[0],Momvector[1],Momvector[2]);
      vecP2.SetXYZ(Momvector[0],Momvector[1],Momvector[2]); //positive daughter 

      Double_t MaxOpenAngle = vecP1.Angle(vecN1);
      if(vecP2.Angle(vecN2)>MaxOpenAngle) MaxOpenAngle = vecP2.Angle(vecN2);

      TLorentzVector tl1, tl2, tlp;
      tl1.SetXYZM(vecN1[0],vecN1[1],vecN1[2],cElectronMass);
      tl2.SetXYZM(vecP1[0],vecP1[1],vecP1[2],cElectronMass);
      tlp=tl1+tl2;

      Double_t Maxphotonmass = tlp.M();

      tl1.SetXYZM(vecN2[0],vecN2[1],vecN2[2],cElectronMass);
      tl2.SetXYZM(vecP2[0],vecP2[1],vecP2[2],cElectronMass);
      tlp=tl1+tl2;

      if(tlp.M()>Maxphotonmass) Maxphotonmass = tlp.M();

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
          if(!KFPi0CD.CheckDaughter(KFPhoton2)) continue;

          KFParticle KFPi0(KFPhoton1,KFPhoton2);
          KFPi0.TransportToDecayVertex();

          Double_t KFPi0DCAPV = TMath::Sqrt(KFPi0.GetX()*KFPi0.GetX()+KFPi0.GetY()*KFPi0.GetY());  
          if(isReallyPi0){ 
            FillHistogram("fHistPi0VertexMC",MCPi0DCAPV);
            FillHistogram("fHistPi0VertexvsMC",KFPi0DCAPV-MCPi0DCAPV);
          }

          KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
          KFSigmaPlusCD.AddDaughter(KFProton);
          if(!KFSigmaPlusCD.CheckDaughter(KFPi0)) continue;
    
          KFParticle KFSigmaPlus(KFProton,KFPi0);
          KFSigmaPlus.TransportToDecayVertex();
          Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));

          Bool_t isReallySigma = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){if(TMath::Abs(ProtonPart->GetPdgCode())==2212){
              if(TMath::Abs(ProtonPart->GetMother())==TMath::Abs(Pi0MotherLabel)){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",3);}
                else if (isReallyPi0fromDelta) FillHistogram("fHistMCdeltamass",sigmaplusmass); 
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
          if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
          if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          if(isReallySigma&&isPrimary) FillHistogram("fHistMCPrimSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          //Propagate the Proton to the proper Vertex. 3 methods implemented: AliExternalTrackParam, AliTrackerBase (with dEdx Corr.!), and KFParticle
          AliExternalTrackParam ETPProton1; 
          AliExternalTrackParam* ETPProton2 = new AliExternalTrackParam();
          ETPProton1.CopyFromVTrack(prot); //Copy properties from proton track
          ETPProton2->CopyFromVTrack(prot); 

          //Create VVertex for the decay vertex from KF    
          AliAODVertex *Sigmavtx3D = new AliAODVertex();  
          Sigmavtx3D->SetX(KFSigmaPlus.GetX());
          Sigmavtx3D->SetY(KFSigmaPlus.GetY());
          Sigmavtx3D->SetZ(KFSigmaPlus.GetZ());
          Double_t Sigmacovmatrix[6];
          Sigmacovmatrix[0] = KFSigmaPlus.GetCovariance(0);
          Sigmacovmatrix[1] = KFSigmaPlus.GetCovariance(1);
          Sigmacovmatrix[2] = KFSigmaPlus.GetCovariance(2);
          Sigmacovmatrix[3] = KFSigmaPlus.GetCovariance(3);
          Sigmacovmatrix[4] = KFSigmaPlus.GetCovariance(4);
          Sigmacovmatrix[5] = KFSigmaPlus.GetCovariance(5);          
          Sigmavtx3D->SetCovMatrix(Sigmacovmatrix);
          Double_t dztemp[2], covartemp[3];
          //Propagate the Proton to the decay vertex. Safety margin: 250. Might be exceedingly large, but ok
          ETPProton1.PropagateToDCA(Sigmavtx3D,Bz,250,dztemp,covartemp);

          AliTrackerBase* ATP = new AliTrackerBase();
          ATP->PropagateTrackTo(ETPProton2, KFSigmaPlus.GetX(), cSigmaMass, 0.01, kTRUE, 0, prot->Charge(), kTRUE, kFALSE);
          float sigposxyz[3] = {KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()};
          KFProton.TransportToPoint(sigposxyz);

          //Calculate momentum of the Sigma with the propagated momentum of the proton
          TLorentzVector trackProtonm1, trackProtonm2, trackProtonm3; 
          TLorentzVector trackSigmaplusm1, trackSigmaplusm2, trackSigmaplusm3; 
          trackProtonm1.SetXYZM(ETPProton1.Px(),ETPProton1.Py(),ETPProton1.Pz(),cProtonMass);
          trackProtonm2.SetXYZM(ETPProton2->Px(),ETPProton2->Py(),ETPProton2->Pz(),cProtonMass);
          trackProtonm3.SetXYZM(KFProton.Px(),KFProton.Py(),KFProton.Pz(),cProtonMass);
          trackSigmaplusm1 = trackPi0 + trackProtonm1;
          trackSigmaplusm2 = trackPi0 + trackProtonm2;
          trackSigmaplusm3 = trackPi0 + trackProtonm3;

          //Check the deviation from MC
          if(isReallySigma){
            FillHistogram("fHistSigmaPxResnoprop",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResnoprop",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResnoprop",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResprop",1000*(trackSigmaplusm1.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResprop",1000*(trackSigmaplusm1.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResprop",1000*(trackSigmaplusm1.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxResdEdxcorr",1000*(trackSigmaplusm2.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyResdEdxcorr",1000*(trackSigmaplusm2.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzResdEdxcorr",1000*(trackSigmaplusm2.Pz()-MCSigmaMom.Pz()));              
            FillHistogram("fHistSigmaPxRespropKF",1000*(trackSigmaplusm3.Px()-MCSigmaMom.Px()));            
            FillHistogram("fHistSigmaPyRespropKF",1000*(trackSigmaplusm3.Py()-MCSigmaMom.Py()));            
            FillHistogram("fHistSigmaPzRespropKF",1000*(trackSigmaplusm3.Pz()-MCSigmaMom.Pz()));              
          }

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {
            if(q==k) continue;
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;
            pairs++;

            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 deltapvec=sigmamomentum-protonmomentum;
            Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
            Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
            Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
            Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
            Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));

            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar);
            if(prot->Charge()>0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistSigmaProtonkstar",kstar);                                    
            if(prot->Charge()<0&&pi0mass>0.1&&pi0mass<0.16&&SigmaPointingAngle<0.04&&TMath::Abs(DCAxy)>0.01&&sigmaplusmass>1.17&&sigmaplusmass<1.2) FillHistogram("fHistAntiSigmaProtonkstar",kstar);                                    
            if(kstar<flowkstar) pairslowkstar++;
            if(kstar<fverylowkstar) pairsverylowkstar++;
            if(kstar<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          if(pi0mass>fMinCorrPi0Mass&&pi0mass<fMaxCorrPi0Mass&&SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
          else fIsGoodCandidate = kFALSE;
          fIsV01fromFinder = kFALSE;
          fIsV02fromFinder = kFALSE;
          fIsV01Onthefly = kFALSE;
          fIsV02Onthefly = kFALSE;
          fHas4DiffIDs = hasdiffindices;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = EventTriggers;
          fSigMCLabel = Pi0MotherLabel;
          fSigProtonID = prot->GetID();
          fSigProtonStatus = prot->GetStatus();
          fSigEventID = fGlobalEventID;
          fSigCentrality = Centrality;
          fSigRefMultComb05 = fRefMultComb05;
          fSigRefMultComb08 = fRefMultComb08;
          fSigRefMultComb10 = fRefMultComb10;
          fSigBField = Bz;
          fInvSigMass = sigmaplusmass; 
          fSigPA = SigmaPointingAngle; 
          fSigCharge = prot->Charge(); 
          fSigPx = trackSigmaplus.Px(); 
          fSigPy = trackSigmaplus.Py(); 
          fSigPz = trackSigmaplus.Pz(); 
          fPrimVertX = primaryVtxPosX; 
          fPrimVertY = primaryVtxPosY; 
          fPrimVertZ = primaryVtxPosZ; 
          fSigDecayVertX = KFSigmaPlus.GetX(); 
          fSigDecayVertY = KFSigmaPlus.GetY(); 
          fSigDecayVertZ = KFSigmaPlus.GetZ();
          fSigDecayVertXMC = MCPi0decayX;
          fSigDecayVertYMC = MCPi0decayY;
          fSigDecayVertZMC = MCPi0decayZ;
          fSigPxMC = MCSigmaMom.Px();        
          fSigPyMC = MCSigmaMom.Py();        
          fSigPzMC = MCSigmaMom.Pz();        
          fPhoton1Radius = TMath::Sqrt(v0_1.Xv()*v0_1.Xv()+v0_1.Yv()*v0_1.Yv());
          fPhoton2Radius = TMath::Sqrt(v0_2.Xv()*v0_2.Xv()+v0_2.Yv()*v0_2.Yv());
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2;
          fPhotonsMinCluster   = nMinTPCClustDaught;
          fPhotonsMaxalpha     = MaxAlpha;
          fPhotonsMaxqt        = MaxQt;
          fPhotonsMaxOpenAngle = MaxOpenAngle;
          fPhotonsMaxinvmass   = Maxphotonmass; 
          fPhotonsMaxNSigTPC   = nMaxNsigTPCDaught;
          fPhotonsMaxChi2      = -999;
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
          fProtonpropPx = ETPProton1.Px();
          fProtonpropPy = ETPProton1.Py();
          fProtonpropPz = ETPProton1.Pz();
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
          fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
          fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
          fProtonNCluster = prot->GetTPCNcls();
          fProtonChi2 = prot->GetTPCchi2();  
          fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
          fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
          fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
          fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
          fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
          fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);
          fnPair = pairs;
          fnPairlowkstar = pairslowkstar; 
          fnPairverylowkstar = pairsverylowkstar; 
          fnPairveryverylowkstar = pairsveryverylowkstar; 
          fSigmaCandTree->Fill();

        }//End of Proton loop

      /************************End of Sigma+ reconstruction*****************************/        
        
      }//End of Photon 2 loop
    }//End of Photon 1 loop

return;

} //End of ReconstructParticlesOff2()

//________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesOneGamma() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TLorentzVector trackPhoton, trackProton, trackSigmaplus;
  KFParticle KFElectron, KFPositron, KFProton; 

  const Int_t nConvPhoton = fConvPhotonArray.size();
  const Int_t nProton = fProtonArray.size();

  for(Int_t i=0; i<nConvPhoton; i++) {

    AliAODv0* v0 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray.at(i));
    if(!v0) continue;

    // Get daughter tracks      
    AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));

    if(!track1 || !track2) {
    AliWarning("ERROR: Could not retrieve all AOD tracks in Pi0 reconstruction!");
    continue;
    }

    trackPhoton.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),0);

    Bool_t isPhotonfromSigma = kFALSE;
    Bool_t isPrimary = kFALSE;
    Int_t SigmaLabel = -1;
    TLorentzVector MCSigmaMom;

    AliAODMCParticle* PhotDaught1 = NULL;
    AliAODMCParticle* PhotDaught2 = NULL;
    AliAODMCParticle* PhotonPart = NULL;    
    AliAODMCParticle* SigmaPart = NULL;
    
    if(isMonteCarlo){
      PhotDaught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
      PhotDaught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
      
      if(PhotDaught1&&PhotDaught2){
        if(PhotDaught1->GetMother()!=-1&&PhotDaught2->GetMother()!=-1&&PhotDaught1->GetMother()==PhotDaught2->GetMother()){
          PhotonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(PhotDaught1->GetMother())));        
        }  
      }    
              
      if(PhotonPart){
        if(PhotonPart->GetMother()!=-1&&PhotonPart->GetPdgCode()==22){
          SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(PhotonPart->GetMother())));        
        }
      }

      if(SigmaPart){
        if(TMath::Abs(SigmaPart->GetPdgCode())==3222){
          isPhotonfromSigma = kTRUE;
          SigmaLabel = SigmaPart->GetLabel();
          if(SigmaPart->IsPrimary()||SigmaPart->IsPhysicalPrimary()) isPrimary = kTRUE;
          MCSigmaMom.SetXYZM(SigmaPart->Px(),SigmaPart->Py(),SigmaPart->Pz(),cSigmaMass);
          FillHistogram("fHistSigmaCounter",4);
        }
      }
    }//End of isMonteCarlo  

    for(Int_t j=0; j<nProton; j++) {

      AliAODTrack *prot;
      prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(j));
      if(!prot) continue;

      trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);
      trackSigmaplus = trackPhoton + trackProton; 
      Float_t sigmaplusmass = trackSigmaplus.M();
      if(sigmaplusmass>fMaxSigmaMass) continue;  //Coarse Mass Cut 

      Bool_t isReallySigma = kFALSE;
      if(isMonteCarlo&&isPhotonfromSigma){
        AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
        if(ProtonPart){
          if(ProtonPart->GetMother()!=-1&&TMath::Abs(ProtonPart->GetPdgCode())==2212){
            if(TMath::Abs(ProtonPart->GetMother())==TMath::Abs(SigmaPart->GetLabel())) isReallySigma = kTRUE;
          }
        }
      }//End of isMonteCarlo

      //KF calculations
      // Set up KFParticle
      trackparams[0] = v0->DecayVertexV0X();
      trackparams[1] = v0->DecayVertexV0Y();
      trackparams[2] = v0->DecayVertexV0Z();
      trackparams[3] = v0->MomPosX();
      trackparams[4] = v0->MomPosY();
      trackparams[5] = v0->MomPosZ();
      if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron.Create(trackparams,covMatrix,1,cElectronMass);

      // Repeat for all other particles
      trackparams[0] = v0->DecayVertexV0X();
      trackparams[1] = v0->DecayVertexV0Y();
      trackparams[2] = v0->DecayVertexV0Z();
      trackparams[3] = v0->MomNegX();
      trackparams[4] = v0->MomNegY();
      trackparams[5] = v0->MomNegZ();
      if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron.Create(trackparams,covMatrix,-1,cElectronMass);

      //Reconstruct the Photon
      KFParticle KFPhoton(KFElectron,KFPositron); //Reconstruct with default method in KF

      //Transport Photon to Conversion Point
      KFPhoton.TransportToDecayVertex();

      prot->GetXYZ(trackxyz);      
      prot->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      prot->GetCovarianceXYZPxPyPz(covMatrix);
      if(prot->Charge()>0) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
      else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);
      
      Double_t ProtPhotDCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPhoton));

      //Calculating DCA of Photon to PV 
      TVector3 PV(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);                   //Prim. Vertex
      TVector3 CV(v0->DecayVertexV0X(),v0->DecayVertexV0Y(),v0->DecayVertexV0Z()); //Conv. Vertex
      TVector3 p(v0->Px(),v0->Py(),v0->Pz());                   //Momentum vectors of the photons
      Double_t DCAPV = (p.Cross((CV-PV))).Mag()/p.Mag();               //DCA to PV of Photons

      KFParticleCD KFSigmaPlusCD; //Check Daughters to avoid floating point exceptions. See .h-file
      KFSigmaPlusCD.AddDaughter(KFProton);
      if(!KFSigmaPlusCD.CheckDaughter(KFPhoton)) continue;
      KFParticle KFSigmaPlus(KFProton,KFPhoton);
      
      KFSigmaPlus.TransportToDecayVertex();

      Float_t  DCAxy = -999., DCAz = -999.;
      prot->GetImpactParameters(DCAxy,DCAz);

      if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
      if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&DCAxy!=-999) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAz&&DCAz!=-999) continue;
      if(TMath::Abs(DCAxy)>fMaxProtonDCAxy&&fRequireDCACut) continue; 
      if(TMath::Abs(DCAz)>fMaxProtonDCAz&&fRequireDCACut) continue;

      TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
      TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
      Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
      FillHistogram("fHistOneGammaSigmaPA",SigmaPointingAngle);
      if(isReallySigma) FillHistogram("fHistMCOneGammaSigmaPA",SigmaPointingAngle);
      if(isReallySigma&&isPrimary) FillHistogram("fHistMCPrimOneGammaSigmaPA",SigmaPointingAngle);
      if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

      // Save Photon quality
      Int_t nMinTPCClustDaught = track1->GetTPCNcls();
      if(track2->GetTPCNcls()<nMinTPCClustDaught) nMinTPCClustDaught = track2->GetTPCNcls();
      Double_t nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1,AliPID::kElectron));
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron))>nMaxNsigTPCDaught) nMaxNsigTPCDaught = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron));
      Double_t MaxAlpha = TMath::Abs(v0->AlphaV0());
      Double_t MaxQt = TMath::Abs(v0->PtArmV0());
      Double_t MaxOpenAngle = v0->OpenAngleV0();
      TLorentzVector tl1, tl2, tlp;
      tl1.SetXYZM(v0->MomNegX(),v0->MomNegY(),v0->MomNegZ(),cElectronMass);
      tl2.SetXYZM(v0->MomPosX(),v0->MomPosY(),v0->MomPosZ(),cElectronMass);
      tlp=tl1+tl2;
      Double_t Maxphotonmass = tlp.M();
      Double_t nMaxPhotchi2 = v0->Chi2V0();

      //Propagate the Proton to the proper Vertex. 3 methods implemented: AliExternalTrackParam, AliTrackerBase (with dEdx Corr.!), and KFParticle
      AliExternalTrackParam ETPProton1; 
      AliExternalTrackParam* ETPProton2 = new AliExternalTrackParam();
      ETPProton1.CopyFromVTrack(prot); //Copy properties from proton track
      ETPProton2->CopyFromVTrack(prot); 

      //Create VVertex for the decay vertex from KF    
      AliAODVertex *Sigmavtx3D = new AliAODVertex();  
      Sigmavtx3D->SetX(KFSigmaPlus.GetX());
      Sigmavtx3D->SetY(KFSigmaPlus.GetY());
      Sigmavtx3D->SetZ(KFSigmaPlus.GetZ());
      Double_t Sigmacovmatrix[6];
      Sigmacovmatrix[0] = KFSigmaPlus.GetCovariance(0);
      Sigmacovmatrix[1] = KFSigmaPlus.GetCovariance(1);
      Sigmacovmatrix[2] = KFSigmaPlus.GetCovariance(2);
      Sigmacovmatrix[3] = KFSigmaPlus.GetCovariance(3);
      Sigmacovmatrix[4] = KFSigmaPlus.GetCovariance(4);
      Sigmacovmatrix[5] = KFSigmaPlus.GetCovariance(5);          
      Sigmavtx3D->SetCovMatrix(Sigmacovmatrix);
      Double_t dztemp[2], covartemp[3];
      //Propagate the Proton to the decay vertex. Safety margin: 250. Might be exceedingly large, but ok
      ETPProton1.PropagateToDCA(Sigmavtx3D,Bz,250,dztemp,covartemp);

      AliTrackerBase* ATP = new AliTrackerBase();
      ATP->PropagateTrackTo(ETPProton2, KFSigmaPlus.GetX(), cSigmaMass, 0.01, kTRUE, 0, prot->Charge(), kTRUE, kFALSE);
      float sigposxyz[3] = {KFSigmaPlus.GetX(),KFSigmaPlus.GetY(),KFSigmaPlus.GetZ()};
      KFProton.TransportToPoint(sigposxyz);

      //Calculate momentum of the Sigma with the propagated momentum of the proton
      TLorentzVector trackProtonm1, trackProtonm2, trackProtonm3; 
      TLorentzVector trackSigmaplusm1, trackSigmaplusm2, trackSigmaplusm3; 
      trackProtonm1.SetXYZM(ETPProton1.Px(),ETPProton1.Py(),ETPProton1.Pz(),cProtonMass);
      trackProtonm2.SetXYZM(ETPProton2->Px(),ETPProton2->Py(),ETPProton2->Pz(),cProtonMass);
      trackProtonm3.SetXYZM(KFProton.Px(),KFProton.Py(),KFProton.Pz(),cProtonMass);
      trackSigmaplusm1 = trackPhoton + trackProtonm1;
      trackSigmaplusm2 = trackPhoton + trackProtonm2;
      trackSigmaplusm3 = trackPhoton + trackProtonm3;

      //Check the deviation from MC
      if(isReallySigma){
        FillHistogram("fHistSigmaPxResnoprop",1000*(trackSigmaplus.Px()-MCSigmaMom.Px()));            
        FillHistogram("fHistSigmaPyResnoprop",1000*(trackSigmaplus.Py()-MCSigmaMom.Py()));            
        FillHistogram("fHistSigmaPzResnoprop",1000*(trackSigmaplus.Pz()-MCSigmaMom.Pz()));              
        FillHistogram("fHistSigmaPxResprop",1000*(trackSigmaplusm1.Px()-MCSigmaMom.Px()));            
        FillHistogram("fHistSigmaPyResprop",1000*(trackSigmaplusm1.Py()-MCSigmaMom.Py()));            
        FillHistogram("fHistSigmaPzResprop",1000*(trackSigmaplusm1.Pz()-MCSigmaMom.Pz()));              
        FillHistogram("fHistSigmaPxResdEdxcorr",1000*(trackSigmaplusm2.Px()-MCSigmaMom.Px()));            
        FillHistogram("fHistSigmaPyResdEdxcorr",1000*(trackSigmaplusm2.Py()-MCSigmaMom.Py()));            
        FillHistogram("fHistSigmaPzResdEdxcorr",1000*(trackSigmaplusm2.Pz()-MCSigmaMom.Pz()));              
        FillHistogram("fHistSigmaPxRespropKF",1000*(trackSigmaplusm3.Px()-MCSigmaMom.Px()));            
        FillHistogram("fHistSigmaPyRespropKF",1000*(trackSigmaplusm3.Py()-MCSigmaMom.Py()));            
        FillHistogram("fHistSigmaPzRespropKF",1000*(trackSigmaplusm3.Pz()-MCSigmaMom.Pz()));              
      }

      Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
      for(Int_t q=0; q<nProton; q++) {
        if(q==j) continue;
        AliAODTrack *pairprot;
        pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
        if(!pairprot) continue;
        if(pairprot->Charge()!=prot->Charge()) continue;
        pairs++;

        TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
        TVector3 deltapvec=sigmamomentum-protonmomentum;
        Double_t SigmaE = TMath::Sqrt(cSigmaMass*cSigmaMass+sigmamomentum.Mag()*sigmamomentum.Mag());
        Double_t ProtonE = TMath::Sqrt(cProtonMass*cProtonMass+protonmomentum.Mag()*protonmomentum.Mag());
        Double_t qinv2 = deltapvec.Mag2() - (SigmaE-ProtonE)*(SigmaE-ProtonE);
        Double_t vara = (qinv2+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass)/2;
        Double_t kstar = TMath::Sqrt((vara*vara-cProtonMass*cProtonMass*cSigmaMass*cSigmaMass)/(2*vara+cProtonMass*cProtonMass+cSigmaMass*cSigmaMass));

        if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar);
        if(kstar<flowkstar) pairslowkstar++;
        if(kstar<fverylowkstar) pairsverylowkstar++;
        if(kstar<fveryverylowkstar) pairsveryverylowkstar++;
      }        

      // Fill the Extra Sigma Candidate Trees
      fIsMCSigma = kFALSE; 
      if(isReallySigma) fIsMCSigma = kTRUE;
      fIsMCPrimary = kFALSE;
      if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
      if(SigmaPointingAngle<fMaxCorrSigmaPA&&TMath::Abs(DCAxy)>fMinCorrProtonDCAxy&&sigmaplusmass>fMinCorrSigmaMass&&sigmaplusmass<fMaxCorrSigmaMass) fIsGoodCandidate = kTRUE;
      else fIsGoodCandidate = kFALSE;
      fIsV0fromFinder = kTRUE;
      fIsV0Onthefly = v0->GetOnFlyStatus();
      fSigRunnumber = aodEvent->GetRunNumber();
      fSigTriggerMask = EventTriggers;
      fSigMCLabel = SigmaLabel;
      fSigProtonID = prot->GetID();
      fSigProtonStatus = prot->GetStatus();
      fSigEventID = fGlobalEventID;
      fSigCentrality = Centrality;
      fSigRefMultComb05 = fRefMultComb05;
      fSigRefMultComb08 = fRefMultComb08;
      fSigRefMultComb10 = fRefMultComb10;
      fSigBField = Bz;
      fInvSigMass = sigmaplusmass; 
      fSigPA = SigmaPointingAngle; 
      fSigCharge = prot->Charge(); 
      fSigPx = trackSigmaplus.Px(); 
      fSigPy = trackSigmaplus.Py(); 
      fSigPz = trackSigmaplus.Pz(); 
      fPrimVertX = primaryVtxPosX; 
      fPrimVertY = primaryVtxPosY; 
      fPrimVertZ = primaryVtxPosZ; 
      fSigDecayVertX = KFSigmaPlus.GetX(); 
      fSigDecayVertY = KFSigmaPlus.GetY(); 
      fSigDecayVertZ = KFSigmaPlus.GetZ();
      fSigDecayVertXMC = -999; //Not implemented
      fSigDecayVertYMC = -999;
      fSigDecayVertZMC = -999;
      fSigPxMC = MCSigmaMom.Px();        
      fSigPyMC = MCSigmaMom.Py();        
      fSigPzMC = MCSigmaMom.Pz();        
      fPhotonRadius = TMath::Sqrt(v0->DecayVertexV0X()*v0->DecayVertexV0X()+v0->DecayVertexV0Y()*v0->DecayVertexV0Y()); 
      fPhotonDCAPV = DCAPV;
      fPhotonsMinCluster   = nMinTPCClustDaught;
      fPhotonsMaxalpha     = MaxAlpha;
      fPhotonsMaxqt        = MaxQt;
      fPhotonsMaxOpenAngle = MaxOpenAngle;
      fPhotonsMaxinvmass   = Maxphotonmass;
      fPhotonsMaxNSigTPC   = nMaxNsigTPCDaught;
      fPhotonsMaxChi2      = nMaxPhotchi2;
      fPhotonPx = trackPhoton.Px(); 
      fPhotonPy = trackPhoton.Py(); 
      fPhotonPz = trackPhoton.Pz(); 
      fExtPhotProtDCA = ProtPhotDCA;
      fProtonPx = prot->Px(); 
      fProtonPy = prot->Py(); 
      fProtonPz = prot->Pz();
      fProtonpropPx = ETPProton1.Px();
      fProtonpropPy = ETPProton1.Py();
      fProtonpropPz = ETPProton1.Pz();
      fProtonDCAtoPVxy = DCAxy; 
      fProtonDCAtoPVz = DCAz; 
      fProtonNSigTPC = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kProton);
      fProtonNSigTOF = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kProton);
      fProtonNCluster = prot->GetTPCNcls();
      fProtonChi2 = prot->GetTPCchi2();  
      fProtonNSigTPCPion = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kPion);
      fProtonNSigTPCKaon = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kKaon);
      fProtonNSigTPCElec = fPIDResponse->NumberOfSigmasTPC(prot,AliPID::kElectron);
      fProtonNSigTOFPion = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kPion);
      fProtonNSigTOFKaon = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kKaon);
      fProtonNSigTOFElec = fPIDResponse->NumberOfSigmasTOF(prot,AliPID::kElectron);
      fnPair = pairs;
      fnPairlowkstar = pairslowkstar; 
      fnPairverylowkstar = pairsverylowkstar; 
      fnPairveryverylowkstar = pairsveryverylowkstar; 
      fSigmaCandTreeExtra->Fill();

    }//End of Proton Loop
  }//End of Photon Loop
  
return;

} //End of ReconstructParticles4Part()

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
