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

#include <thread>         // std::this_thread::sleep_for()
#include <chrono>         // std::chrono::seconds()

class TTree;
class TParticle;
class TVector3;
class AliAODVertex;
class AliAODv0;

class AliAnalysisTaskSigmaPlus;    // analysis class

using std::cout;
using std::map;            
using std::make_pair;            
using std::pair;            

ClassImp(AliAnalysisTaskSigmaPlus) // classimp: necessary for root

AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus() : AliAnalysisTaskSE(), 
fOutputList(0), aodEvent(0x0), mcEvent(0x0),        
AODMCTrackArray(0x0), fPIDResponse(0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaCandTreeExtra(0x0),
cElectronMass(0), cProtonMass(0), cPionMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(0), primaryVtxPosY(0), primaryVtxPosZ(0), nTracks(0), Centrality(0),
fDebug(kFALSE),

fProcessProtons(kTRUE), 
fProcessMCParticles(kTRUE), 
fProcessV0s(kTRUE), 
fProcessAddPhoton(kTRUE), 
fProcessElectrons(kTRUE), 
fProcessReco(kTRUE), 
fProcessRecoOff(kTRUE), 
fProcessRecoOff2(kTRUE), 
fProcessOneGamma(kFALSE), 

fMaxProtEta(0.9),    
fMinTPCClustProt(70),
fMaxNsigProtTPC(3),
fMaxNsigProtTOF(3),  
fMaxpOnlyTPCPID(0.75),
fMinProtpt(0), 
fMaxProtpt(15),   

fMaxMCEta(0.9),

fMaxDaughtEta(0.9),    
fMinTPCClustDaught(60),
fMaxNsigDaughtTPC(3),
fMaxalpha(0.9),
fMaxqt(0.04),
fMaxopenangle(0.2),
fMaxdeltatheta(0.1),
fMinV0CPA(0.85),
fMinV0Radius(3), 
fMaxV0Radius(200),
fMaxphotonmass(0.05),

fMinDCADaughtPV(0.1),
fMaxDCADaught(5),

fMaxElecEta(0.9),    
fMinTPCClustElec(40),
fMaxNsigElecTPC(3),
fMinNsigHadronTPC(0),
fMaxNsigElecTOF(3),  
fMaxElecpt(5), 

fRejectGammaAutoCorr(0.02),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.2),  
fMaxSigmaMass(1.7),
fMinProtonDCAxy(0),
fMinProtonDCAz(0), 
flowkstar(0.6),
fverylowkstar(0.3),
fveryverylowkstar(0.1),

fIsMCSigma(kFALSE),
fIsMCPrimary(kFALSE),
fIsV01fromFinder(kFALSE),
fIsV02fromFinder(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(999),
fSigCentrality(-999),
fSigBField(-999),
fInvSigMass(-999),
fSigPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
fPhoton1Radius(-999),
fPhoton2Radius(-999),
fPhoton1DCAPV(-999),
fPhoton2DCAPV(-999),
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
fProtonDCAtoPVxy(-999),
fProtonDCAtoPVz(-999),
fProtonPi0DCA(-999),
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
fExtPhotProtDCA(999)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus(const char* name) : AliAnalysisTaskSE(name),
fOutputList(0), aodEvent(0x0), mcEvent(0x0),        
AODMCTrackArray(0x0), fPIDResponse(0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), fSigmaCandTreeExtra(0x0),
cElectronMass(0), cProtonMass(0), cPionMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(0), primaryVtxPosY(0), primaryVtxPosZ(0), nTracks(0), Centrality(0),
fDebug(kFALSE),

fProcessProtons(kTRUE), 
fProcessMCParticles(kTRUE), 
fProcessV0s(kTRUE), 
fProcessAddPhoton(kTRUE), 
fProcessElectrons(kTRUE), 
fProcessReco(kTRUE), 
fProcessRecoOff(kTRUE), 
fProcessRecoOff2(kTRUE), 
fProcessOneGamma(kFALSE), 

fMaxProtEta(0.9),    
fMinTPCClustProt(70),
fMaxNsigProtTPC(3),
fMaxNsigProtTOF(3),  
fMaxpOnlyTPCPID(0.75),
fMinProtpt(0), 
fMaxProtpt(15),   

fMaxMCEta(0.9),

fMaxDaughtEta(0.9),
fMinTPCClustDaught(60),
fMaxNsigDaughtTPC(3),
fMaxalpha(0.9),
fMaxqt(0.04),
fMaxopenangle(0.2),
fMaxdeltatheta(0.1),
fMinV0CPA(0.85),
fMinV0Radius(3), 
fMaxV0Radius(200),
fMaxphotonmass(0.05),

fMinDCADaughtPV(0.1),
fMaxDCADaught(5),

fMaxElecEta(0.9),    
fMinTPCClustElec(40),
fMaxNsigElecTPC(3),
fMinNsigHadronTPC(0),
fMaxNsigElecTOF(3),  
fMaxElecpt(5), 

fRejectGammaAutoCorr(0.02),
fMinPi0Mass(0.06), 
fMaxPi0Mass(0.19),  
fMaxSigmaPA(0.2),  
fMaxSigmaMass(1.7),
fMinProtonDCAxy(0),
fMinProtonDCAz(0), 
flowkstar(0.6),
fverylowkstar(0.3),
fveryverylowkstar(0.1),

fIsMCSigma(kFALSE),
fIsMCPrimary(kFALSE),
fIsV01fromFinder(kFALSE),
fIsV02fromFinder(kFALSE),
fIsV01Onthefly(kFALSE),
fIsV02Onthefly(kFALSE),
fSigRunnumber(-999),
fSigTriggerMask(999),
fSigCentrality(-999),
fSigBField(-999),
fInvSigMass(-999),
fSigPA(-999),
fSigCharge(-999),
fSigPx(-999),
fSigPy(-999),
fSigPz(-999),
fPrimVertX(-999),
fPrimVertY(-999),
fPrimVertZ(-999),
fSigDecayVertX(-999),
fSigDecayVertY(-999),
fSigDecayVertZ(-999),
fPhoton1Radius(-999),
fPhoton2Radius(-999),
fPhoton1DCAPV(-999),
fPhoton2DCAPV(-999),
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
fProtonDCAtoPVxy(-999),
fProtonDCAtoPVz(-999),
fProtonPi0DCA(-999),
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
fExtPhotProtDCA(999)
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
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::~AliAnalysisTaskSigmaPlus()
{
    // destructor
    // delete objects from memory at the end of the task
    if(fOutputList) delete fOutputList;     

    if(fSigmaCandTree) delete fSigmaCandTree;

    if(fSigmaCandTreeExtra) delete fSigmaCandTreeExtra;
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

    // Create TTree of Sigma Candidates
    fSigmaCandTree = new TTree("fSigmaCandTree","Tree of Sigma Candidates");
    fSigmaCandTree->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaCandTree->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaCandTree->Branch("fIsV01fromFinder",&fIsV01fromFinder,"fIsV01fromFinder/O");
    fSigmaCandTree->Branch("fIsV02fromFinder",&fIsV02fromFinder,"fIsV02fromFinder/O");
    fSigmaCandTree->Branch("fIsV01Onthefly",&fIsV01Onthefly,"fIsV01Onthefly/O");
    fSigmaCandTree->Branch("fIsV02Onthefly",&fIsV02Onthefly,"fIsV02Onthefly/O");
    fSigmaCandTree->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTree->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/g");
    fSigmaCandTree->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
    fSigmaCandTree->Branch("fSigBField",&fSigBField,"fSigBField/F");
    fSigmaCandTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
    fSigmaCandTree->Branch("fSigPA",&fSigPA,"fSigPA/F");
    fSigmaCandTree->Branch("fSigCharge",&fSigCharge,"fSigCharge/F");
    fSigmaCandTree->Branch("fSigPx",&fSigPx,"fSigPx/F");
    fSigmaCandTree->Branch("fSigPy",&fSigPy,"fSigPy/F");
    fSigmaCandTree->Branch("fSigPz",&fSigPz,"fSigPz/F");
    fSigmaCandTree->Branch("fPrimVertX",&fPrimVertX,"fPrimVertX/F");
    fSigmaCandTree->Branch("fPrimVertY",&fPrimVertY,"fPrimVertY/F");
    fSigmaCandTree->Branch("fPrimVertZ",&fPrimVertZ,"fPrimVertZ/F");
    fSigmaCandTree->Branch("fSigDecayVertX",&fSigDecayVertX,"fSigDecayVertX/F");
    fSigmaCandTree->Branch("fSigDecayVertY",&fSigDecayVertY,"fSigDecayVertY/F");
    fSigmaCandTree->Branch("fSigDecayVertZ",&fSigDecayVertZ,"fSigDecayVertZ/F");
    fSigmaCandTree->Branch("fPhoton1Radius",&fPhoton1Radius,"fPhoton1Radius/F");
    fSigmaCandTree->Branch("fPhoton2Radius",&fPhoton2Radius,"fPhoton2Radius/F");
    fSigmaCandTree->Branch("fPhoton1DCAPV",&fPhoton1DCAPV,"fPhoton1DCAPV/F");
    fSigmaCandTree->Branch("fPhoton2DCAPV",&fPhoton2DCAPV,"fPhoton2DCAPV/F");
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
    fSigmaCandTree->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTree->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTree->Branch("fProtonPi0DCA",&fProtonPi0DCA,"fProtonPi0DCA/F");
    fSigmaCandTree->Branch("fnPair",&fnPair,"fnPair/S");
    fSigmaCandTree->Branch("fnPairlowkstar",&fnPairlowkstar,"fnPairlowkstar/S");
    fSigmaCandTree->Branch("fnPairverylowkstar",&fnPairverylowkstar,"fnPairverylowkstar/S");
    fSigmaCandTree->Branch("fnPairveryverylowkstar",&fnPairveryverylowkstar,"fnPairveryverylowkstar/S");

    // Create TTree of Extra Sigma Candidates
    fSigmaCandTreeExtra = new TTree("fSigmaCandTreeExtra","Tree of Extra Sigma Candidates");
    fSigmaCandTreeExtra->Branch("fIsMCSigma",&fIsMCSigma,"fIsMCSigma/O");
    fSigmaCandTreeExtra->Branch("fIsMCPrimary",&fIsMCPrimary,"fIsMCPrimary/O");
    fSigmaCandTreeExtra->Branch("fIsV0fromFinder",&fIsV0fromFinder,"fIsV0fromFinder/O");
    fSigmaCandTreeExtra->Branch("fIsV0Onthefly",&fIsV0Onthefly,"fIsV0Onthefly/O");
    fSigmaCandTreeExtra->Branch("fSigRunnumber",&fSigRunnumber,"fSigRunnumber/I");
    fSigmaCandTreeExtra->Branch("fSigTriggerMask",&fSigTriggerMask,"fSigTriggerMask/g");
    fSigmaCandTreeExtra->Branch("fSigCentrality",&fSigCentrality,"fSigCentrality/F");
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
    fSigmaCandTreeExtra->Branch("fPhotonPx",&fPhotonPx,"fPhotonPx/F");
    fSigmaCandTreeExtra->Branch("fPhotonPy",&fPhotonPy,"fPhotonPy/F");
    fSigmaCandTreeExtra->Branch("fPhotonPz",&fPhotonPz,"fPhotonPz/F");
    fSigmaCandTreeExtra->Branch("fPhotonRadius",&fPhotonRadius,"fPhotonRadius/F");
    fSigmaCandTreeExtra->Branch("fExtPhotProtDCA",&fExtPhotProtDCA,"fExtPhotProtDCA/F");
    fSigmaCandTreeExtra->Branch("fProtonPx",&fProtonPx,"fProtonPx/F");
    fSigmaCandTreeExtra->Branch("fProtonPy",&fProtonPy,"fProtonPy/F");
    fSigmaCandTreeExtra->Branch("fProtonPz",&fProtonPz,"fProtonPz/F");
    fSigmaCandTreeExtra->Branch("fProtonDCAtoPVxy",&fProtonDCAtoPVxy,"fProtonDCAtoPVxy/F");
    fSigmaCandTreeExtra->Branch("fProtonDCAtoPVz",&fProtonDCAtoPVz,"fProtonDCAtoPVz/F");
    fSigmaCandTreeExtra->Branch("fnPair",&fnPair,"fnPair/S");
    fSigmaCandTreeExtra->Branch("fnPairlowkstar",&fnPairlowkstar,"fnPairlowkstar/S");
    fSigmaCandTreeExtra->Branch("fnPairverylowkstar",&fnPairverylowkstar,"fnPairverylowkstar/S");
    fSigmaCandTreeExtra->Branch("fnPairveryverylowkstar",&fnPairveryverylowkstar,"fnPairveryverylowkstar/S");

    //Save Particle Masses and other constants for later use
    cElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();      
    cProtonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    
    cPionMass     = TDatabasePDG::Instance()->GetParticle(211)->Mass();     
    cPi0Mass      = TDatabasePDG::Instance()->GetParticle(111)->Mass();     
    c = 2.99792457999999984e-02; // [cm/ps]

/**************************Histograms********************************/

    //Book Keeper for used Cuts 
    TH1D* fHistCutBookKeeper           = new TH1D("fHistCutBookKeeper", "Book Keeper for used Cuts", 38, 0.5, 38.5);

    //Event related                    
    TH1F* fHistVertexZ                 = new TH1F("fHistVertexZ", "Z Vertex Position;z [cm];Counts/mm", 400, -20, 20);
    TH1F* fHistCentrality              = new TH1F("fHistCentrality", "Centrality Percentile;Centrality [%];Counts/%", 100, 0, 100);
    TH1F* fHistEventCounter            = new TH1F("fHistEventCounter", "Event Counter", 1, 0.5, 1.5);
    TH1F* fHistTrackMultiplicity       = new TH1F("fHistTrackMultiplicity","Number of Tracks per Event",20000,0,20000);
                                       
    //Counters                         
    TH1F* fHistMCCounter               = new TH1F("fHistMCCounter", "Monte Carlo Counter", 18, 0.5, 18.5);
    TH1F* fHistV0Statistics            = new TH1F("fHistV0Statistics", "V0 Counter;Stage;Counts",14,0.5,14.5);
    TH1F* fHistV0StatisticsMC          = new TH1F("fHistV0StatisticsMC", "V0 Counter MC;Stage;Counts",14,0.5,14.5);
    TH1F* fHistV0StatisticsSigmaMC     = new TH1F("fHistV0StatisticsSigmaMC", "V0 Counter Sigma MC;Stage;Counts",14,0.5,14.5);
    TH1F* fHistProtonStatistics        = new TH1F("fHistProtonStatistics", "Proton Counter;Stage;Counts",7,0.5,7.5);
    TH1F* fHistProtonStatisticsMC      = new TH1F("fHistProtonStatisticsMC", "Proton Counter MC;Stage;Counts",7,0.5,7.5);
    TH1F* fHistProtonStatisticsSigmaMC = new TH1F("fHistProtonStatisticsSigmaMC", "Proton Counter Sigma MC;Stage;Counts",7,0.5,7.5);
    TH1F* fHistAddV0Statistics         = new TH1F("fHistAddV0Statistics", "Additional V0 Counter;Stage;Counts",11,0,11);
    TH1F* fHistAddV0StatisticsMC       = new TH1F("fHistAddV0StatisticsMC", "Additional V0 Counter MC;Stage;Counts",11,0,11);
    TH1F* fHistAddV0StatisticsSigmaMC  = new TH1F("fHistAddV0StatisticsSigmaMC", "Additional V0 Counter Sigma MC;Stage;Counts",11,0,11);
    TH1F* fHistSigmaCounter            = new TH1F("fHistSigmaCounter", "#Sigma Counter", 6, 0.5, 6.5);
                                       
    //MC Information              
    TH1F* fHistMCSigmaPt               = new TH1F("fHistMCSigmaPt","Transverse momentum of MC #Sigma^{+};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiSigmaPt           = new TH1F("fHistMCAntiSigmaPt","Transverse momentum of MC #bar#Sigma^{-};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimSigmaPt           = new TH1F("fHistMCPrimSigmaPt","Transverse momentum of primary MC #Sigma^{+};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiSigmaPt       = new TH1F("fHistMCPrimAntiSigmaPt","Transverse momentum of primary MC #bar#Sigma^{-};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaOrigin           = new TH1F("fHistMCSigmaOrigin","Origin of MC #Sigma^{+};r_{xy} [cm];Counts/mm",5000,0,500);   
    TH1F* fHistMCAntiSigmaOrigin       = new TH1F("fHistMCAntiSigmaOrigin","Origin of MC #bar#Sigma^{-};r_{xy} [cm];Counts/mm",5000,0,500);   
    TH1F* fHistMCDeltaPt               = new TH1F("fHistMCDeltaPt","Transverse momentum of MC #Delta^{+};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiDeltaPt           = new TH1F("fHistMCAntiDeltaPt","Transverse momentum of MC #bar#Delta^{-};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPi0Pt                 = new TH1F("fHistMCPi0Pt","Transverse momentum of MC #pi^{0};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCProtonPt              = new TH1F("fHistMCProtonPt","Transverse momentum of MC Protons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCAntiProtonPt          = new TH1F("fHistMCAntiProtonPt","Transverse momentum of MC Anti-Protons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimProtonPt          = new TH1F("fHistMCPrimProtonPt","Transverse momentum primary of MC Protons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPrimAntiProtonPt      = new TH1F("fHistMCPrimAntiProtonPt","Transverse momentum of primary MC Anti-Protons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCPhotonPt              = new TH1F("fHistMCPhotonPt","Transverse momentum of MC Photons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaPhotonPt         = new TH1F("fHistMCSigmaPhotonPt","Transverse momentum of MC Photons from #Sigma^{+}/#bar#Sigma^{-};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCConvPhotonPt          = new TH1F("fHistMCConvPhotonPt","Transverse momentum of converted MC Photons (R<180cm);#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCSigmaConvPhotonPt     = new TH1F("fHistMCSigmaConvPhotonPt","Transverse momentum of converted MC Photons from #Sigma^{+}/#bar#Sigma^{-} (R<180cm);#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistMCConvRadius            = new TH1F("fHistMCConvRadius","Conversion radius of MC Photons;r_{xy} [cm];Counts/mm",5000,0,500);   
    TH2F* fHistMCConvRadiusvspt        = new TH2F("fHistMCConvRadiusvspt","Conversion radius of MC Photons vs. #it{p}_{T};r_{xy} [cm];#it{p}_{T} [GeV/#it{c}]",500,0,500,500,0,5);   
    TH1F* fHistSigmaMotherPart         = new TH1F("fHistSigmaMotherPart", "#Sigma^{+} Mother Particle", 25, 0.5, 25.5); 
    TH1F* fHistAntiSigmaMotherPart     = new TH1F("fHistAntiSigmaMotherPart", "#bar#Sigma^{-} Mother Particle", 25, 0.5, 25.5); 

    //Track Quality                        
    TH2F* fHistTrackEtaPhi             = new TH2F("fHistTrackEtaPhi","#eta vs. #phi of Tracks;#eta;#phi [rad]",300,-1.5,1.5,300,0,2*TMath::Pi());                                
    TH1F* fHistTrackChi2               = new TH1F("fHistTrackChi2","#chi^{2}/NDF of Tracks;#chi^{2}/NDF",200,0,100);
    TH1F* fHistTrackTPCCluster         = new TH1F("fHistTrackTPCCluster","Number of TPC Clusters of Tracks",160,0,160);
    TH2F* fHistTrackpvsdEdx            = new TH2F("fHistTrackpvsdEdx","Momentum vs. TPC dE/dx;p [GeV/#it{c}];TPC #frac{dE}{dx} (a.u.)",500,0,10,500,0,500);
    TH2F* fHistTrackpvsbeta            = new TH2F("fHistTrackpvsbeta","Momentum vs. #beta;p [GeV/#it{c}];#beta",500,0,10,500,0,1.2);
    TH1F* fHistTrackpt                 = new TH1F("fHistTrackpt","Transverse momentum of Tracks;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);

    //Proton QA                        
    TH2F* fHistProtonEtaPhiMC          = new TH2F("fHistProtonEtaPhiMC","#eta vs. #phi of Protons MC;#eta;#phi [rad]",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistProtonChi2MC            = new TH1F("fHistProtonChi2MC","#chi^{2}/NDF of Protons MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistProtonTPCClusterMC      = new TH1F("fHistProtonTPCClusterMC","Number of TPC Clusters of Protons MC",160,0,160);
    TH2F* fHistProtonpvsNSigmaTPC      = new TH2F("fHistProtonpvsNSigmaTPC","Momentum vs N sigma TPC Proton;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTPCMC    = new TH2F("fHistProtonpvsNSigmaTPCMC","Momentum vs N sigma TPC Proton MC;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTOF      = new TH2F("fHistProtonpvsNSigmaTOF","Momentum vs N sigma TOF Proton;p [GeV/#it{c}];n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH2F* fHistProtonpvsNSigmaTOFMC    = new TH2F("fHistProtonpvsNSigmaTOFMC","Momentum vs N sigma TOF Proton MC;p [GeV/#it{c}];n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH1F* fHistProtonDCAxy             = new TH1F("fHistProtonDCAxy","DCA_{xy} of Protons;DCA_{xy} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAz              = new TH1F("fHistProtonDCAz","DCA_{z} of Protons;DCA_{z} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAxyMC           = new TH1F("fHistProtonDCAxyMC","DCA_{xy} of Protons MC;DCA_{xy} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAzMC            = new TH1F("fHistProtonDCAzMC","DCA_{z} of Protons MC;DCA_{z} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAxyMCSigma      = new TH1F("fHistProtonDCAxyMCSigma","DCA_{xy} of Protons from MC #Sigma^{+}/#bar#Sigma^{-};DCA_{xy} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAzMCSigma       = new TH1F("fHistProtonDCAzMCSigma","DCA_{z} of Protons from MC #Sigma^{+}/#bar#Sigma^{-};DCA_{z} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAxyMCPrimSig    = new TH1F("fHistProtonDCAxyMCPrimSig","DCA_{xy} of Protons from primary MC #Sigma^{+}/#bar#Sigma^{-};DCA_{xy} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonDCAzMCPrimSig     = new TH1F("fHistProtonDCAzMCPrimSig","DCA_{z} of Protons from primary MC #Sigma^{+}/#bar#Sigma^{-};DCA_{z} [cm];Counts/(0.1 mm))",1000,0,10);
    TH1F* fHistProtonptMC              = new TH1F("fHistProtonptMC","Transverse momentum of MC Protons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistProtonptwCutsMC         = new TH1F("fHistProtonptwCutsMC","Transverse momentum of MC Protons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistProtonptwCuts           = new TH1F("fHistProtonptwCuts","Transverse momentum of Protons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);   

    //V0 QA                            
    TH2F* fHistV0OnflyvsOffline        = new TH2F("fHistV0OnflyvsOffline", "Number of V0s found by On-the-fly Finder vs. Offline Finder", 201, -0.5, 200.5, 201, -0.5, 200.5);
    TH2F* fHistV0DaughtEtaPhi          = new TH2F("fHistV0DaughtEtaPhi","#eta vs. #phi of V0 Daughters;#eta;#phi [rad]",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH2F* fHistV0DaughtEtaPhiMC        = new TH2F("fHistV0DaughtEtaPhiMC","#eta vs. #phi of V0 Daughters MC;#eta;#phi [rad]",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistV0DaughtChi2            = new TH1F("fHistV0DaughtChi2","#chi^{2}/NDF of V0 Daughters;#chi^{2}/NDF",200,0,100);
    TH1F* fHistV0DaughtChi2MC          = new TH1F("fHistV0DaughtChi2MC","#chi^{2}/NDF of V0 Daughters MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistV0DaughtTPCClust        = new TH1F("fHistV0DaughtTPCClust","Number of TPC Clusters of V0 Daughters",160,0,160);
    TH1F* fHistV0DaughtTPCClustMC      = new TH1F("fHistV0DaughtTPCClustMC","Number of TPC Clusters of V0 Daughters MC",160,0,160);
    TH2F* fHistV0DaughtpvsNSigmaTPC    = new TH2F("fHistV0DaughtpvsNSigmaTPC","Momentum vs N sigma TPC Electron;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH1F* fHistV0DaughtDCAtoPV         = new TH1F("fHistV0DaughtDCAtoPV","DCA to PV of V0 Daughters;DCA [cm];Counts/mm",200,0,20);
    TH1F* fHistV0DaughtDCAtoPVMC       = new TH1F("fHistV0DaughtDCAtoPVMC","DCA to PV of V0 Daughters MC;DCA [cm];Counts/mm",200,0,20);
    TH1F* fHistV0DaughtDCA             = new TH1F("fHistV0DaughtDCA","DCA between V0 Daughters;DCA [cm];Counts/(0.5mm)",200,0,10);
    TH1F* fHistV0DaughtDCAMC           = new TH1F("fHistV0DaughtDCAMC","DCA between V0 Daughters MC;DCA [cm];Counts/(0.5mm)",200,0,10);
    TH1F* fHistV0CPA                   = new TH1F("fHistV0CPA","Cosine of Pointing Angle of V0s;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0CPAMC                 = new TH1F("fHistV0CPAMC","Cosine of Pointing Angle of V0s MC;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0CPAMCSigma            = new TH1F("fHistV0CPAMCSigma","Cosine of Pointing Angle of V0s MC (#gamma from #Sigma);cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistV0Radius                = new TH1F("fHistV0Radius","Radius of V0s;r_{xy} [cm];Counts/mm",3000,0,300);   
    TH1F* fHistV0RadiusMC              = new TH1F("fHistV0RadiusMC","Radius of V0s MC;r_{xy} [cm];Counts/mm",3000,0,300);   
    TH2F* fHistV0Position2D            = new TH2F("fHistV0Position2D","Position of V0s;x [cm];y [cm]",1000,-180,180,1000,-180,180);   
    TH2F* fHistV0Position2DMC          = new TH2F("fHistV0Position2DMC","Position of V0s MC;x [cm];y [cm]",1000,-180,180,1000,-180,180);   
    TH1F* fHistV0PhotonDCAPV           = new TH1F("fHistV0PhotonDCAPV","DCA to PV of Photons;DCA [cm];Counts/mm)",100,0,10);
    TH1F* fHistV0PhotonDCAPVMC         = new TH1F("fHistV0PhotonDCAPVMC","DCA to PV of MC Photons;DCA [cm];Counts/mm)",100,0,10);
    TH1F* fHistV0PhotonDCAPVMCSigma    = new TH1F("fHistV0PhotonDCAPVMCSigma","DCA to PV of MC Photons from #Sigma^{+}/#bar#Sigma^{-};DCA [cm];Counts/mm)",100,0,10);
    TH2F* fHistV0ArmPod                = new TH2F("fHistV0ArmPod", "Armenteros-Podolanski Plot of V0s;#alpha;q_{t} [GeV/#it{c}]",200,-1.0,1.0,200,-0.05,0.4);       
    TH2F* fHistV0ArmPodMC              = new TH2F("fHistV0ArmPodMC", "Armenteros-Podolanski Plot of V0s MC;#alpha;q_{t} [GeV/#it{c}]",200,-1.0,1.0,200,-0.05,0.4);
    TH1F* fHistV0OpenAngle             = new TH1F("fHistV0OpenAngle","Total Opening Angle of V0s;#xi [rad];Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistV0OpenAngleMC           = new TH1F("fHistV0OpenAngleMC","Total Opening Angle of V0s MC;#xi [rad];Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistV0DeltaTheta            = new TH1F("fHistV0DeltaTheta","#Delta#Theta of V0 Daughters;#Delta#Theta [rad];Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistV0DeltaThetaMC          = new TH1F("fHistV0DeltaThetaMC","#Delta#Theta of V0 Daughters MC;#Delta#Theta [rad];Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistV0InvMass               = new TH1F("fHistV0InvMass", "Invariant mass of V0s;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistV0InvMassMC             = new TH1F("fHistV0InvMassMC", "Invariant mass of V0s MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistV0ptMC                  = new TH1F("fHistV0ptMC","Transverse momentum of MC V0 Photons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistV0ptwCutsMC             = new TH1F("fHistV0ptwCutsMC","Transverse momentum of MC V0 Photons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistV0SigmaptMC             = new TH1F("fHistV0SigmaptMC","Transverse momentum of MC V0 Photons from #Sigma^{+}/#bar#Sigma^{-};#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistV0SigmaptwCutsMC        = new TH1F("fHistV0SigmaptwCutsMC","Transverse momentum of MC V0 Photons from #Sigma^{+}/#bar#Sigma^{-} with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistV0ptwCuts               = new TH1F("fHistV0ptwCuts","Transverse momentum of V0 Photons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);    

    //Additional Photon QA                    
    TH1F* fHistPairDCAtoPV             = new TH1F("fHistPairDCAtoPV","DCA to PV of Pair Daughters;DCA [cm];Counts/mm",200,0,20); 
    TH1F* fHistPairDCA                 = new TH1F("fHistPairDCA","DCA between Pair Daughters;DCA [cm];Counts/(0.5mm)",200,0,10); 
    TH1F* fHistPairDCAMC               = new TH1F("fHistPairDCAMC","DCA between Pair Daughters MC;DCA [cm];Counts/(0.5mm)",200,0,10);
    TH1F* fHistPairRadius              = new TH1F("fHistPairRadius","Radius of Pairs;r_{xy} [cm];Counts/mm",3000,0,300);      
    TH1F* fHistPairRadiusMC            = new TH1F("fHistPairRadiusMC","Radius of Pairs MC;r_{xy} [cm];Counts/mm",3000,0,300);   
    TH1F* fHistPairPhotonDCAPV         = new TH1F("fHistPairPhotonDCAPV","DCA to PV of Photons;DCA [cm];Counts/mm)",100,0,10);
    TH1F* fHistPairPhotonDCAPVMC       = new TH1F("fHistPairPhotonDCAPVMC","DCA to PV of MC Photons;DCA [cm];Counts/mm)",100,0,10);
    TH1F* fHistPairPhotonDCAPVMCSigma  = new TH1F("fHistPairPhotonDCAPVMCSigma","DCA to PV of MC Photons from #Sigma^{+}/#bar#Sigma^{-};DCA [cm];Counts/mm)",100,0,10);
    TH1F* fHistPairCPA                 = new TH1F("fHistPairCPA","Cosine of Pointing Angle of Pairs;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistPairCPAMC               = new TH1F("fHistPairCPAMC","Cosine of Pointing Angle of Pairs MC;cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH1F* fHistPairCPAMCSigma          = new TH1F("fHistPairCPAMCSigma","Cosine of Pointing Angle of Pairs MC (#gamma from #Sigma);cos(PA);Counts/0.0001",2000,0.85,1.05);
    TH2F* fHistPairArmPod              = new TH2F("fHistPairArmPod", "Armenteros-Podolanski Plot of Pairs;#alpha;q_{t} [GeV/#it{c}]",200,-1.0,1.0,200,-0.05,0.4);       
    TH2F* fHistPairArmPodMC            = new TH2F("fHistPairArmPodMC", "Armenteros-Podolanski Plot of Pairs MC;#alpha;q_{t} [GeV/#it{c}]",200,-1.0,1.0,200,-0.05,0.4);
    TH1F* fHistPairOpenAngle           = new TH1F("fHistPairOpenAngle","Total Opening Angle of Pairs;#xi [rad];Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistPairOpenAngleMC         = new TH1F("fHistPairOpenAngleMC","Total Opening Angle of Pairs MC;#xi [rad];Counts/(0.00314)",1000,0,TMath::Pi());
    TH1F* fHistPairDeltaTheta          = new TH1F("fHistPairDeltaTheta","#Delta#Theta of Pair Daughters;#Delta#Theta [rad];Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistPairDeltaThetaMC        = new TH1F("fHistPairDeltaThetaMC","#Delta#Theta of Pair Daughters MC;#Delta#Theta [rad];Counts/(0.00314)",1000,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistPairInvMass             = new TH1F("fHistPairInvMass", "Invariant mass of Pairs;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistPairInvMassMC           = new TH1F("fHistPairInvMassMC", "Invariant mass of Pairs MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 2000, 0, 2);
    TH1F* fHistPairptMC                = new TH1F("fHistPairptMC","Transverse momentum of MC Pair Photons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistPairptwCutsMC           = new TH1F("fHistPairptwCutsMC","Transverse momentum of MC Pair Photons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistPairptwCuts             = new TH1F("fHistPairptwCuts","Transverse momentum of Pair Photons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);   

    //Electron QA                      
    TH2F* fHistElectronEtaPhiMC        = new TH2F("fHistElectronEtaPhiMC","#eta vs. #phi of Electrons MC;#eta;#phi [rad]",300,-1.5,1.5,300,0,2*TMath::Pi());
    TH1F* fHistElectronChi2MC          = new TH1F("fHistElectronChi2MC","#chi^{2}/NDF of Electrons MC;#chi^{2}/NDF",200,0,100);
    TH1F* fHistElectronTPCClusterMC    = new TH1F("fHistElectronTPCClusterMC","Number of TPC Clusters of Electrons MC",160,0,160);
    TH2F* fHistElectronpvsNSigmaTPC    = new TH2F("fHistElectronpvsNSigmaTPC","Momentum vs N sigma TPC Electrons;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigTPCwCuts = new TH2F("fHistElectronpvsNSigTPCwCuts","Momentum vs N sigma TPC Electrons with Cuts;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTPCMC  = new TH2F("fHistElectronpvsNSigmaTPCMC","Momentum vs N sigma TPC Electrons MC;p [GeV/#it{c}];n_{#sigma,TPC}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTOF    = new TH2F("fHistElectronpvsNSigmaTOF","Momentum vs N sigma TOF Electrons;p [GeV/#it{c}];n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH2F* fHistElectronpvsNSigmaTOFMC  = new TH2F("fHistElectronpvsNSigmaTOFMC","Momentum vs N sigma TOF Electrons MC;p [GeV/#it{c}];n_{#sigma,TOF}",500,0,10,500,-10,10);
    TH1F* fHistElectronptMC            = new TH1F("fHistElectronptMC","Transverse momentum of MC Electrons;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);          
    TH1F* fHistElectronptwCutsMC       = new TH1F("fHistElectronptwCutsMC","Transverse momentum of MC Electrons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);
    TH1F* fHistElectronptwCuts         = new TH1F("fHistElectronptwCuts","Transverse momentum of Electrons with Cuts;#it{p}_{T} [GeV/#it{c}];Counts/(MeV/#it{c})",12000,0,12);   

    //Gamma Gamma QA                      
    TH1F* fHistGammaPairInvMass        = new TH1F("fHistGammaPairInvMass", "Invariant mass of Photon Pairs;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMC      = new TH1F("fHistGammaPairInvMassMC", "Invariant mass of Photon Pairs MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnfly   = new TH1F("fHistGammaPairInvMassOnfly", "Invariant mass of On-the-fly Photon Pairs;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassMCOnfly = new TH1F("fHistGammaPairInvMassMCOnfly", "Invariant mass of On-the-fly Photon Pairs MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAdd  = new TH1F("fHistGammaPairInvMassOneAdd", "Invariant mass of Photon Pairs, one additional Photon;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOneAddMC= new TH1F("fHistGammaPairInvMassOneAddMC", "Invariant mass of On-fly Photon Pairs, one additional Photon MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAdd = new TH1F("fHistGammaPairInvMassOnlyAdd", "Invariant mass of Photon Pairs, only additional Photons;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairInvMassOnlyAddMC= new TH1F("fHistGammaPairInvMassOnlyAddMC", "Invariant mass of Photon Pairs, only additional Photons MC;m_{inv} [GeV/#it{c}^{2}];Counts/(MeV/#it{c}^{2})", 1400, 0, 0.7);
    TH1F* fHistGammaPairDCA            = new TH1F("fHistGammaPairDCA","DCA between Photon Pairs;DCA [cm];Counts/(0.5mm)",200,0,10); 
    TH1F* fHistGammaPairDCAMC          = new TH1F("fHistGammaPairDCAMC","DCA between Photon Pairs MC;DCA [cm];Counts/(0.5mm)",200,0,10);

    //MC Sigma Topology
    TH1F* fHistPi0VertexvsMC           = new TH1F("fHistPi0VertexvsMC","#pi^{0} KF Decay Radius - MC Decay Radius;#Deltar [cm];Counts/(mm)", 400, -20, 20);
    TH1F* fHistMCSigmaPA               = new TH1F("fHistMCSigmaPA","#Sigma^{+}/#bar#Sigma^{-} PA;PA [rad];Counts/(0.005)",200,0,1);
    TH1F* fHistSigmaPA                 = new TH1F("fHistSigmaPA","#Sigma^{+}/#bar#Sigma^{-} PA;PA [rad];Counts/(0.005)",200,0,1);
    TH1F* fHistInvSigmaMass            = new TH1F("fHistInvSigmaMass","Invariant mass of #Sigma^{+}/#bar#Sigma^{-} Candidates;m_{inv} [GeV/#it{c}^{2}];Counts/(10 MeV/#it{c}^{2})", 300, 0, 3);
    TH1F* fHistMCOneGammaSigmaPA       = new TH1F("fHistMCOneGammaSigmaPA","#Sigma^{+} PA;PA [rad];Counts/(0.005)",200,0,1);
    TH1F* fHistPi0VertexMC             = new TH1F("fHistPi0VertexMC","#pi^{0} MC Decay Radius;r [cm];Counts/(mm)", 200, 0, 20);
    TH1F* fHistMCSigmaProtonkstar      = new TH1F("fHistMCSigmaProtonkstar","MC #Sigma^{+}/#bar#Sigma^{-}-p k*;k* [GeV/#it{c}];Counts/(5 MeV/#it{c})",1000,0,5);

    //Miscellaneous
    /*** EMPTY ***/

    //Alphanumeric Histogram Labels
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(1,"fMaxProtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(2,"fMinTPCClustProt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(3,"fMaxNsigProtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(4,"fMaxNsigProtTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(5,"fMaxpOnlyTPCPID");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(6,"fMinProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(7,"fMaxProtpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(8,"fMaxMCEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(9,"fMaxDaughtEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(10,"fMinTPCClustDaught");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(11,"fMaxNsigDaughtTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(12,"fMaxalpha");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(13,"fMaxqt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(14,"fMaxopenangle");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(15,"fMaxdeltatheta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(16,"fMinV0CPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(17,"fMinV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(18,"fMaxV0Radius");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(19,"fMaxphotonmass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(20,"fMinDCADaughtPV");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(21,"fMaxDCADaught");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(22,"fMaxElecEta");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(23,"fMinTPCClustElec");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(24,"fMaxNsigElecTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(25,"fMinNsigHadronTPC");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(26,"fMaxNsigElecTOF");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(27,"fMaxElecpt");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(28,"fMinPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(29,"fRejectGammaAutoCorr");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(30,"fMaxPi0Mass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(31,"fMaxSigmaPA");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(32,"fMaxSigmaMass");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(33,"fMinProtonDCAxy");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(34,"fMinProtonDCAz");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(35,"low kstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(36,"very low kstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(37,"very very low kstar");
    fHistCutBookKeeper->GetXaxis()->SetBinLabel(38,"Number of Fills");

    fHistCutBookKeeper->SetBinContent(1,fMaxProtEta);
    fHistCutBookKeeper->SetBinContent(2,fMinTPCClustProt);
    fHistCutBookKeeper->SetBinContent(3,fMaxNsigProtTPC);
    fHistCutBookKeeper->SetBinContent(4,fMaxNsigProtTOF);
    fHistCutBookKeeper->SetBinContent(5,fMaxpOnlyTPCPID);
    fHistCutBookKeeper->SetBinContent(6,fMinProtpt);
    fHistCutBookKeeper->SetBinContent(7,fMaxProtpt);
    fHistCutBookKeeper->SetBinContent(8,fMaxMCEta);
    fHistCutBookKeeper->SetBinContent(9,fMaxDaughtEta);
    fHistCutBookKeeper->SetBinContent(10,fMinTPCClustDaught);
    fHistCutBookKeeper->SetBinContent(11,fMaxNsigDaughtTPC);
    fHistCutBookKeeper->SetBinContent(12,fMaxalpha);
    fHistCutBookKeeper->SetBinContent(13,fMaxqt);
    fHistCutBookKeeper->SetBinContent(14,fMaxopenangle);
    fHistCutBookKeeper->SetBinContent(15,fMaxdeltatheta);
    fHistCutBookKeeper->SetBinContent(16,fMinV0CPA);
    fHistCutBookKeeper->SetBinContent(17,fMinV0Radius);
    fHistCutBookKeeper->SetBinContent(18,fMaxV0Radius);
    fHistCutBookKeeper->SetBinContent(19,fMaxphotonmass);
    fHistCutBookKeeper->SetBinContent(20,fMinDCADaughtPV);
    fHistCutBookKeeper->SetBinContent(21,fMaxDCADaught);
    fHistCutBookKeeper->SetBinContent(22,fMaxElecEta);
    fHistCutBookKeeper->SetBinContent(23,fMinTPCClustElec);
    fHistCutBookKeeper->SetBinContent(24,fMaxNsigElecTPC);
    fHistCutBookKeeper->SetBinContent(25,fMinNsigHadronTPC);
    fHistCutBookKeeper->SetBinContent(26,fMaxNsigElecTOF);
    fHistCutBookKeeper->SetBinContent(27,fMaxElecpt);
    fHistCutBookKeeper->SetBinContent(28,fRejectGammaAutoCorr);
    fHistCutBookKeeper->SetBinContent(29,fMinPi0Mass);
    fHistCutBookKeeper->SetBinContent(30,fMaxPi0Mass);
    fHistCutBookKeeper->SetBinContent(31,fMaxSigmaPA);
    fHistCutBookKeeper->SetBinContent(32,fMaxSigmaMass);
    fHistCutBookKeeper->SetBinContent(33,fMinProtonDCAxy);
    fHistCutBookKeeper->SetBinContent(34,fMinProtonDCAz);
    fHistCutBookKeeper->SetBinContent(35,flowkstar);
    fHistCutBookKeeper->SetBinContent(36,fverylowkstar);
    fHistCutBookKeeper->SetBinContent(37,fveryverylowkstar);
    fHistCutBookKeeper->SetBinContent(38,1);

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
    fHistMCCounter->GetXaxis()->SetBinLabel(13,"#gamma from pi^{0}");
    fHistMCCounter->GetXaxis()->SetBinLabel(14,"#gamma from #Sigma^{+}/#bar#Sigma^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(15,"#gamma from #Sigma^{+}/#bar#Sigma^{-} via pi^{0}");
    fHistMCCounter->GetXaxis()->SetBinLabel(16,"#gamma from pi^{0} conv.");
    fHistMCCounter->GetXaxis()->SetBinLabel(17,"#gamma from #Sigma^{+}/#bar#Sigma^{-} conv.");
    fHistMCCounter->GetXaxis()->SetBinLabel(18,"#gamma from #Sigma^{+}/#bar#Sigma^{-} via pi^{0} conv.");

    fHistEventCounter->GetXaxis()->SetBinLabel(1,"Events");

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

    fHistSigmaCounter->GetXaxis()->SetBinLabel(1, "Both from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(2, "One from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(3, "Additional #gammas");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(4, "4 Particle Reco");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(5, "#Sigma->p#gamma from Finder");
    fHistSigmaCounter->GetXaxis()->SetBinLabel(6, "#Sigma->p#gamma additional");

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
    fOutputList->Add(fHistVertexZ);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistEventCounter);
    fOutputList->Add(fHistTrackMultiplicity);
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
    fOutputList->Add(fHistSigmaCounter);
    //MC Kinematics              
    fOutputList->Add(fHistMCSigmaPt);
    fOutputList->Add(fHistMCAntiSigmaPt);
    fOutputList->Add(fHistMCPrimSigmaPt);
    fOutputList->Add(fHistMCPrimAntiSigmaPt);
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
    fOutputList->Add(fHistTrackpvsdEdx);
    fOutputList->Add(fHistTrackpvsbeta);
    fOutputList->Add(fHistTrackpt);
    //Proton QA
    fOutputList->Add(fHistProtonEtaPhiMC);
    fOutputList->Add(fHistProtonChi2MC);
    fOutputList->Add(fHistProtonTPCClusterMC);
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
    fOutputList->Add(fHistV0DaughtTPCClustMC);
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
    //MC Sigma Topology
    fOutputList->Add(fHistPi0VertexvsMC);     
    fOutputList->Add(fHistMCSigmaPA);
    fOutputList->Add(fHistSigmaPA);
    fOutputList->Add(fHistInvSigmaMass);         
    fOutputList->Add(fHistMCOneGammaSigmaPA);
    fOutputList->Add(fHistPi0VertexMC);
    fOutputList->Add(fHistMCSigmaProtonkstar);

/**************************************************************************/
    PostData(1, fOutputList);         
    PostData(2, fSigmaCandTree);         
    PostData(3, fSigmaCandTreeExtra);         

}//end of UserCreateOutputObjects()

//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlus::UserExec(Option_t *)
{  

  // Main loop. Called once for each event

  // Load the Input Event and check it
  aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aodEvent) {
    AliWarning("ERROR: AOD Event not available!");
    return;
  }

  mcEvent = MCEvent();           // Get MC event (called mcEvent) from the input file 
  if(mcEvent) isMonteCarlo = kTRUE; 
  else isMonteCarlo = kFALSE;

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
  
  Double_t vertexZ = aodEvent->GetPrimaryVertex()->GetZ();
  FillHistogram("fHistVertexZ",vertexZ);
  if(std::abs(vertexZ)>10) return;   //Return if vertex z position >10cm!                  
  
  primaryVtxPosX=primaryVtxPos[0];
  primaryVtxPosY=primaryVtxPos[1];
  primaryVtxPosZ=primaryVtxPos[2];

  //Magnetic Field
  Bz = aodEvent->GetMagneticField();    
  // Set Magnetic field for ALL KFParticles
  KFParticle::SetField(Bz);

  // Centrality
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(aodEvent->FindListObject("MultSelection"));   //Get Mult Selection..
  if(multSelection) Centrality = multSelection->GetMultiplicityPercentile("V0M");                                //..and retrieve centrality
  FillHistogram("fHistCentrality",Centrality);

  FillHistogram("fHistMCCounter",1); //Event Counter 
  FillHistogram("fHistEventCounter",1); //Event Counter 

  //Fill V0 maps
  fONMap.clear(); //clear the maps
  fOFFMap.clear();
  fV0ParticleIDArray.clear();

  Int_t nV0 = aodEvent->GetNumberOfV0s(); //Number of V0s in the event
  if(nV0 == 0) return;      //Return if there is no V0 to be processed

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
    if (trackN->GetSign() == trackP->GetSign()) continue;
    if (trackN->Charge() > 0){ //Check correct charge 
    trackN = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
    trackP = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
    }

    if(aodV0->GetOnFlyStatus()){
      //map convention: negative track first, positive track second
      fONMap.insert(make_pair(make_pair(trackN->GetID(), trackP->GetID()), iV0));
    }
    else if(!aodV0->GetOnFlyStatus()){
      //map convention: negative track first, positive track second
      fOFFMap.insert(make_pair(make_pair(trackN->GetID(), trackP->GetID()), iV0));
    }
  }//End of V0 Loop. Finished preparing the maps

/************************Start Event Processing**************************************/

  if(fDebug) cout << "Start of Event\n";

  //Process Protons
  if(fDebug) cout << "Processing Protons\n";
  if(fProcessProtons) FillProtonArray();

  //Process MC Particles
  if(fDebug&&isMonteCarlo) cout << "Processing MC Particles\n";
  if(isMonteCarlo && fProcessMCParticles) ProcessMCParticles();

  //Process V0s
  if(fDebug) cout << "Processing V0s\n";
  if(fProcessV0s) FillV0PhotonArray();

  //Find Additional Photons
  if(fDebug) cout << "Searching for additional Photons\n";
  if(fProcessAddPhoton) FindAddPhotons();

  //Process Electrons
  if(fDebug) cout << "Processing Electrons\n";
  if(fProcessElectrons) FillElectronArray(); 

  //Reconstruct Pi0 and Sigma+
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas\n";
  if(fProcessReco) ReconstructParticles();

  //Reconstruct Pi0 and Sigma+
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas. One offline Photon\n";
  if(fProcessRecoOff) ReconstructParticlesOff();

  //Reconstruct Pi0 and Sigma+
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas. Two offline Photons\n";
  if(fProcessRecoOff2) ReconstructParticlesOff2();

  //Reconstruct Sigma+
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas. Exotic decay channel\n";
  if(fProcessOneGamma) ReconstructParticlesOneGamma();

  if(fDebug) cout << "End of Event\n";

/************************End of Event Processing**************************************/

    PostData(1, fOutputList); // stream the analysis results of the current event to output manager        
    PostData(2, fSigmaCandTree);
    PostData(3, fSigmaCandTreeExtra);                  
    return;

}//end of UserExec()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillProtonArray() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fProtonArray.clear();
  Int_t countProton = 0;
  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  //Loop for Proton Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {

    //Initialisation of local variables
    Bool_t isReallyProton = kFALSE; //From MC
    Bool_t isProtonfromSigma = kFALSE; //From MC
    Bool_t isPrimarySigma = kFALSE; //From MC

    Bool_t   isTPCProton = kFALSE;
    Bool_t   isTOFProton = kFALSE;

    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack) {
      AliWarning("No AOD Track!");
      continue;
    }

    FillHistogram("fHistProtonStatistics",1);
    FillHistogram("fHistProtonStatisticsMC",1);
    FillHistogram("fHistProtonStatisticsSigmaMC",1);

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

    FillHistogram("fHistProtonStatistics",2);
    FillHistogram("fHistProtonStatisticsMC",2);
    FillHistogram("fHistProtonStatisticsSigmaMC",2);

    Double_t eta    = aodTrack->Eta();
    Double_t phi    = aodTrack->Phi();
    Double_t dEdx   = aodTrack->GetTPCsignal();
    Double_t p      = aodTrack->GetTPCmomentum();
    Double_t pt     = aodTrack->Pt();
    Int_t    nTPCClustProt = aodTrack->GetTPCNcls();
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

    //Check MC Truth
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        if(mcPart->GetPdgCode()==2212 || mcPart->GetPdgCode()==-2212){
          isReallyProton = kTRUE;
          if(mcPart->GetMother()!=-1){
            AliAODMCParticle* ProtonMotherPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mcPart->GetMother()));
            if(ProtonMotherPart){
              if(ProtonMotherPart->GetPdgCode()==3222 || ProtonMotherPart->GetPdgCode()==-3222){ 
                isProtonfromSigma = kTRUE;
                if(ProtonMotherPart->IsPrimary()||ProtonMotherPart->IsPhysicalPrimary()) isPrimarySigma = kTRUE;
              }//Is really Sigma
            }//MC Mother exists
          }//Proton has Mother
        }//is Proton or Anti-Proton
      }//MC Particle exists
    }//MC treatment

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

    // Cluster Cut
    if(nTPCClustProt<fMinTPCClustProt) continue;

    FillHistogram("fHistProtonStatistics",4);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",4);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",4);

    // TPC PID Cut
    if(!isTPCProton) continue;

    FillHistogram("fHistProtonStatistics",5);
    if(isReallyProton) FillHistogram("fHistProtonStatisticsMC",5);
    if(isProtonfromSigma) FillHistogram("fHistProtonStatisticsSigmaMC",5);

    // TOF PID Cut
    if(p > fMaxpOnlyTPCPID && !isTOFProton) continue;

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
    countProton++;

  } // End of proton track loop

return;

}//End of FillProtonArray  

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ProcessMCParticles() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  if(!isMonteCarlo) return;

  Int_t nMCTracks = mcEvent->GetNumberOfTracks();

  //loop over all MC tracks
  for(Int_t iMCtrack = 1; iMCtrack < nMCTracks; iMCtrack++){

    AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iMCtrack));
    if(!mcPart) {
      AliWarning("Could not retrieve AOD MC Particle!");
      continue;
    }

    if(TMath::Abs(mcPart->Eta())>fMaxMCEta) continue; //Acceptance Cut

    Int_t MCPartPDGCode = mcPart->PdgCode(); 

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

      AliAODMCParticle* Pi0Part = NULL;
      AliAODMCParticle* SigmaPart = NULL;
      Bool_t isgammafrompi0 = kFALSE;
      Bool_t isgammafromSigma = kFALSE;
      Bool_t isgammadirectlyfromSigma = kFALSE;

      if(mcPart->GetMother()!=-1){
        Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      }
      if(Pi0Part){
        if(TMath::Abs(Pi0Part->GetPdgCode())==3222){isgammadirectlyfromSigma = kTRUE; FillHistogram("fHistMCCounter",14);}
        if(Pi0Part->GetPdgCode()==111){isgammafrompi0 = kTRUE; FillHistogram("fHistMCCounter",13);}
        if(Pi0Part->GetPdgCode()==111&&Pi0Part->GetMother()!=-1) SigmaPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
      }
      if(SigmaPart){
        if(TMath::Abs(SigmaPart->GetPdgCode())==3222){
          FillHistogram("fHistMCSigmaPhotonPt",mcPart->Pt());
          FillHistogram("fHistMCCounter",15);
          isgammafromSigma = kTRUE;
        }
      }

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
      
      if(prdvtx.Perp()<180){ 
        FillHistogram("fHistMCCounter",12); FillHistogram("fHistMCConvPhotonPt",mcPart->Pt());
        if(isgammafromSigma) FillHistogram("fHistMCSigmaConvPhotonPt",mcPart->Pt());
        if(isgammafrompi0) FillHistogram("fHistMCCounter",16);
        if(isgammadirectlyfromSigma) FillHistogram("fHistMCCounter",17);
        if(isgammafromSigma) FillHistogram("fHistMCCounter",18);
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
      map<pair<int,int>, int>::iterator iter = fONMap.find(make_pair(trackN->GetID(),trackP->GetID()));      
      if(iter != fONMap.end()){
        Int_t lEquivalentONV0 = (*iter).second; // or iter->second;
        AliAODv0 *v0_on = aodEvent->GetV0(lEquivalentONV0);
        if(v0_on) continue;  //No further processing if On-fly V0 exists!
      }
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

    // Get kinematic values
    Double_t ptV0   = aodV0->Pt();
    Double_t pPos   = trackP->P();
    Double_t pNeg   = trackN->P();
    Double_t thetaPos = trackP->Theta();
    Double_t thetaNeg = trackN->Theta();
    Double_t alpha=aodV0->AlphaV0();
    Double_t qt=aodV0->PtArmV0();
    Double_t totangle=aodV0->OpenAngleV0();

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
            
      FillHistogram("fHistAddV0Statistics",0.5); //number of considered pairs
      FillHistogram("fHistAddV0StatisticsMC",0.5);
      FillHistogram("fHistAddV0StatisticsSigmaMC",0.5);
            
      //discard the ones obtained with on-the-fly finding
      map<pair<int,int>, int>::iterator iteron = fONMap.find(make_pair(ntrk->GetID(),ptrk->GetID()));
      map<pair<int,int>, int>::iterator iteroff = fOFFMap.find(make_pair(ntrk->GetID(),ptrk->GetID()));
      
      if(iteron != fONMap.end()){
        Int_t lEquivalentONV0 = (*iteron).second; // or iter->second;
        AliAODv0 *v0_on = aodEvent->GetV0(lEquivalentONV0);
        if(v0_on) continue;  //No further processing if On-fly V0 exists!
      }

      if(iteroff != fOFFMap.end()){
        Int_t lEquivalentOFFV0 = (*iteroff).second; // or iter->second;
        AliAODv0 *v0_off = aodEvent->GetV0(lEquivalentOFFV0);
        if(v0_off) continue; //No further processing if Offline V0 exists!
      }

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
      FillHistogram("fHistAddV0Statistics",1.5); //number of pairs not found by the Finders
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",1.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",1.5);

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
            
      FillHistogram("fHistAddV0Statistics",2.5); //passed dca cut
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",2.5);      
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",2.5);      

      if ((xn+xp) > 2*fMaxV0Radius) continue;
      if ((xn+xp) < 2*fMinV0Radius) continue;
            
      FillHistogram("fHistAddV0Statistics",3.5); //passed X within R2D cut
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",3.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",3.5);

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
            
      FillHistogram("fHistAddV0Statistics",4.5); //passed radius cut
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",4.5);      
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",4.5);      

      Float_t cpa=vertex.GetV0CosineOfPointingAngle(primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ);

      FillHistogram("fHistPairCPA",cpa);
      if(isReallyPhoton) FillHistogram("fHistPairCPAMC",cpa);
      if(isPhotonfromSigma) FillHistogram("fHistPairCPAMCSigma",cpa);  

      //Simple cosine cut (no pt dependence for now)
      if (cpa < fMinV0CPA) continue;
            
      FillHistogram("fHistAddV0Statistics",5.5); //passed cosPA
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",5.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",5.5);

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
      FillHistogram("fHistAddV0Statistics",6.5); //passed alpha
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",6.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",6.5);

      if(TMath::Abs(qt) > fMaxqt) continue;
      FillHistogram("fHistAddV0Statistics",7.5); //passed qt
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",7.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",7.5);

      // Angle Cut
      if(TMath::Abs(totangle) > fMaxopenangle) continue;
      FillHistogram("fHistAddV0Statistics",8.5); //passed openangle
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",8.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",8.5);

      if(TMath::Abs(deltatheta) > fMaxdeltatheta) continue;
      FillHistogram("fHistAddV0Statistics",9.5); //passed deltatheta
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",9.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",9.5);

      // Inv. Mass Cut
      if(photonmass > fMaxphotonmass) continue;
      FillHistogram("fHistAddV0Statistics",10.5); //passed inv mass 
      if(isReallyPhoton) FillHistogram("fHistAddV0StatisticsMC",10.5);
      if(isPhotonfromSigma) FillHistogram("fHistAddV0StatisticsSigmaMC",10.5);

      if(isReallyPhoton) FillHistogram("fHistPairptwCutsMC",photon.Perp());
      FillHistogram("fHistPairptwCuts",photon.Perp());

      //Store ESD IDs of used electrons
      fV0ParticleIDArray.push_back(ntrk->GetID());
      fV0ParticleIDArray.push_back(ptrk->GetID());

      //Store the AOD! Track indices as pairs
      PairIndexArray.push_back(make_pair(nidx,pidx));  

    }//End of Positive Particle Track Loop
  }//End of Negative Particle Track Loop
  
  return;

}//End of FindAddPhotons() const

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillElectronArray() {
  
  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fElectronArray.clear();
  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  //Loop for Proton Selection
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
  Int_t countPi0 = 0;

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
      Bool_t isReallyPi0 = kFALSE; 
      Bool_t isReallyPi0fromSigma = kFALSE;
      Bool_t isReallyPi0fromDelta = kFALSE;
      Bool_t isPrimary = kFALSE;
      Int_t Pi0MotherLabel=-1;

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
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22&&TMath::Abs(V0Part1->GetLabel())!=TMath::Abs(V0Part2->GetLabel())){
                if(V0Part1->GetMother()!=-1&&V0Part1->GetMother()==V0Part2->GetMother()){

                  AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                  if(Pi0Part){if(Pi0Part->GetPdgCode()==111){
                  
                    isReallyPi0=kTRUE;

                    TVector3 prdvtx1(V0Part1->Xv(),V0Part1->Yv(),V0Part1->Zv());   
                    TVector3 prdvtx2(V0Part2->Xv(),V0Part2->Yv(),V0Part2->Zv());   
                    TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
                    MCPi0DCAPV = prdvtx.Perp();

                    AliAODMCParticle* Pi0Mother = NULL;
                    if(Pi0Part->GetMother()!=-1){
                      Pi0MotherLabel=Pi0Part->GetMother(); 
                      Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                    }
                    if(Pi0Mother){
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) {isReallyPi0fromSigma = kTRUE; if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE;}
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

      if(pi0mass<fRejectGammaAutoCorr&&(track1->GetID()==track3->GetID()||track2->GetID()==track4->GetID())) continue; //Clean up auto-correlations!

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

      if(isReallyPi0){
  	    FillHistogram("fHistGammaPairInvMassMC",pi0mass);
  	    FillHistogram("fHistGammaPairDCAMC",PhotPhotDCA);
  	    if(v0_1->GetOnFlyStatus()&&v0_2->GetOnFlyStatus()) FillHistogram("fHistGammaPairInvMassMCOnfly",pi0mass);
      }
      
      if(pi0mass>0.12 && pi0mass<0.15) countPi0++; //Just for counting purposes here

      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

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
              if(ProtonPart->GetMother()==Pi0MotherLabel){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",1);}
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {
            if(q==k) continue;
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;
            pairs++;
            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 kstar=sigmamomentum-protonmomentum;
            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar.Mag());                        
            if(kstar.Mag()<flowkstar) pairslowkstar++;
            if(kstar.Mag()<fverylowkstar) pairsverylowkstar++;
            if(kstar.Mag()<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          fIsV01fromFinder = kTRUE;
          fIsV02fromFinder = kTRUE;
          fIsV01Onthefly = v0_1->GetOnFlyStatus();
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = (ULong_t)aodEvent->GetTriggerMask();
          fSigCentrality = Centrality;
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
          fPhoton1Radius = TMath::Sqrt(v0_1->DecayVertexV0X()*v0_1->DecayVertexV0X()+v0_1->DecayVertexV0Y()*v0_1->DecayVertexV0Y());
          fPhoton2Radius = TMath::Sqrt(v0_2->DecayVertexV0X()*v0_2->DecayVertexV0X()+v0_2->DecayVertexV0Y()*v0_2->DecayVertexV0Y()); 
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2;
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
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
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
  Int_t countPi0 = 0;

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
      Bool_t isReallyPi0 = kFALSE; 
      Bool_t isReallyPi0fromSigma = kFALSE;
      Bool_t isReallyPi0fromDelta = kFALSE;
      Bool_t isPrimary = kFALSE;
      Int_t Pi0MotherLabel=-1;

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
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22&&TMath::Abs(V0Part1->GetLabel())!=TMath::Abs(V0Part2->GetLabel())){
                if(V0Part1->GetMother()!=-1&&V0Part1->GetMother()==V0Part2->GetMother()){

                  AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                  if(Pi0Part){if(Pi0Part->GetPdgCode()==111){
                  
                    isReallyPi0=kTRUE;

                    TVector3 prdvtx1(V0Part1->Xv(),V0Part1->Yv(),V0Part1->Zv());   
                    TVector3 prdvtx2(V0Part2->Xv(),V0Part2->Yv(),V0Part2->Zv());   
                    TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
                    MCPi0DCAPV = prdvtx.Perp();

                    AliAODMCParticle* Pi0Mother = NULL;
                    if(Pi0Part->GetMother()!=-1){
                      Pi0MotherLabel=Pi0Part->GetMother(); 
                      Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                    }
                    if(Pi0Mother){
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) {isReallyPi0fromSigma = kTRUE; if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE;}
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

      if(pi0mass<fRejectGammaAutoCorr&&(track1->GetID()==track3->GetID()||track2->GetID()==track4->GetID())) continue; //Clean up auto-correlations!

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

  	  FillHistogram("fHistGammaPairInvMassOneAdd",pi0mass);
      if(isReallyPi0){FillHistogram("fHistGammaPairInvMassOneAddMC",pi0mass);}

      if(pi0mass>0.12 && pi0mass<0.15) countPi0++; //Just for counting purposes here

      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

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
              if(ProtonPart->GetMother()==Pi0MotherLabel){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",2);}
                else if (isReallyPi0fromDelta) FillHistogram("fHistMCdeltamass",sigmaplusmass); 
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {
            if(q==k) continue;
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;
            pairs++;
            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 kstar=sigmamomentum-protonmomentum;
            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar.Mag());
            if(kstar.Mag()<flowkstar) pairslowkstar++;
            if(kstar.Mag()<fverylowkstar) pairsverylowkstar++;
            if(kstar.Mag()<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          fIsV01fromFinder = kFALSE;
          fIsV02fromFinder = kTRUE;
          fIsV01Onthefly = kFALSE;
          fIsV02Onthefly = v0_2->GetOnFlyStatus();
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = (ULong_t)aodEvent->GetTriggerMask();
          fSigCentrality = Centrality;
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
          fPhoton1Radius = TMath::Sqrt(v0_1.Xv()*v0_1.Xv()+v0_1.Yv()*v0_1.Yv());
          fPhoton2Radius = TMath::Sqrt(v0_2->DecayVertexV0X()*v0_2->DecayVertexV0X()+v0_2->DecayVertexV0Y()*v0_2->DecayVertexV0Y());
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2; 
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
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
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
  Int_t countPi0 = 0;

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
      Bool_t isReallyPi0 = kFALSE; 
      Bool_t isReallyPi0fromSigma = kFALSE;
      Bool_t isReallyPi0fromDelta = kFALSE;
      Bool_t isPrimary = kFALSE;
      Int_t Pi0MotherLabel=-1;

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
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22&&TMath::Abs(V0Part1->GetLabel())!=TMath::Abs(V0Part2->GetLabel())){
                if(V0Part1->GetMother()!=-1&&V0Part1->GetMother()==V0Part2->GetMother()){

                  AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part1->GetMother())));
                  if(Pi0Part){if(Pi0Part->GetPdgCode()==111){
                  
                    isReallyPi0=kTRUE;

                    TVector3 prdvtx1(V0Part1->Xv(),V0Part1->Yv(),V0Part1->Zv());   
                    TVector3 prdvtx2(V0Part2->Xv(),V0Part2->Yv(),V0Part2->Zv());   
                    TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
                    MCPi0DCAPV = prdvtx.Perp();

                    AliAODMCParticle* Pi0Mother = NULL;
                    if(Pi0Part->GetMother()!=-1){
                      Pi0MotherLabel=Pi0Part->GetMother(); 
                      Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                    }
                    if(Pi0Mother){
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) {isReallyPi0fromSigma = kTRUE; if(Pi0Mother->IsPrimary()||Pi0Mother->IsPhysicalPrimary()) isPrimary = kTRUE;}
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

      if(pi0mass<fRejectGammaAutoCorr&&(track1->GetID()==track3->GetID()||track2->GetID()==track4->GetID())) continue; //Clean up auto-correlations!

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

  	  FillHistogram("fHistGammaPairInvMassOnlyAdd",pi0mass);
      if(isReallyPi0){FillHistogram("fHistGammaPairInvMassOnlyAddMC",pi0mass);}

      if(pi0mass>0.12 && pi0mass<0.15) countPi0++; //Just for counting purposes here

      if(pi0mass<fMinPi0Mass || pi0mass>fMaxPi0Mass) continue;  //Coarse mass cut to reduce combinatorics!

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
              if(ProtonPart->GetMother()==Pi0MotherLabel){
                if(isReallyPi0fromSigma) {isReallySigma = kTRUE; FillHistogram("fHistSigmaCounter",3);}
                else if (isReallyPi0fromDelta) FillHistogram("fHistMCdeltamass",sigmaplusmass); 
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          if(TMath::Abs(DCAxy)<fMinProtonDCAxy) continue; //Coarse cut on DCA to PV to reduce the tree size
          if(TMath::Abs(DCAz)<fMinProtonDCAz) continue;

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
          Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
          if(isReallySigma) FillHistogram("fHistMCSigmaPA",SigmaPointingAngle);
          FillHistogram("fHistSigmaPA",SigmaPointingAngle);
          if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

          Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
          for(Int_t q=0; q<nProton; q++) {
            if(q==k) continue;
            AliAODTrack *pairprot;
            pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
            if(!pairprot) continue;
            if(pairprot->Charge()!=prot->Charge()) continue;
            pairs++;
            TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
            TVector3 kstar=sigmamomentum-protonmomentum;
            if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar.Mag());
            if(kstar.Mag()<flowkstar) pairslowkstar++;
            if(kstar.Mag()<fverylowkstar) pairsverylowkstar++;
            if(kstar.Mag()<fveryverylowkstar) pairsveryverylowkstar++;
          }        

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; 
          if(isReallySigma) fIsMCSigma = kTRUE;
          fIsMCPrimary = kFALSE;
          if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
          fIsV01fromFinder = kFALSE;
          fIsV02fromFinder = kFALSE;
          fIsV01Onthefly = kFALSE;
          fIsV02Onthefly = kFALSE;
          fSigRunnumber = aodEvent->GetRunNumber();
          fSigTriggerMask = (ULong_t)aodEvent->GetTriggerMask();
          fSigCentrality = Centrality;
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
          fPhoton1Radius = TMath::Sqrt(v0_1.Xv()*v0_1.Xv()+v0_1.Yv()*v0_1.Yv());
          fPhoton2Radius = TMath::Sqrt(v0_2.Xv()*v0_2.Xv()+v0_2.Yv()*v0_2.Yv());
          fPhoton1DCAPV = DCAPV1;
          fPhoton2DCAPV = DCAPV2; 
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
          fProtonDCAtoPVxy = DCAxy; 
          fProtonDCAtoPVz = DCAz; 
          fProtonPi0DCA = ProtPi0DCA;
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
          if(SigmaPart->IsPrimary()||SigmaPart->IsPhysicalPrimary()) isPrimary = kTRUE;
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

      TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
      TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);
      Float_t SigmaPointingAngle = sigmamomentum.Angle(sigmavertex);
      if(isReallySigma) FillHistogram("fHistMCOneGammaSigmaPA",SigmaPointingAngle);
      if(SigmaPointingAngle>fMaxSigmaPA) continue;  //Coarse cut on the Pointing Angle to reduce the tree size

      Short_t pairs = 0, pairslowkstar = 0, pairsverylowkstar = 0, pairsveryverylowkstar = 0;
      for(Int_t q=0; q<nProton; q++) {
        if(q==j) continue;
        AliAODTrack *pairprot;
        pairprot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray.at(q));
        if(!pairprot) continue;
        if(pairprot->Charge()!=prot->Charge()) continue;
        pairs++;
        TVector3 protonmomentum(pairprot->Px(),pairprot->Py(),pairprot->Pz());
        TVector3 kstar=sigmamomentum-protonmomentum;
        if(isReallySigma) FillHistogram("fHistMCSigmaProtonkstar",kstar.Mag());
        if(kstar.Mag()<flowkstar) pairslowkstar++;
        if(kstar.Mag()<fverylowkstar) pairsverylowkstar++;
        if(kstar.Mag()<fveryverylowkstar) pairsveryverylowkstar++;
      }        

      // Fill the Extra Sigma Candidate Trees
      fIsMCSigma = kFALSE; 
      if(isReallySigma) fIsMCSigma = kTRUE;
      fIsMCPrimary = kFALSE;
      if(isReallySigma&&isPrimary) fIsMCPrimary = kTRUE;
      fIsV0fromFinder = kTRUE;
      fIsV0Onthefly = v0->GetOnFlyStatus();
      fSigRunnumber = aodEvent->GetRunNumber();
      fSigTriggerMask = (ULong_t)aodEvent->GetTriggerMask();
      fSigCentrality = Centrality;
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
      fPhotonRadius = TMath::Sqrt(v0->DecayVertexV0X()*v0->DecayVertexV0X()+v0->DecayVertexV0Y()*v0->DecayVertexV0Y()); 
      fPhotonDCAPV = DCAPV;
      fPhotonPx = trackPhoton.Px(); 
      fPhotonPy = trackPhoton.Py(); 
      fPhotonPz = trackPhoton.Pz(); 
      fExtPhotProtDCA = ProtPhotDCA;
      fProtonPx = prot->Px(); 
      fProtonPy = prot->Py(); 
      fProtonPz = prot->Pz(); 
      fProtonDCAtoPVxy = DCAxy; 
      fProtonDCAtoPVz = DCAz; 
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
