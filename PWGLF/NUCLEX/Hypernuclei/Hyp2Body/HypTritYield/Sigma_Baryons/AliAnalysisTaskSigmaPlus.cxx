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

class TTree;
class TParticle;
class TVector3;
class AliAODVertex;
class AliAODv0;

class AliAnalysisTaskSigmaPlus;    // analysis class

using std::cout;            

ClassImp(AliAnalysisTaskSigmaPlus) // classimp: necessary for root

AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus() : AliAnalysisTaskSE(), 
fOutputList(0), aodEvent(0x0), mcEvent(0x0),        
AODMCTrackArray(0x0), fPIDResponse(0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), f4PartSigmaCandTree(0x0), f3PartPi0CandTree(0x0), 
fProtonArray(0x0), fAntiProtonArray(0x0), fConvPhotonArray(0x0),
fElectronArray(0x0), fPositronArray(0x0), fElectronPairArray(0x0), fPositronPairArray(0x0),
cElectronMass(0), cProtonMass(0), cPionMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(0), primaryVtxPosY(0), primaryVtxPosZ(0), nTracks(0),

fMaxProtEta(0.9),    
fMinTPCClustProt(70),
fMaxNsigProtTPC(3),
fMaxNsigProtTOF(3),  
fMaxpOnlyTPCPID(0.75),
fMinProtpt(0.5), 
fMaxProtpt(10),   

fMaxElecEta(0.9),
fMinTPCClustElec(60),
fMaxNsigElecTPC(3),
fMaxNsigElecTOF(3),
fNsigHadronreject(2),

fMaxMCEta(0.9),

fMaxDaughtEta(0.9),    
fMinTPCClustDaught(60),
fMaxNsigDaughtTPC(3),
fMaxalpha(0.9),
fMaxqt(0.04),
fMaxopenangle(0.2),
fMaxdeltatheta(0.1),
fMaxphotonmass(0.05),

fMaxPairalpha(0.9),
fMaxPairqt(0.04),
fMaxPairopenangle(0.2),
fMaxPairDCA(2),
fMaxPairphotonmass(0.06),

fMinSingleElectronPt(0.25),
fMinPi0mass3Part(0.12),
fMaxPi0mass3Part(0.15),
fMaxElecPt3Part(0.9),
fMaxSigmaPA3Part(0.175),
fMinPi0DCAtoPV3Part(1),
fMaxProtPi0DCA3Part(2),
fMaxPhotElecDCA(2),

fMinPi0mass(0.12),
fMaxPi0mass(0.15),
fMaxElecPt(0.9),
fMaxSigmaPA(0.175),
fMinPi0DCAtoPV(1),
fMaxProtPi0DCA(2),
fMaxPhotPhotDCA(2),
fMindcazProt(0),
fMindcaxyProt(0), 

fIsMCSigma(kFALSE),
fInvSigMass(-999),
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
fInvPi0Mass(-999),
fPi0Px(-999),
fPi0Py(-999),
fPi0Pz(-999),
fPi0DaughtPt1(-999),
fPi0DaughtPt2(-999),
fPi0DaughtPt3(-999),
fPi0DaughtPt4(-999),
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

f4PartIsMCSigma(kFALSE),
f4PartInvSigMass(-999),
f4PartSigCharge(-999),
f4PartSigPx(-999),
f4PartSigPy(-999),
f4PartSigPz(-999),
f4PartPrimVertX(-999),
f4PartPrimVertY(-999),
f4PartPrimVertZ(-999),
f4PartSigDecayVertX(-999),
f4PartSigDecayVertY(-999),
f4PartSigDecayVertZ(-999),
f4PartInvPi0Mass(-999),
f4PartPi0Px(-999),
f4PartPi0Py(-999),
f4PartPi0Pz(-999),
f4PartPi0DaughtPt1(-999),
f4PartPi0DaughtPt2(-999),
f4PartPi0DaughtPt3(-999),
f4PartPi0DecayVertX(-999),
f4PartPi0DecayVertY(-999),
f4PartPi0DecayVertZ(-999),
f4PartPi0PhotPhotDCA(-999),
f4PartProtonPx(-999),
f4PartProtonPy(-999),
f4PartProtonPz(-999),
f4PartProtonDCAtoPVxy(-999),
f4PartProtonDCAtoPVz(-999),
f4PartProtonPi0DCA(-999),

fInv3PartPi0Mass(-999),
f3PartPi0Px(-999),
f3PartPi0Py(-999),
f3PartPi0Pz(-999),
f3PartPi0DaughtPt1(-999),
f3PartPi0DaughtPt2(-999),
f3PartPi0DaughtPt3(-999),
f3PartPi0DecayVertX(-999),
f3PartPi0DecayVertY(-999),
f3PartPi0DecayVertZ(-999),
f3PartPi0PhotPhotDCA(-999),
f3PartPi0PhotConvx(-999),
f3PartPi0PhotConvy(-999),
f3PartPi0PhotConvz(-999),
f3PartPrimVertX(-999),
f3PartPrimVertY(-999),
f3PartPrimVertZ(-999),
f3PartPi0ElecPosx(-999),
f3PartPi0ElecPosy(-999),
f3PartPi0ElecPosz(-999)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::AliAnalysisTaskSigmaPlus(const char* name) : AliAnalysisTaskSE(name),
fOutputList(0), aodEvent(0x0), mcEvent(0x0),        
AODMCTrackArray(0x0), fPIDResponse(0), isMonteCarlo(kFALSE),
fSigmaCandTree(0x0), f4PartSigmaCandTree(0x0), f3PartPi0CandTree(0x0), 
fProtonArray(0x0), fAntiProtonArray(0x0), fConvPhotonArray(0x0),
fElectronArray(0x0), fPositronArray(0x0), fElectronPairArray(0x0), fPositronPairArray(0x0),
cElectronMass(0), cProtonMass(0), cPionMass(0), cPi0Mass(0), c(0), Bz(0),
primaryVtxPosX(0), primaryVtxPosY(0), primaryVtxPosZ(0), nTracks(0),

fMaxProtEta(0.9),    
fMinTPCClustProt(70),
fMaxNsigProtTPC(3),
fMaxNsigProtTOF(3),  
fMaxpOnlyTPCPID(0.75),
fMinProtpt(0.5), 
fMaxProtpt(10),   

fMaxElecEta(0.9),
fMinTPCClustElec(60),
fMaxNsigElecTPC(3),
fMaxNsigElecTOF(3),
fNsigHadronreject(2),

fMaxMCEta(0.9),

fMaxDaughtEta(0.9),    
fMinTPCClustDaught(60),
fMaxNsigDaughtTPC(3),
fMaxalpha(0.9),
fMaxqt(0.04),
fMaxopenangle(0.2),
fMaxdeltatheta(0.1),
fMaxphotonmass(0.05),

fMaxPairalpha(0.9),
fMaxPairqt(0.04),
fMaxPairopenangle(0.2),
fMaxPairDCA(2),
fMaxPairphotonmass(0.06),

fMinSingleElectronPt(0.25),
fMinPi0mass3Part(0.12),
fMaxPi0mass3Part(0.15),
fMaxElecPt3Part(0.9),
fMaxSigmaPA3Part(0.175),
fMinPi0DCAtoPV3Part(1),
fMaxProtPi0DCA3Part(2),
fMaxPhotElecDCA(2),

fMinPi0mass(0.12),
fMaxPi0mass(0.15),
fMaxElecPt(0.9),
fMaxSigmaPA(0.175),
fMinPi0DCAtoPV(1),
fMaxProtPi0DCA(2),
fMaxPhotPhotDCA(2),
fMindcazProt(0),
fMindcaxyProt(0), 

fIsMCSigma(kFALSE),
fInvSigMass(-999),
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
fInvPi0Mass(-999),
fPi0Px(-999),
fPi0Py(-999),
fPi0Pz(-999),
fPi0DaughtPt1(-999),
fPi0DaughtPt2(-999),
fPi0DaughtPt3(-999),
fPi0DaughtPt4(-999),
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

f4PartIsMCSigma(kFALSE),
f4PartInvSigMass(-999),
f4PartSigCharge(-999),
f4PartSigPx(-999),
f4PartSigPy(-999),
f4PartSigPz(-999),
f4PartPrimVertX(-999),
f4PartPrimVertY(-999),
f4PartPrimVertZ(-999),
f4PartSigDecayVertX(-999),
f4PartSigDecayVertY(-999),
f4PartSigDecayVertZ(-999),
f4PartInvPi0Mass(-999),
f4PartPi0Px(-999),
f4PartPi0Py(-999),
f4PartPi0Pz(-999),
f4PartPi0DaughtPt1(-999),
f4PartPi0DaughtPt2(-999),
f4PartPi0DaughtPt3(-999),
f4PartPi0DecayVertX(-999),
f4PartPi0DecayVertY(-999),
f4PartPi0DecayVertZ(-999),
f4PartPi0PhotPhotDCA(-999),
f4PartProtonPx(-999),
f4PartProtonPy(-999),
f4PartProtonPz(-999),
f4PartProtonDCAtoPVxy(-999),
f4PartProtonDCAtoPVz(-999),
f4PartProtonPi0DCA(-999),

fInv3PartPi0Mass(-999),
f3PartPi0Px(-999),
f3PartPi0Py(-999),
f3PartPi0Pz(-999),
f3PartPi0DaughtPt1(-999),
f3PartPi0DaughtPt2(-999),
f3PartPi0DaughtPt3(-999),
f3PartPi0DecayVertX(-999),
f3PartPi0DecayVertY(-999),
f3PartPi0DecayVertZ(-999),
f3PartPi0PhotPhotDCA(-999),
f3PartPi0PhotConvx(-999),
f3PartPi0PhotConvy(-999),
f3PartPi0PhotConvz(-999),
f3PartPrimVertX(-999),
f3PartPrimVertY(-999),
f3PartPrimVertZ(-999),
f3PartPi0ElecPosx(-999),
f3PartPi0ElecPosy(-999),
f3PartPi0ElecPosz(-999)
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
    DefineOutput(3, TTree::Class());    //Additional Output: 4Part Sigma Candidate Tree
    DefineOutput(4, TTree::Class());    //Additional Output: 3Part Pi0 Candidate Tree
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlus::~AliAnalysisTaskSigmaPlus()
{
    // destructor
    // delete objects from memory at the end of the task
    if(fOutputList) delete fOutputList;     

    if(fSigmaCandTree) delete fSigmaCandTree;
    if(f4PartSigmaCandTree) delete f4PartSigmaCandTree;
    if(f3PartPi0CandTree) delete f3PartPi0CandTree;

    if(fProtonArray)       delete fProtonArray;      
    if(fAntiProtonArray)   delete fAntiProtonArray;  
    if(fConvPhotonArray)   delete fConvPhotonArray;  
    if(fElectronArray)     delete fElectronArray;    
    if(fPositronArray)     delete fPositronArray;    
    if(fElectronPairArray) delete fElectronPairArray;
    if(fPositronPairArray) delete fPositronPairArray;
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
    fSigmaCandTree->Branch("fInvSigMass",&fInvSigMass,"fInvSigMass/F");
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
    fSigmaCandTree->Branch("fInvPi0Mass",&fInvPi0Mass,"fInvPi0Mass/F");
    fSigmaCandTree->Branch("fPi0Px",&fPi0Px,"fPi0Px/F");
    fSigmaCandTree->Branch("fPi0Py",&fPi0Py,"fPi0Py/F");
    fSigmaCandTree->Branch("fPi0Pz",&fPi0Pz,"fPi0Pz/F");
    fSigmaCandTree->Branch("fPi0DaughtPt1",&fPi0DaughtPt1,"fPi0DaughtPt1/F");
    fSigmaCandTree->Branch("fPi0DaughtPt2",&fPi0DaughtPt2,"fPi0DaughtPt2/F");
    fSigmaCandTree->Branch("fPi0DaughtPt3",&fPi0DaughtPt3,"fPi0DaughtPt3/F");
    fSigmaCandTree->Branch("fPi0DaughtPt4",&fPi0DaughtPt4,"fPi0DaughtPt4/F");
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

    f4PartSigmaCandTree = new TTree("f4PartSigmaCandTree","Tree of Sigma Candidates (3 Electron Method)");
    f4PartSigmaCandTree->Branch("f4PartIsMCSigma",&f4PartIsMCSigma,"f4PartIsMCSigma/O");
    f4PartSigmaCandTree->Branch("f4PartInvSigMass",&f4PartInvSigMass,"f4PartInvSigMass/F");
    f4PartSigmaCandTree->Branch("f4PartSigCharge",&f4PartSigCharge,"f4PartSigCharge/F");
    f4PartSigmaCandTree->Branch("f4PartSigPx",&f4PartSigPx,"f4PartSigPx/F");
    f4PartSigmaCandTree->Branch("f4PartSigPy",&f4PartSigPy,"f4PartSigPy/F");
    f4PartSigmaCandTree->Branch("f4PartSigPz",&f4PartSigPz,"f4PartSigPz/F");
    f4PartSigmaCandTree->Branch("f4PartPrimVertX",&f4PartPrimVertX,"f4PartPrimVertX/F");
    f4PartSigmaCandTree->Branch("f4PartPrimVertY",&f4PartPrimVertY,"f4PartPrimVertY/F");
    f4PartSigmaCandTree->Branch("f4PartPrimVertZ",&f4PartPrimVertZ,"f4PartPrimVertZ/F");
    f4PartSigmaCandTree->Branch("f4PartSigDecayVertX",&f4PartSigDecayVertX,"f4PartSigDecayVertX/F");
    f4PartSigmaCandTree->Branch("f4PartSigDecayVertY",&f4PartSigDecayVertY,"f4PartSigDecayVertY/F");
    f4PartSigmaCandTree->Branch("f4PartSigDecayVertZ",&f4PartSigDecayVertZ,"f4PartSigDecayVertZ/F");
    f4PartSigmaCandTree->Branch("f4PartInvPi0Mass",&f4PartInvPi0Mass,"f4PartInvPi0Mass/F");
    f4PartSigmaCandTree->Branch("f4PartPi0Px",&f4PartPi0Px,"f4PartPi0Px/F");
    f4PartSigmaCandTree->Branch("f4PartPi0Py",&f4PartPi0Py,"f4PartPi0Py/F");
    f4PartSigmaCandTree->Branch("f4PartPi0Pz",&f4PartPi0Pz,"f4PartPi0Pz/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DaughtPt1",&f4PartPi0DaughtPt1,"f4PartPi0DaughtPt1/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DaughtPt2",&f4PartPi0DaughtPt2,"f4PartPi0DaughtPt2/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DaughtPt3",&f4PartPi0DaughtPt3,"f4PartPi0DaughtPt3/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DecayVertX",&f4PartPi0DecayVertX,"f4PartPi0DecayVertX/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DecayVertY",&f4PartPi0DecayVertY,"f4PartPi0DecayVertY/F");
    f4PartSigmaCandTree->Branch("f4PartPi0DecayVertZ",&f4PartPi0DecayVertZ,"f4PartPi0DecayVertZ/F");
    f4PartSigmaCandTree->Branch("f4PartPi0PhotPhotDCA",&f4PartPi0PhotPhotDCA,"f4PartPi0PhotPhotDCA/F");
    f4PartSigmaCandTree->Branch("f4PartProtonPx",&f4PartProtonPx,"f4PartProtonPx/F");
    f4PartSigmaCandTree->Branch("f4PartProtonPy",&f4PartProtonPy,"f4PartProtonPy/F");
    f4PartSigmaCandTree->Branch("f4PartProtonPz",&f4PartProtonPz,"f4PartProtonPz/F");
    f4PartSigmaCandTree->Branch("f4PartProtonDCAtoPVxy",&f4PartProtonDCAtoPVxy,"f4PartProtonDCAtoPVxy/F");
    f4PartSigmaCandTree->Branch("f4PartProtonDCAtoPVz",&f4PartProtonDCAtoPVz,"f4PartProtonDCAtoPVz/F");
    f4PartSigmaCandTree->Branch("f4PartProtonPi0DCA",&f4PartProtonPi0DCA,"f4PartProtonPi0DCA/F");

    f3PartPi0CandTree = new TTree("f3PartPi0CandTree","Tree of Pi0 Candidates (3 Electron Method)");
    f3PartPi0CandTree->Branch("fInv3PartPi0Mass",&fInv3PartPi0Mass,"fInv3PartPi0Mass/F");
    f3PartPi0CandTree->Branch("f3PartPi0Px",&f3PartPi0Px,"f3PartPi0Px/F");
    f3PartPi0CandTree->Branch("f3PartPi0Py",&f3PartPi0Py,"f3PartPi0Py/F");
    f3PartPi0CandTree->Branch("f3PartPi0Pz",&f3PartPi0Pz,"f3PartPi0Pz/F");
    f3PartPi0CandTree->Branch("f3PartPi0DaughtPt1",&f3PartPi0DaughtPt1,"f3PartPi0DaughtPt1/F");
    f3PartPi0CandTree->Branch("f3PartPi0DaughtPt2",&f3PartPi0DaughtPt2,"f3PartPi0DaughtPt2/F");
    f3PartPi0CandTree->Branch("f3PartPi0DaughtPt3",&f3PartPi0DaughtPt3,"f3PartPi0DaughtPt3/F");
    f3PartPi0CandTree->Branch("f3PartPi0DecayVertX",&f3PartPi0DecayVertX,"f3PartPi0DecayVertX/F");
    f3PartPi0CandTree->Branch("f3PartPi0DecayVertY",&f3PartPi0DecayVertY,"f3PartPi0DecayVertY/F");
    f3PartPi0CandTree->Branch("f3PartPi0DecayVertZ",&f3PartPi0DecayVertZ,"f3PartPi0DecayVertZ/F");
    f3PartPi0CandTree->Branch("f3PartPi0PhotPhotDCA",&f3PartPi0PhotPhotDCA,"f3PartPi0PhotPhotDCA/F");
    f3PartPi0CandTree->Branch("f3PartPi0PhotConvx",&f3PartPi0PhotConvx,"f3PartPi0PhotConvx/F");
    f3PartPi0CandTree->Branch("f3PartPi0PhotConvy",&f3PartPi0PhotConvy,"f3PartPi0PhotConvy/F");
    f3PartPi0CandTree->Branch("f3PartPi0PhotConvz",&f3PartPi0PhotConvz,"f3PartPi0PhotConvz/F");
    f3PartPi0CandTree->Branch("f3PartPrimVertX",&f3PartPrimVertX,"f3PartPrimVertX/F");
    f3PartPi0CandTree->Branch("f3PartPrimVertY",&f3PartPrimVertY,"f3PartPrimVertY/F");
    f3PartPi0CandTree->Branch("f3PartPrimVertZ",&f3PartPrimVertZ,"f3PartPrimVertZ/F");
    f3PartPi0CandTree->Branch("f3PartPi0ElecPosx",&f3PartPi0ElecPosx,"f3PartPi0ElecPosx/F");
    f3PartPi0CandTree->Branch("f3PartPi0ElecPosy",&f3PartPi0ElecPosy,"f3PartPi0ElecPosy/F");
    f3PartPi0CandTree->Branch("f3PartPi0ElecPosz",&f3PartPi0ElecPosz,"f3PartPi0ElecPosz/F");

    // Create Particle Arrays to store selected Protons and Photons
    fProtonArray       = new std::vector<int>;
    fAntiProtonArray   = new std::vector<int>;
    fConvPhotonArray   = new std::vector<int>;
    fElectronArray     = new std::vector<int>;
    fPositronArray     = new std::vector<int>;
    fElectronPairArray = new std::vector<int>;
    fPositronPairArray = new std::vector<int>;
    
    //Save Particle Masses and other constants for later use
    cElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();      
    cProtonMass   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    
    cPionMass     = TDatabasePDG::Instance()->GetParticle(211)->Mass();     
    cPi0Mass      = TDatabasePDG::Instance()->GetParticle(111)->Mass();     
    c = 2.99792457999999984e-02; // [cm/ps]

/**************************Histograms********************************/

    //Event related    
    TH1F* fHistVertexZ             = new TH1F("fHistVertexZ", "Z Vertex Position;z [cm];Counts/mm", 400, -20, 20);
    TH1F* fHistCentrality          = new TH1F("fHistCentrality", "Centrality Percentile;Centrality [%];Counts/%", 100, 0, 100);
    TH1F* fHistEventCounter        = new TH1F("fHistEventCounter", "Event Counter", 1, 0.5, 1.5);

    //Counters
    TH1F* fHistMCCounter           = new TH1F("fHistMCCounter", "Monte Carlo Counter", 27, 0.5, 27.5);
    TH2F* fHistCandperEvent        = new TH2F("fHistCandperEvent", "Particle Candidates per Event;#pi^{0};p",21,-0.5,20.5,21,-0.5,20.5);

    //MC Kinematics
    TH1F* fHistMCPhotonPt          = new TH1F("fHistMCPhotonPt","Transverse momentum of MC Photons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistMCSigmaPt           = new TH1F("fHistMCSigmaPt","Transverse momentum of MC #Sigma^{+};p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistV0PhotonPt          = new TH1F("fHistV0PhotonPt","Transverse momentum of V0 Photons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistMCElectronPt        = new TH1F("fHistMCElectronPt","Transverse momentum of MC Electrons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistRecoElectronPt      = new TH1F("fHistRecoElectronPt","Transverse momentum of reconstructed Electrons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistMCSigmaPhotonPt     = new TH1F("fHistMCSigmaPhotonPt","Transverse momentum of MC Photons from Sigma;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH1F* fHistV0SigmaPhotonPt     = new TH1F("fHistV0SigmaPhotonPt","Transverse momentum of V0 Photons from Sigma;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,20);
    TH2F* fHistSigmaptvsept        = new TH2F("fHistSigmaptvsept","p_{t} of Sigma vs p_{t} of electrons;p_{T} [GeV/c];p_{T} [GeV/c]",200,0,20,200,0,20);
    TH2F* fHistSigmaptvspi0pt      = new TH2F("fHistSigmaptvspi0pt","p_{t} of Sigma vs p_{t} of #pi^{0};#Sigma^{+} p_{T} [GeV/c];#pi^{0} p_{T} [GeV/c]",200,0,20,200,0,20);

    //Proton Track Quality
    TH1F* fHistNClusterProt        = new TH1F("fHistNClusterProt","Number of TPC Clusters",160,0,160);

    //Proton PID
    TH2F* fHistPvsBeta             = new TH2F("fHistPvsBeta","momentum vs TOF #beta;p [GeV/c];TOF #beta",500,0,5,500,0.1,1.1);                                
    TH2F* fHistPvsdEdx             = new TH2F("fHistPvsdEdx","momentum vs TPC dE/dx;p [GeV/c];TPC #frac{dE}{dx} (a.u.)",500,0,10,500,0,500);
    TH2F* fHistPvsBetaBothCuts     = new TH2F("fHistPvsBetaBothCuts","momentum vs TOF #beta with TPC and TOF Cut;p [GeV/c];TOF #beta",500,0,5,500,0.1,1.1);              
    TH2F* fHistPvsdEdxBothCuts     = new TH2F("fHistPvsdEdxBothCuts","momentum vs TPC dE/dx with TPC and TOF Cut;p [GeV/c];TPC #frac{dE}{dx} (a.u.)",500,0,10,500,0,500);

    //Proton Kinematics and Topology
    TH1F* fHistRealProtonPt        = new TH1F("fHistRealProtonPt","Transverse momentum of MC Protons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistProtonfromSigmaPt   = new TH1F("fHistProtonfromSigmaPt","Transverse momentum of MC Protons from #Sigma^{+};p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistSelectProtonPt      = new TH1F("fHistSelectProtonPt","Transverse momentum of selected Protons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistRealSelectProtonPt  = new TH1F("fHistRealSelectProtonPt","Transverse momentum of selected MC Protons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistSelectProtfromSigPt = new TH1F("fHistSelectProtfromSigPt","Transverse momentum of selected MC Protons from #Sigma^{+};p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistRealProtonDCAxy     = new TH1F("fHistRealProtonDCAxy","Proton DCAxy;DCA_{xy} [cm];Counts/(25#mum)",8000,-10,10);
    TH1F* fHistRealProtonDCAz      = new TH1F("fHistRealProtonDCAz","Proton DCAz;DCA_{z} [cm];Counts/(25#mum)",8000,-10,10);
    TH1F* fHistProtonfromSigDCAxy  = new TH1F("fHistProtonfromSigDCAxy","Proton from #Sigma^{+} DCAxy;DCA_{xy} [cm];Counts/(25#mum)",8000,-10,10);
    TH1F* fHistProtonfromSigDCAz   = new TH1F("fHistProtonfromSigDCAz","Proton from #Sigma^{+} DCAz;DCA_{z} [cm];Counts/(25#mum)",8000,-10,10);

    //V0 kinematics
    TH1F* fHistSelectedPhotPt      = new TH1F("fHistSelectedPhotPt","Transverse momentum of selected Photons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistRealPhotonPt        = new TH1F("fHistRealPhotonPt","Transverse momentum of real Photons;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistSelectRealPhotonPt  = new TH1F("fHistSelectRealPhotonPt","Transverse momentum of real Photons with Cuts;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistV0DaughtPt          = new TH1F("fHistV0DaughtPt","Transverse momentum of V0 Daughters;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);
    TH1F* fHistSelectedV0DaughtPt  = new TH1F("fHistSelectedV0DaughtPt","Transverse momentum of V0 Photon Daughters;p_{T} [GeV/c];Counts/(50 MeV/c)",200,0,10);

    //V0 topology
  	TH1F* fHistCosPA               = new TH1F("fHistCosPA","Cosine of pointing angle;cos(PA);Counts/0.0005",60,0.98,1.01);
  	TH1F* fHistDCAV0Daught         = new TH1F("fHistDCAV0Daught","DCA between photon daughters;DCA [cm];Counts/(0.5mm)",200,0,10);
    TH1F* fHistDCADaughPV          = new TH1F("fHistDCADaughPV","DCA to PV;DCA [cm];Counts/mm",200,0,20);
    TH1F* fHistdeltatheta          = new TH1F("fHistdeltatheta","#Delta#Theta;#Delta#Theta;Counts/(0.0157)",200,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHisttotangle            = new TH1F("fHisttotangle","Total V0 Opening Angle;#xi;Counts/(0.0157)",200,0,TMath::Pi());
    TH2F* fHistVertexPos           = new TH2F("fHistVertexPos","Position of photon conversion vertex;x [cm];y [cm];",480,-120,120,480,-120,120);    
    TH1F* fHistTransverseRadius    = new TH1F("fHistTransverseRadius","Conversion Radius of Photon;r_{xy} [cm];Counts/(5 mm)",600,0,300);

    //MC V0 topology
	  TH1F* fHistCosPASigma          = new TH1F("fHistCosPASigma","Cosine of pointing angle, Photons from #Sigma^{+};cos(PA);Counts/0.0005",60,0.98,1.01);
    TH1F* fHistMCPhotondeltatheta  = new TH1F("fHistMCPhotondeltatheta","#Delta#Theta of MC Photons;#Delta#Theta;Counts/(0.0157)",200,-TMath::Pi()/2,TMath::Pi()/2);
    TH1F* fHistMCPhotontotangle    = new TH1F("fHistMCPhotontotangle","Total V0 Opening Angle of MC Photons;#xi;Counts/(0.0157)",200,0,TMath::Pi());
    TH1F* fHistConvRadius          = new TH1F("fHistConvRadius","Conversion Radius;r [cm];Counts/(5 mm)",1000,0,500);
    TH2F* fHistpVsMCConvRadius     = new TH2F("fHistpVsMCConvRadius","Conversion Radius of MC Photon vs momentum;r_{xy} [cm];p [GeV/c]",1000,0,500,500,0,5);

    //V0 track quality
    TH1F* fHistNClusterV0          = new TH1F("fHistNClusterV0","Number of TPC Clusters of V0 Daughters",160,0,160);
    
    //TPC PID
    TH2F* fHistPdEdxDaught         = new TH2F("fHistPdEdxDaught", "p VS. TPC dEdx Daughters;p [GeV/c];TPC #frac{dE}{dx} (a.u.)",500,0,5,500,0,200);
    TH2F* fHistElecNsigTPC         = new TH2F("fHistElecNsigTPC","Electron momentum vs N sigma TPC;p [GeV/c];n_{#sigma,TPC}",500,0,5,500,-10,10);
    
    //Armenteros-Podolanski
    TH2F* fHistArmPod              = new TH2F("fHistArmPod", "Armenteros-Podolanski Plot;#alpha;q_{t} [GeV/c]",200,-1.0,1.0,200,-0.05,0.4);
    TH2F* fHistArmPodPhotons       = new TH2F("fHistArmPodPhotons", "Armenteros-Podolanski Plot Photons;#alpha;q_{t} [GeV/c]",200,-1.0,1.0,200,-0.01,0.3);

    //Electron Track Quality
    TH1F* fHistNClusterElec        = new TH1F("fHistNClusterElec","Number of TPC Clusters",160,0,160);

    //Electron PID
    TH2F* fHistElecNsigTPCall      = new TH2F("fHistElecNsigTPCall","e^{-} p vs n_{#sigma,TPC}, No PID;p [GeV/c];n_{#sigma,TPC}",500,0,5,500,-10,10);
    TH2F* fHistElecNsigTPCTPCsamp  = new TH2F("fHistElecNsigTPCTPCsamp","e^{-} p vs n_{#sigma,TPC}, TPC Sample;p [GeV/c];n_{#sigma,TPC}",500,0,5,500,-10,10);
    TH2F* fHistElecNsigTPCTOFsamp  = new TH2F("fHistElecNsigTPCTOFsamp","e^{-} p vs n_{#sigma,TPC}, TOF Sample;p [GeV/c];n_{#sigma,TPC}",500,0,5,500,-10,10);
    TH2F* fHistElecNsigTPCcomb     = new TH2F("fHistElecNsigTPCcomb","e^{-} p vs n_{#sigma,TPC}, Combined Sample;p [GeV/c];n_{#sigma,TPC}",500,0,5,500,-10,10);

    //invariant masses
    TH1F* fHistV0photonmass        = new TH1F("fHistV0photonmass", "Photon Mass;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 1000, 0, 1);
    TH1F* fHistpi0mass             = new TH1F("fHistpi0mass", "#pi^{0} Mass;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHistpi0massKF           = new TH1F("fHistpi0massKF", "#pi^{0} Mass with KFParticle;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHistMCpi0mass           = new TH1F("fHistMCpi0mass", "MC #pi^{0} Mass;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHistMCsigmamass         = new TH1F("fHistMCsigmamass", "MC #Sigma^{+} Mass;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 500, 1, 1.5);
    TH1F* fHistMCdeltamass         = new TH1F("fHistMCdeltamass", "MC #Delta^{+} Mass;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 500, 1, 1.5);
    TH1F* fHist3Partpi0mass        = new TH1F("fHist3Partpi0mass", "#pi^{0} Mass (3 electrons);m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHist3Partpi0massKF      = new TH1F("fHist3Partpi0massKF", "#pi^{0} Mass (3 electrons) with KFParticle;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHist3Partpi0massMC      = new TH1F("fHist3Partpi0massMC", "#pi^{0} Mass (3 electrons);m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHist3Partpi0massMChighpt= new TH1F("fHist3Partpi0massMChighpt", "#pi^{0} Mass (3 electrons);m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH2F* fHist3Partpi0massvspt    = new TH2F("fHist3Partpi0massvspt", "#pi^{0} Mass (3 electrons) vs. single electron pt;m_{inv} [GeV/c^{2}];p_{t} [GeV/c]", 700, 0, 0.7, 100, 0, 1);

    //KFParticle
    TH1F* fHisteephotmassTL        = new TH1F("fHisteephotmassTL", "Photon Mass with TLorentzvector;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 1000, 0, 1);
    TH1F* fHisteephotmassKF        = new TH1F("fHisteephotmassKF", "Photon Mass with KFParticle;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 1000, 0, 1);
    TH1F* fHisteephotmassKFMC      = new TH1F("fHisteephotmassKFMC", "MC Photon Mass with KFParticle;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 1000, 0, 1);
    TH1F* fHisteephotmassKFNoCut   = new TH1F("fHisteephotmassKFNoCut", "Photon Mass with KFParticle without Cuts;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 1000, 0, 1);
    TH1F* fHisteepi0massTL         = new TH1F("fHisteepi0massTL", "#pi^{0} Mass with TL from e^+-e^--pairs;m_{inv} [GeV/c^{2}];Counts/(MeV/c^{2})", 700, 0, 0.7);
    TH1F* fHisteeopenangle         = new TH1F("fHisteeopenangle","Total Opening Angle of e^+-e^--pair;#xi;Counts/(0.0157)",200,0,TMath::Pi());
    TH2F* fHistArmPodKF            = new TH2F("fHistArmPodKF", "Armenteros-Podolanski Plot with KFParticle;#alpha;q_{t} [GeV/c]",200,-1.0,1.0,200,-0.05,0.4);
	  TH1F* fHistKFeeDCA             = new TH1F("fHistKFeeDCA","DCA between photon daughters with KF;DCA [cm];Counts/(0.5mm)",200,0,10);
    TH1F* fHisteeopenangleMC       = new TH1F("fHisteeopenangleMC","Total Opening Angle of e^+-e^--pair;#xi;Counts/(0.0157)",200,0,TMath::Pi());
    TH2F* fHistArmPodKFMC          = new TH2F("fHistArmPodKFMC", "Armenteros-Podolanski Plot with KFParticle;#alpha;q_{t} [GeV/c]",200,-1.0,1.0,200,-0.05,0.4);
  	TH1F* fHistKFeeDCAMC           = new TH1F("fHistKFeeDCAMC","DCA between photon daughters with KF;DCA [cm];Counts/(0.5mm)",200,0,10);

    //MC Sigma Topologyc
    TH1F* fHistPi0VertexvsMC       = new TH1F("fHistPi0VertexvsMC","#pi^{0} KF Decay Radius - MC Decay Radius;#Deltar [cm];Counts/(mm)", 400, -20, 20);
    TH1F* fHist3PartPi0VertexvsMC  = new TH1F("fHist3PartPi0VertexvsMC","#pi^{0} KF Decay Radius - MC Decay Radius;#Deltar [cm];Counts/(mm)", 400, -20, 20);
    TH1F* fHistMCSigmaPA           = new TH1F("fHistMCSigmaPA","#Sigma^{+} PA;PA [rad];Counts/(0.0157)",200,0,TMath::Pi());
    TH1F* fHistMC4PartSigmaPA      = new TH1F("fHistMC4PartSigmaPA","#Sigma^{+} PA;PA [rad];Counts/(0.0157)",200,0,TMath::Pi());
    TH1F* fHistPi0VertexMC         = new TH1F("fHistPi0VertexMC","#pi^{0} MC Decay Radius;r [cm];Counts/(mm)", 200, 0, 20);

    //Miscellaneous
    /*** EMPTY ***/

    //Alphanumeric Histogram Labels
    fHistMCCounter->GetXaxis()->SetBinLabel(1,"Events");
    fHistMCCounter->GetXaxis()->SetBinLabel(2,"p");
    fHistMCCounter->GetXaxis()->SetBinLabel(3,"#barp");
    fHistMCCounter->GetXaxis()->SetBinLabel(4,"#Delta^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(5,"#bar#Delta^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(6,"#Sigma^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(7,"#bar#Sigma^{-}");    
    fHistMCCounter->GetXaxis()->SetBinLabel(8,"pi^{0}");
    fHistMCCounter->GetXaxis()->SetBinLabel(9,"#gamma");
    fHistMCCounter->GetXaxis()->SetBinLabel(10,"#gamma from pi^{0}");
    fHistMCCounter->GetXaxis()->SetBinLabel(11,"#gamma from #Sigma^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(12,"#gamma->e^{+}e^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(13,"#gamma from pi^{0}->e^{+}e^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(14,"#gamma from #Sigma^{+}->e^{+}e^{-}");
    fHistMCCounter->GetXaxis()->SetBinLabel(15,"V0");
    fHistMCCounter->GetXaxis()->SetBinLabel(16,"On-fly V0");
    fHistMCCounter->GetXaxis()->SetBinLabel(17,"V0 is #gamma");
    fHistMCCounter->GetXaxis()->SetBinLabel(18,"#gamma is from #Sigma^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(19,"e+-e--pairs");
    fHistMCCounter->GetXaxis()->SetBinLabel(20,"pair is #gamma");
    fHistMCCounter->GetXaxis()->SetBinLabel(21,"#gamma is from #Sigma^{+}");
    fHistMCCounter->GetXaxis()->SetBinLabel(22,"Reconstructed MC #pi^{0} V0");
    fHistMCCounter->GetXaxis()->SetBinLabel(23,"Reconstructed MC #Sigma V0");
    fHistMCCounter->GetXaxis()->SetBinLabel(24,"MC Electrons");
    fHistMCCounter->GetXaxis()->SetBinLabel(25,"MC Positrons");
    fHistMCCounter->GetXaxis()->SetBinLabel(26,"Reconstructed Electrons");
    fHistMCCounter->GetXaxis()->SetBinLabel(27,"Reconstructed Positrons");

    for(Int_t i=0; i<21; i++){
    fHistCandperEvent->GetXaxis()->SetBinLabel(i+1,Form("%d", i));
    fHistCandperEvent->GetYaxis()->SetBinLabel(i+1,Form("%d", i));
    }

    fHistEventCounter->GetXaxis()->SetBinLabel(1,"Events");

    //Add Histograms to Output List
    //Event related
    fOutputList->Add(fHistVertexZ);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistEventCounter);
    //Counters
    fOutputList->Add(fHistMCCounter);
    fOutputList->Add(fHistCandperEvent);
    //MC Kinematics
    fOutputList->Add(fHistMCPhotonPt);
    fOutputList->Add(fHistMCSigmaPt);
    fOutputList->Add(fHistV0PhotonPt);
    fOutputList->Add(fHistMCElectronPt);
    fOutputList->Add(fHistRecoElectronPt);
    fOutputList->Add(fHistMCSigmaPhotonPt);
    fOutputList->Add(fHistV0SigmaPhotonPt);
    fOutputList->Add(fHistSigmaptvsept);
    fOutputList->Add(fHistSigmaptvspi0pt);
    //Proton Track Quality
    fOutputList->Add(fHistNClusterProt);
    //Proton PID
    fOutputList->Add(fHistPvsBeta);
    fOutputList->Add(fHistPvsdEdx);
    fOutputList->Add(fHistPvsBetaBothCuts);
    fOutputList->Add(fHistPvsdEdxBothCuts);
    //Proton Kinematics and Topology
    fOutputList->Add(fHistRealProtonPt);
    fOutputList->Add(fHistProtonfromSigmaPt);
    fOutputList->Add(fHistSelectProtonPt);
    fOutputList->Add(fHistRealSelectProtonPt);
    fOutputList->Add(fHistSelectProtfromSigPt);
    fOutputList->Add(fHistRealProtonDCAxy);
    fOutputList->Add(fHistRealProtonDCAz);
    fOutputList->Add(fHistProtonfromSigDCAxy);           
    fOutputList->Add(fHistProtonfromSigDCAz);       
    //V0 kinematics
    fOutputList->Add(fHistSelectedPhotPt);
    fOutputList->Add(fHistRealPhotonPt);
    fOutputList->Add(fHistSelectRealPhotonPt);
    fOutputList->Add(fHistV0DaughtPt);
    fOutputList->Add(fHistSelectedV0DaughtPt);
    //V0 topology
    fOutputList->Add(fHistCosPA);
    fOutputList->Add(fHistDCAV0Daught);
    fOutputList->Add(fHistDCADaughPV);
    fOutputList->Add(fHistdeltatheta);
    fOutputList->Add(fHisttotangle);
    fOutputList->Add(fHistVertexPos);
    fOutputList->Add(fHistTransverseRadius);
    //MC V0 topology
    fOutputList->Add(fHistCosPASigma);
    fOutputList->Add(fHistMCPhotondeltatheta);
    fOutputList->Add(fHistMCPhotontotangle);
    fOutputList->Add(fHistConvRadius);
    fOutputList->Add(fHistpVsMCConvRadius);
    //V0 track quality
    fOutputList->Add(fHistNClusterV0);
    //TPC PID
    fOutputList->Add(fHistPdEdxDaught);
    fOutputList->Add(fHistElecNsigTPC);
    //Armenteros-Podolanski
    fOutputList->Add(fHistArmPod);
    fOutputList->Add(fHistArmPodPhotons);
    //Electron Track Quality
    fOutputList->Add(fHistNClusterElec);
    //Electron PID
    fOutputList->Add(fHistElecNsigTPCall);
    fOutputList->Add(fHistElecNsigTPCTPCsamp);
    fOutputList->Add(fHistElecNsigTPCTOFsamp);
    fOutputList->Add(fHistElecNsigTPCcomb);
    //invariant masses
    fOutputList->Add(fHistV0photonmass);
    fOutputList->Add(fHistpi0mass);
    fOutputList->Add(fHistpi0massKF);
    fOutputList->Add(fHistMCpi0mass);
    fOutputList->Add(fHistMCsigmamass);
    fOutputList->Add(fHistMCdeltamass);
    fOutputList->Add(fHist3Partpi0mass);
    fOutputList->Add(fHist3Partpi0massKF);
    fOutputList->Add(fHist3Partpi0massMC);
    fOutputList->Add(fHist3Partpi0massMChighpt);
    fOutputList->Add(fHist3Partpi0massvspt);
    //KFParticle
    fOutputList->Add(fHisteephotmassTL);
    fOutputList->Add(fHisteephotmassKF);
    fOutputList->Add(fHisteephotmassKFMC);
    fOutputList->Add(fHisteephotmassKFNoCut);
    fOutputList->Add(fHisteepi0massTL);
    fOutputList->Add(fHisteeopenangle);
    fOutputList->Add(fHistArmPodKF);
    fOutputList->Add(fHistKFeeDCA);
    fOutputList->Add(fHisteeopenangleMC);
    fOutputList->Add(fHistArmPodKFMC);
    fOutputList->Add(fHistKFeeDCAMC);
    //MC Sigma Topology
    fOutputList->Add(fHistPi0VertexvsMC);     
    fOutputList->Add(fHist3PartPi0VertexvsMC);
    fOutputList->Add(fHistMCSigmaPA);         
    fOutputList->Add(fHistMC4PartSigmaPA);    
    fOutputList->Add(fHistPi0VertexMC);

/**************************************************************************/
    PostData(1, fOutputList);         
    PostData(2, fSigmaCandTree);         
    PostData(3, f4PartSigmaCandTree);         
    PostData(4, f3PartPi0CandTree);         

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
  Float_t centrality = 0;
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(aodEvent->FindListObject("MultSelection"));   //Get Mult Selection..
  if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");                                //..and retrieve centrality
  FillHistogram("fHistCentrality",centrality);

  FillHistogram("fHistMCCounter",1);  //Event Counter 
  FillHistogram("fHistEventCounter",1);    //Event Counter 

/************************Start Event Processing**************************************/

  Bool_t fDebug = kFALSE;

  if(fDebug) cout << "Start of Event\n";

  //Process Protons
  if(fDebug) cout << "Processing Protons\n";
  FillProtonArray();

  //Process Electrons
  if(fDebug) cout << "Processing Electrons\n";
  //FillElectronArray();

  //Process MC Particles
  if(fDebug&&isMonteCarlo) cout << "Processing MC Particles\n";
  //if(isMonteCarlo) ProcessMCParticles();

  //Process V0s
  if(fDebug) cout << "Processing V0s\n";
  FillV0PhotonArray();

  //Pair selected Electrons
  if(fDebug) cout << "Pairing Electrons\n";
  //PairElectrons();

  //Pair Photons and Electrons
  if(fDebug) cout << "Pairing Photons and Electrons\n";
  //PairPhotonandElectron();

  //Reconstruct Pi0 and Sigma+
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas\n";
  ReconstructParticles();

  //Reconstruct Pi0 from ee pairs
  if(fDebug) cout << "Reconstructing Pi0s and Sigmas from Electrons\n";
  //ReconstructParticlesfromee();

  if(fDebug) cout << "End of Event\n";

/************************End of Event Processing**************************************/

    PostData(1, fOutputList); // stream the analysis results of the current event to output manager        
    PostData(2, fSigmaCandTree);         
    PostData(3, f4PartSigmaCandTree);
    PostData(4, f3PartPi0CandTree);         
    return;

}//end of UserExec()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillProtonArray() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fProtonArray->clear();
  fAntiProtonArray->clear();
  Int_t countProton = 0;
  Int_t countAntiProton = 0;
  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  //Loop for Proton Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {

    //Initialisation of local variables
    Bool_t isReallyProton = kFALSE; //From MC
    Bool_t isProtonfromSigma = kFALSE; //From MC

    Bool_t   isTPCProton = kFALSE;
    Bool_t   isTOFProton = kFALSE;

    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack) {
      AliWarning("No AOD Track!");
      continue;
    }

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

    //Acceptance Cut
    if(TMath::Abs(aodTrack->Eta()) > fMaxProtEta) continue; 

    Double_t charge = aodTrack->Charge();
    Double_t dEdx   = aodTrack->GetTPCsignal();
    Double_t p      = aodTrack->GetTPCmomentum();
    Double_t pt     = aodTrack->Pt();
    Int_t    nTPCClustProt = aodTrack->GetTPCNcls();

    const Double_t len = aodTrack->GetIntegratedLength();
    const Double_t tim = aodTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(aodTrack->GetTPCmomentum());
    Double_t beta = -1;
    if(tim != 0.) beta = len / (tim * c);

    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton))<fMaxNsigProtTPC) isTPCProton = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton))<fMaxNsigProtTOF) isTOFProton = kTRUE;

    Float_t  DCAxy = -999., DCAz = -999.;
    aodTrack->GetImpactParameters(DCAxy,DCAz);

    //Check MC Truth
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        if(mcPart->GetPdgCode()==2212 || mcPart->GetPdgCode()==-2212){
          isReallyProton = kTRUE;
          if(mcPart->GetMother()!=-1){
            AliAODMCParticle* ProtonMotherPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mcPart->GetMother()));
            if(ProtonMotherPart){
                if(ProtonMotherPart->GetPdgCode()==3222 || ProtonMotherPart->GetPdgCode()==-3222) isProtonfromSigma = kTRUE;
            }//MC Mother exists
          }//Proton has Mother
        }//is Proton or Anti-Proton
      }//MC Particle exists
    }//MC treatment

    if(isReallyProton){
        FillHistogram("fHistRealProtonPt",pt);  
        FillHistogram("fHistRealProtonDCAxy",DCAxy);   
        FillHistogram("fHistRealProtonDCAz",DCAz);    
    }
    if(isProtonfromSigma){
        FillHistogram("fHistProtonfromSigmaPt",pt);  
        FillHistogram("fHistProtonfromSigDCAxy",DCAxy);   
        FillHistogram("fHistProtonfromSigDCAz",DCAz);    
    }

    // Check Track Quality Histograms
    FillHistogram("fHistNClusterProt",nTPCClustProt);
    if(nTPCClustProt<fMinTPCClustProt) continue;

    // Filling PID Histograms
    FillHistogram("fHistPvsBeta",p,beta);        
    FillHistogram("fHistPvsdEdx",p,dEdx);        

    // Proton PID Cut
    if(!isTPCProton) continue;
    if(p > fMaxpOnlyTPCPID && !isTOFProton) continue;

    // Filling PID Histograms after PID Cut
    FillHistogram("fHistPvsBetaBothCuts",p,beta);
    FillHistogram("fHistPvsdEdxBothCuts",p,dEdx);

    FillHistogram("fHistSelectProtonPt",pt);           
    if(isReallyProton) FillHistogram("fHistRealSelectProtonPt",pt);  

    if(isProtonfromSigma) FillHistogram("fHistSelectProtfromSigPt",pt);  

    // Proton Kinematic Cut
    if(pt < fMinProtpt || pt > fMaxProtpt) continue;

    // Store (Anti-)Proton candidates after selection
    if(charge > 0.) { // Protons
      fProtonArray->push_back(iTrack);
      countProton++;
    }// End of Charge > 0
    else if(charge < 0.) { // Anti-Protons
      fAntiProtonArray->push_back(iTrack);
      countAntiProton++;
    } // End of Charge < 0

  } // End of proton track loop

return;

}//End of FillProtonArray  

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillElectronArray() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear Proton and Antiproton Arrays and reset counters
  fElectronArray->clear();
  fPositronArray->clear();
  Int_t countElectrons = 0;
  Int_t countPositrons = 0;

  //Save IDs of Tracks to avoid double counting
  std::vector<int> IDvector;
  IDvector.clear();

  //Loop for Proton Selection
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) {
  
    //Initialisation of local Bools
    Bool_t isTPCElectron = kFALSE;
    Bool_t isTPCMuon     = kFALSE;
    Bool_t isTPCPion     = kFALSE;
    Bool_t isTPCKaon     = kFALSE;
    Bool_t isTPCProton   = kFALSE;
    Bool_t isTOFElectron = kFALSE;
    Bool_t isGoodTPCElectron = kFALSE;
    Bool_t isGoodTOFElectron = kFALSE;
    Bool_t isReallyElectron = kFALSE;
        
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!aodTrack) {
      AliWarning("No AOD Track!");
      continue;
    }

    if(TMath::Abs(aodTrack->Eta()) > fMaxElecEta) continue; // Acceptance Cut

    //Check MC Truth
    if(isMonteCarlo){
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(aodTrack->GetLabel())));
      if(mcPart){
        if(mcPart->GetPdgCode()==11 || mcPart->GetPdgCode()==-11){isReallyElectron = kTRUE; FillHistogram("fHistRecoElectronPt",aodTrack->Pt());}
      }//MC Particle exists
    }//MC treatment

    // Track Quality Check
    FillHistogram("fHistNClusterElec",TMath::Abs(aodTrack->GetTPCNcls()));
    if(TMath::Abs(aodTrack->GetTPCNcls()) < fMinTPCClustElec) continue;

    //Check for double counted Tracks if no filterbit is used
    Int_t nIDs = IDvector.size();
    Bool_t isdouble = kFALSE;
    for(Int_t iID = 0; iID<nIDs; iID++){
      if(aodTrack->GetID()==IDvector[iID]) isdouble = kTRUE; 
    }
    if(isdouble) continue;
    else IDvector.push_back(aodTrack->GetID());

    Double_t charge = aodTrack->Charge();
    Double_t p      = aodTrack->GetTPCmomentum();
    Double_t pt     = aodTrack->GetTPCmomentum();
    
    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kElectron);
    Double_t nSigmaTOFelectron = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kElectron));
    Double_t nSigmaTPCmuon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kMuon));
    Double_t nSigmaTPCpion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
    Double_t nSigmaTPCkaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
    Double_t nSigmaTPCproton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));

    //Electron PID
    if(TMath::Abs(nSigmaTPCelectron)<fMaxNsigElecTPC) isTPCElectron = kTRUE;
    if(TMath::Abs(nSigmaTOFelectron)<fMaxNsigElecTOF) isTOFElectron = kTRUE;

    //Hadron + Muon Rejection
    //if(nSigmaTPCmuon<fNsigHadronreject) isTPCMuon = kTRUE; //was -1.75 .. 1
    if(nSigmaTPCpion<fNsigHadronreject) isTPCPion = kTRUE;   //was -2.5 .. 4
    //if(nSigmaTPCkaon<fNsigHadronreject) isTPCKaon = kTRUE;   //was -2.25 .. 3
    //if(nSigmaTPCproton<fNsigHadronreject) isTPCProton = kTRUE;   //was -2.25 .. 3.5

    if(isTPCElectron&&!isTPCPion&&!isTPCKaon&&!isTPCProton) isGoodTPCElectron = kTRUE; 
    if(isTOFElectron&&isTPCElectron&&!isTPCPion) isGoodTOFElectron = kTRUE; 

    // Filling TPC N Sigma Histograms
    FillHistogram("fHistElecNsigTPCall",p,nSigmaTPCelectron);
    if(isGoodTPCElectron) FillHistogram("fHistElecNsigTPCTPCsamp",p,nSigmaTPCelectron);        
    if(isGoodTOFElectron) FillHistogram("fHistElecNsigTPCTOFsamp",p,nSigmaTPCelectron);
    if(isGoodTPCElectron||isGoodTOFElectron) FillHistogram("fHistElecNsigTPCcomb",p,nSigmaTPCelectron);

    // Electron PID
    if(!isGoodTPCElectron&&!isGoodTOFElectron) continue;

    // Store Electron and Positron candidates after selection
    if(charge < 0.) { // Electrons
      fElectronArray->push_back(iTrack);
      FillHistogram("fHistMCCounter",26);
      countElectrons++;
    }// End of Charge > 0
    else if(charge > 0.) { // Positrons
      fPositronArray->push_back(iTrack);
      countPositrons++;
      FillHistogram("fHistMCCounter",27);
    } // End of Charge < 0

  } // End of electron track loop

return;

}//End of FillElectronArray()

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
    if(MCPartPDGCode == 2212) {FillHistogram("fHistMCCounter",2);}
    if(MCPartPDGCode == -2212){FillHistogram("fHistMCCounter",3);}
    if(MCPartPDGCode == 2214) {FillHistogram("fHistMCCounter",4);}
    if(MCPartPDGCode == -2214){FillHistogram("fHistMCCounter",5);}
    if(MCPartPDGCode == 3222) {FillHistogram("fHistMCCounter",6);
    FillHistogram("fHistMCSigmaPt",mcPart->Pt()); 
    }
    if(MCPartPDGCode == -3222){FillHistogram("fHistMCCounter",7);
    FillHistogram("fHistMCSigmaPt",mcPart->Pt());       
    }
    if(MCPartPDGCode == 111)  {FillHistogram("fHistMCCounter",8);}
    if(MCPartPDGCode == 11)  {FillHistogram("fHistMCCounter",24);
    FillHistogram("fHistMCElectronPt",mcPart->Pt()); 
    }
    if(MCPartPDGCode == -11)  {FillHistogram("fHistMCCounter",25);
    FillHistogram("fHistMCElectronPt",mcPart->Pt()); 
    }

    if(MCPartPDGCode == 22){ //Photons

      FillHistogram("fHistMCCounter",9);

      Bool_t isfromPi0 = kFALSE, isfromSigma = kFALSE;
      Double_t Sigmapt = -1;
      AliAODMCParticle* PhotonMother = NULL;
      if(mcPart->GetMother()!=-1) PhotonMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetMother())));
      if(PhotonMother){
        if(PhotonMother->GetPdgCode()==111){
          FillHistogram("fHistMCCounter",10);
          isfromPi0 = kTRUE;
          AliAODMCParticle* Pi0Mother = NULL;
          if(mcPart->GetMother()!=-1) Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(PhotonMother->GetMother())));
          if(Pi0Mother){
            if(TMath::Abs(Pi0Mother->GetPdgCode())==3222){isfromSigma = kTRUE; FillHistogram("fHistMCCounter",11); Sigmapt=Pi0Mother->Pt(); FillHistogram("fHistSigmaptvspi0pt",Sigmapt,PhotonMother->Pt());}
          }
        }
      }

      FillHistogram("fHistMCPhotonPt",mcPart->Pt()); 
      if(isfromSigma) FillHistogram("fHistMCSigmaPhotonPt",mcPart->Pt());

      if(mcPart->GetNDaughters()!=2) continue;  //Check if Photon has 2 Daughters (i.e. if it was converted to e+ e-)
                    
      AliAODMCParticle* FirstDaughtParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetDaughterFirst())));
      AliAODMCParticle* LastDaughtParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(mcPart->GetDaughterLast())));
      if(!FirstDaughtParticle || !LastDaughtParticle) continue;

      //Check if Photon Daughters are e+ and e-
      if(FirstDaughtParticle->GetPdgCode()!=11 && FirstDaughtParticle->GetPdgCode()!=-11) continue; 
      if(LastDaughtParticle->GetPdgCode()!=11 && LastDaughtParticle->GetPdgCode()!=-11) continue;      
        
      //Check conversion vertex (aka production vertex of e+/e-)
      TVector3 prdvtx1(FirstDaughtParticle->Xv(),FirstDaughtParticle->Yv(),FirstDaughtParticle->Zv());   
      TVector3 prdvtx2(FirstDaughtParticle->Xv(),FirstDaughtParticle->Yv(),FirstDaughtParticle->Zv());   
      TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
      FillHistogram("fHistConvRadius",prdvtx.Perp());
      FillHistogram("fHistpVsMCConvRadius",prdvtx.Perp(),mcPart->P());
      if(prdvtx.Perp()<180){
        FillHistogram("fHistMCCounter",12);
        if(isfromPi0) FillHistogram("fHistMCCounter",13);
        if(isfromSigma){ FillHistogram("fHistMCCounter",14);
        FillHistogram("fHistSigmaptvsept",Sigmapt,FirstDaughtParticle->Pt());
        FillHistogram("fHistSigmaptvsept",Sigmapt,LastDaughtParticle->Pt());
        }
      }
                
    } // End of Photon treatment
  } //End of MC Particle loop 

return;

}//End of ProcessMCParticles()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::FillV0PhotonArray() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  //Clear V0 Photon Array and reset counter
  fConvPhotonArray->clear();
  Int_t countPhotons = 0;

  Int_t nV0 = aodEvent->GetNumberOfV0s(); //Number of V0s in the event
  if(nV0 == 0) return;      //Return if there is no V0 to be processed

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

    FillHistogram("fHistMCCounter",15); //V0 Counter
    if(!aodV0->GetOnFlyStatus()) continue;    //Use only On-Fly V0s (better Resolution!)
    FillHistogram("fHistMCCounter",16); //Count on-fly V0

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

    if(TMath::Abs(trackN->Eta()) > fMaxDaughtEta) continue; // Acceptance Cut
    if(TMath::Abs(trackP->Eta()) > fMaxDaughtEta) continue;

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

    //Get reconstructed cartesian momentum
    vecN.SetXYZ(aodV0->MomNegX(),aodV0->MomNegY(),aodV0->MomNegZ()); //negative daughter
    vecP.SetXYZ(aodV0->MomPosX(),aodV0->MomPosY(),aodV0->MomPosZ()); //positive daughter 
    vecM.SetXYZ(aodV0->MomV0X(),aodV0->MomV0Y(),aodV0->MomV0Z());    //mother 

    // Get kinematic values
    Double_t ptV0   = aodV0->Pt();
    Double_t ptPos  = trackP->Pt();
    Double_t ptNeg  = trackN->Pt();
    Double_t pV0    = aodV0->P();
    Double_t pPos   = trackP->P();
    Double_t pNeg   = trackN->P();
    Double_t thetaPos = trackP->Theta();
    Double_t thetaNeg = trackN->Theta();
    Double_t dEdxP = trackP->GetTPCsignal();
    Double_t dEdxN = trackN->GetTPCsignal();
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
            if(V0Mother){if(V0Mother->GetPdgCode()==111){
                      
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

    //Filling MC Counter
    if(isReallyPhoton) FillHistogram("fHistMCCounter",17);
    if(isPhotonfromSigma) FillHistogram("fHistMCCounter",18);

    //Filling kinematic Histograms
    if(isReallyPhoton) FillHistogram("fHistRealPhotonPt",ptV0);
    FillHistogram("fHistV0DaughtPt",ptPos);
    FillHistogram("fHistV0DaughtPt",ptNeg);

    //Filling topological Histograms
    FillHistogram("fHistCosPA",cosPointAngle);              
    FillHistogram("fHistDCAV0Daught",dcaV0Daughters);        
    FillHistogram("fHistDCADaughPV",dcaPosToPrimVtx);     
    FillHistogram("fHistDCADaughPV",dcaNegToPrimVtx);     
    FillHistogram("fHistVertexPos",vtxPosV0[0],vtxPosV0[1]);
    FillHistogram("fHistTransverseRadius",Vtxradius);
    FillHistogram("fHistdeltatheta",deltatheta);        
    FillHistogram("fHisttotangle",totangle);

    if(isReallyPhoton) FillHistogram("fHistMCPhotondeltatheta",deltatheta);  
    if(isReallyPhoton) FillHistogram("fHistMCPhotontotangle",totangle);
    if(isPhotonfromSigma) FillHistogram("fHistCosPASigma",cosPointAngle);   

    // Check Track quality and reject poor qualty tracks
    FillHistogram("fHistNClusterV0",nTPCClustNeg);          
    FillHistogram("fHistNClusterV0",nTPCClustPos);          
    if(nTPCClustNeg < fMinTPCClustDaught) continue;
    if(nTPCClustPos < fMinTPCClustDaught) continue;

    //Filling PID Histograms
    FillHistogram("fHistPdEdxDaught",pPos,dEdxP);
    FillHistogram("fHistPdEdxDaught",pNeg,dEdxN);       
    FillHistogram("fHistElecNsigTPC",pNeg,nSigmaTPCelectron);
    FillHistogram("fHistElecNsigTPC",pPos,nSigmaTPCpositron);       

    //Filling Armenteros-Podolanski Histogram
    FillHistogram("fHistArmPod",alpha,qt);            

    // Armenteros-Podolanski Cuts
    if(TMath::Abs(alpha) > fMaxalpha) continue;
    if(TMath::Abs(qt) > fMaxqt) continue;

    // Angle Cut
    if(TMath::Abs(totangle) > fMaxopenangle) continue;
    if(TMath::Abs(deltatheta) > fMaxdeltatheta) continue;

    // CPA Cut
    //if(cosPointAngle < fMinV0CPA) continue;    //Pointing Angle Cut bad for secondary Photons!

    // PID Cut
    if(!isPhotonTPC) continue;    

    FillHistogram("fHistV0photonmass",photonmass);  

    // Inv. Mass Cut
    if(photonmass > fMaxphotonmass) continue;

    //Filling Armenteros-Podolanski Histogram after Selection
    FillHistogram("fHistArmPodPhotons",alpha,qt);            

    //Filling kinematic Histograms after Selection
    FillHistogram("fHistSelectedPhotPt",ptV0);  
    if(isReallyPhoton) FillHistogram("fHistSelectRealPhotonPt",ptV0);  

    if(isReallyPhoton) FillHistogram("fHistV0PhotonPt",ptV0); 
    if(isPhotonfromSigma) FillHistogram("fHistV0SigmaPhotonPt",ptV0);

    FillHistogram("fHistSelectedV0DaughtPt",ptPos);
    FillHistogram("fHistSelectedV0DaughtPt",ptNeg);

    // Store Photon candidates after selection
    fConvPhotonArray->push_back(iV0);
    countPhotons++;

  }//End of V0 Loop

return;

}//End of FillV0PhotonArray()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::PairElectrons() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // PID and MC Bools
  Bool_t   isReallyPhoton = kFALSE;
  Bool_t   isPhotonfromSigma = kFALSE;

  TLorentzVector trackElectron, trackPositron, trackPhoton;
  KFParticle KFElectron, KFPositron;                 //KFParticle for Photon reconstruction

  const Int_t nElectrons = fElectronArray->size();
  const Int_t nPositrons = fPositronArray->size();

  //Clear selected Pair Arrays and reset counter
  fElectronPairArray->clear();
  fPositronPairArray->clear();
  Int_t countGoodPairs = 0;

  for(Int_t i=0; i<nElectrons; i++) {
      
    AliAODTrack *electrack = (AliAODTrack*)aodEvent->GetTrack(fElectronArray->at(i));
    if(!electrack) continue;
    
    for(Int_t j=0; j<nPositrons; j++) {
      
      AliAODTrack *positrack = (AliAODTrack*)aodEvent->GetTrack(fPositronArray->at(j));
      if(!positrack) continue;

      FillHistogram("fHistMCCounter",19);

      // AOD MC treatment
      isReallyPhoton = kFALSE; isPhotonfromSigma = kFALSE;
      if(isMonteCarlo){    
        AliAODMCParticle* V0Part = NULL;
        AliAODMCParticle* NPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(electrack->GetLabel())));
        AliAODMCParticle* PPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(positrack->GetLabel())));
        if(NPart&&PPart){
          if(NPart->GetMother()==PPart->GetMother()&&NPart->GetMother()!=-1){
            V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(NPart->GetMother())));    
            if(V0Part){if(V0Part->GetPdgCode()==22&&NPart->GetPdgCode()==11&&PPart->GetPdgCode()==-11){ 
              isReallyPhoton = kTRUE;
              FillHistogram("fHistMCCounter",20);
              AliAODMCParticle* V0Mother = NULL;
              if(V0Part->GetMother()!=-1) V0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part->GetMother())));
              if(V0Mother){if(V0Mother->GetPdgCode()==111){

                AliAODMCParticle* Pi0Mother = NULL;
                if(V0Mother->GetMother()!=-1) Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Mother->GetMother())));
                if(Pi0Mother){if(TMath::Abs(Pi0Mother->GetPdgCode())==3222){isPhotonfromSigma = kTRUE; FillHistogram("fHistMCCounter",21);}}

              }}//Mother of Photon exists and is a Pi0
            }}//Mother exists and is a Photon
          }//Both Tracks have a common Mother
        }//Both Tracks have matched MC Particle
      }//End of isMonteCarlo

      //First method: TLorentzVector
      trackElectron.SetXYZM(electrack->Px(),electrack->Py(),electrack->Pz(),cElectronMass);
      trackPositron.SetXYZM(positrack->Px(),positrack->Py(),positrack->Pz(),cElectronMass);
      trackPhoton = trackElectron + trackPositron;

      // Get Track parameters
      Double_t trackxyz[3];
      Double_t trackpxpypz[3];
      Double_t trackparams[6];
      Double_t covMatrix[21];

      // Set up KFParticle for the Electron
      electrack->GetXYZ(trackxyz);      
      electrack->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      electrack->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron.Create(trackparams,covMatrix,-1,cElectronMass);

      // Repeat for the Positron
      positrack->GetXYZ(trackxyz);      
      positrack->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      positrack->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron.Create(trackparams,covMatrix,1,cElectronMass); 

      // Combine Electron and Positron to KF Photon Particle
      KFParticle KFPhoton(KFElectron,KFPositron);

      Float_t mass,masserr;
      KFPhoton.GetMass(mass,masserr);

      KFPhoton.TransportToDecayVertex(); //Move Photon to the point of Conversion 
      Float_t VertexGamma[3] = {KFPhoton.GetX(), KFPhoton.GetY(), KFPhoton.GetZ()}; // Create Vertex and 
      KFElectron.TransportToPoint(VertexGamma);                                      // transport electron
      KFPositron.TransportToPoint(VertexGamma);                                      // and positron there  
      Float_t armenterosQtAlfa[2] = {0.};                                            
      KFParticle::GetArmenterosPodolanski(KFPositron,KFElectron,armenterosQtAlfa);

      FillHistogram("fHistKFeeDCA",TMath::Abs(KFElectron.GetDistanceFromParticle(KFPositron)));        
      FillHistogram("fHistArmPodKF",armenterosQtAlfa[1],armenterosQtAlfa[0]);            
      FillHistogram("fHisteeopenangle",KFElectron.GetAngle(KFPositron));

      if(isReallyPhoton){  
      FillHistogram("fHistKFeeDCAMC",TMath::Abs(KFElectron.GetDistanceFromParticle(KFPositron)));        
      FillHistogram("fHistArmPodKFMC",armenterosQtAlfa[1],armenterosQtAlfa[0]);            
      FillHistogram("fHisteeopenangleMC",KFElectron.GetAngle(KFPositron));          
      FillHistogram("fHisteephotmassKFMC",mass); //MC Photon Mass with KF
      } 
      FillHistogram("fHisteephotmassKFNoCut",mass); //Photon Mass with KF

      //Cuts
      if(TMath::Abs(armenterosQtAlfa[0])>fMaxPairqt) continue; //qt
      if(TMath::Abs(armenterosQtAlfa[1])>fMaxPairalpha) continue; //alpha
      if(KFElectron.GetAngle(KFPositron)>fMaxPairopenangle) continue; //openangle
      //if(TMath::Abs(KFElectron.GetDistanceFromParticle(KFPositron))>fMaxPairDCA) continue; //DCA

      FillHistogram("fHisteephotmassTL",trackPhoton.M()); //Photon Mass
      FillHistogram("fHisteephotmassKF",mass); //Photon Mass with KF

      // Inv. Mass Cut
      if(mass > fMaxphotonmass) continue;

      fElectronPairArray->push_back(fElectronArray->at(i));
      fPositronPairArray->push_back(fPositronArray->at(j));
      countGoodPairs++;

    }//End of Positron loop
  }//End of Electron loop

return;

} //End of PairElectrons()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::PairPhotonandElectron() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  Bool_t isReallyPi0 = kFALSE; 
  Bool_t isReallyPi0fromSigma = kFALSE;

  const Int_t nConvPhoton = fConvPhotonArray->size();
  const Int_t nElectron   = fElectronArray->size();
  const Int_t nPositron   = fPositronArray->size();
  const Int_t nProton     = fProtonArray->size();
  const Int_t nAntiProton = fAntiProtonArray->size();

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  KFParticle KFSingleElectron, KFElectron, KFPositron, KFProton; 
  TLorentzVector trackPhoton, trackElectron, trackPi0, trackProton, trackSigmaplus;  

  if(nConvPhoton==0) return;
  if(nElectron==0&&nPositron==0) return;

  for(Int_t i=0; i<nConvPhoton; i++) {

    AliAODv0 *photon = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(i));
    if(!photon) continue;

    for(Int_t j=0; j<nElectron+nPositron; j++) {

      //Check if electron is paired with second electron
      Bool_t ispaired=kFALSE;
      if(j<nElectron) for(UInt_t q=0; q<fElectronPairArray->size(); q++){if(fElectronArray->at(j)==fElectronPairArray->at(q)) ispaired=kTRUE;}
      else for(UInt_t q=0; q<fPositronPairArray->size(); q++){if(fPositronArray->at(j-nElectron)==fPositronPairArray->at(q)) ispaired=kTRUE;}
      if(ispaired) continue;

      AliAODTrack *electron;
      if(j<nElectron) electron = (AliAODTrack*)aodEvent->GetTrack(fElectronArray->at(j));
      else electron = (AliAODTrack*)aodEvent->GetTrack(fPositronArray->at(j-nElectron));
      if(!electron) continue;

      AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(photon->GetDaughter(0));
      AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(photon->GetDaughter(1));

      if(!track1 || !track2) {
        AliWarning("ERROR: Could not retrieve all AOD tracks in Pi0 reconstruction!");
        continue;
      }

      // AOD MC treatment
      Double_t MCPi0DCAPV; //MC Pi0 Decay Vertex; 
      isReallyPi0 = kFALSE; 
      isReallyPi0fromSigma = kFALSE;
      Int_t Pi0MotherLabel=-1;

      Double_t f3PartPi0DecayVert[3];
      Double_t f3PartPi0PhotConv[3];
      Double_t f3PartPi0ElecPos[3];

      if(isMonteCarlo){
        AliAODMCParticle* SingleElec = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(electron->GetLabel())));
        AliAODMCParticle* V0Daught1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track1->GetLabel())));
        AliAODMCParticle* V0Daught2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track2->GetLabel())));
        if(V0Daught1 && V0Daught2 && SingleElec){
          if(TMath::Abs(V0Daught1->GetPdgCode())==11&&TMath::Abs(V0Daught2->GetPdgCode())==11&&TMath::Abs(SingleElec->GetPdgCode())==11){
            if(V0Daught1->GetMother()!=-1&&SingleElec->GetMother()!=-1&&V0Daught1->GetMother()==V0Daught2->GetMother()){
              AliAODMCParticle* V0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Daught1->GetMother())));
              AliAODMCParticle* SingleElecMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(SingleElec->GetMother())));
              if(V0Part&&SingleElecMother){if(V0Part->GetPdgCode()==22&&SingleElecMother->GetPdgCode()==22){
                if(V0Part->GetMother()!=-1&&V0Part->GetMother()==SingleElecMother->GetMother()){

                  AliAODMCParticle* Pi0Part = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(V0Part->GetMother())));
                  if(Pi0Part){if(Pi0Part->GetPdgCode()==111){
                  
                    isReallyPi0=kTRUE; 

                    TVector3 prdvtx1(V0Part->Xv(),V0Part->Yv(),V0Part->Zv());   
                    TVector3 prdvtx2(SingleElecMother->Xv(),SingleElecMother->Yv(),SingleElecMother->Zv());   
                    TVector3 prdvtx=prdvtx1+prdvtx2; prdvtx*=0.5;
                    MCPi0DCAPV = prdvtx.Perp();
                    f3PartPi0DecayVert[0] = prdvtx.X();
                    f3PartPi0DecayVert[1] = prdvtx.Y();
                    f3PartPi0DecayVert[2] = prdvtx.Z();

                    TVector3 prdvtx3(V0Daught1->Xv(),V0Daught1->Yv(),V0Daught1->Zv());   
                    TVector3 prdvtx4(V0Daught2->Xv(),V0Daught2->Yv(),V0Daught2->Zv());   
                    TVector3 prdvtx5=prdvtx3+prdvtx4; prdvtx5*=0.5;
                    f3PartPi0PhotConv[0] = prdvtx5.X();
                    f3PartPi0PhotConv[1] = prdvtx5.Y();
                    f3PartPi0PhotConv[2] = prdvtx5.Z();

                    f3PartPi0ElecPos[0] = SingleElec->Xv();
                    f3PartPi0ElecPos[1] = SingleElec->Yv();
                    f3PartPi0ElecPos[2] = SingleElec->Zv();

                    AliAODMCParticle* Pi0Mother = NULL;
                    if(Pi0Part->GetMother()!=-1){
                      Pi0MotherLabel=Pi0Part->GetMother(); 
                      Pi0Mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(Pi0Part->GetMother())));
                      }
                    if(Pi0Mother){if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) isReallyPi0fromSigma = kTRUE;}

                  }}//Photon Mother exists and is Pi0
                }//Photon and Single electron have common mother
              }}//MC Photon exists
            }//Electrons have common mother
          }//Daughters are all electrons
        }//V0 Daughters and Single Electron exist
      }//End of isMonteCarlo

      trackPhoton.SetXYZM(photon->Px(),photon->Py(),photon->Pz(),0);
      trackElectron.SetXYZM(electron->Px(),electron->Py(),electron->Pz(),0);
      trackPi0 = trackPhoton + trackElectron;
      Double_t pi0mass = trackPi0.M();
      FillHistogram("fHist3Partpi0massvspt",pi0mass,electron->Pt());    
      if(electron->Pt()>0.25) FillHistogram("fHist3Partpi0mass",pi0mass);    
      if(isReallyPi0){
        FillHistogram("fHist3Partpi0massMC",pi0mass);    
        if(electron->Pt()>0.25) FillHistogram("fHist3Partpi0massMChighpt",pi0mass);    
      }

      // Set up KFParticle
      electron->GetXYZ(trackxyz);      
      electron->GetPxPyPz(trackpxpypz);
      for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
      for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
      electron->GetCovarianceXYZPxPyPz(covMatrix);
      KFSingleElectron.Create(trackparams,covMatrix,0,cElectronMass);

      // Repeat for all other particles
      trackparams[0] = photon->DecayVertexV0X();
      trackparams[1] = photon->DecayVertexV0Y();
      trackparams[2] = photon->DecayVertexV0Z();
      trackparams[3] = photon->MomPosX();
      trackparams[4] = photon->MomPosY();
      trackparams[5] = photon->MomPosZ();
      if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFElectron.Create(trackparams,covMatrix,-1,cElectronMass);

      trackparams[0] = photon->DecayVertexV0X();
      trackparams[1] = photon->DecayVertexV0Y();
      trackparams[2] = photon->DecayVertexV0Z();
      trackparams[3] = photon->MomPosX();
      trackparams[4] = photon->MomPosY();
      trackparams[5] = photon->MomPosZ();
      if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
      else track2->GetCovarianceXYZPxPyPz(covMatrix);
      KFPositron.Create(trackparams,covMatrix,1,cElectronMass);

      KFParticle KFPhoton(KFElectron,KFPositron);
      KFPhoton.TransportToDecayVertex();
      KFParticle KFPi0(KFPhoton,KFSingleElectron);
      KFPi0.TransportToDecayVertex();

      //Get the invariant mass to compare with TLorentzvector!
      Float_t mass,masserr;
      KFPi0.GetMass(mass,masserr);

      if(electron->Pt()>0.25) FillHistogram("fHist3Partpi0massKF",mass);
      //End of KF Pi0 calculation

      Double_t PhotElecDCA = TMath::Abs(KFPhoton.GetDistanceFromParticle(KFSingleElectron));

      Double_t KFPi0DCAPV = TMath::Sqrt((KFPi0.GetX()-primaryVtxPosX)*(KFPi0.GetX()-primaryVtxPosX)+(KFPi0.GetY()-primaryVtxPosY)*(KFPi0.GetY()-primaryVtxPosY));  
      if(isReallyPi0) FillHistogram("fHist3PartPi0VertexvsMC",KFPi0DCAPV-MCPi0DCAPV);

      if(isReallyPi0 && pi0mass>0 && pi0mass<0.3){ //Fill MC Pi0 Tree

        fInv3PartPi0Mass = pi0mass;
        f3PartPi0Px = trackPi0.Px(); 
        f3PartPi0Py = trackPi0.Py();
        f3PartPi0Pz = trackPi0.Pz();
        f3PartPi0DaughtPt1 = TMath::Sqrt(photon->MomPosX()*photon->MomPosX()+photon->MomPosY()*photon->MomPosY());
        f3PartPi0DaughtPt2 = TMath::Sqrt(photon->MomNegX()*photon->MomNegX()+photon->MomNegY()*photon->MomNegY());
        f3PartPi0DaughtPt3 = electron->Pt();
        f3PartPi0DecayVertX = f3PartPi0DecayVert[0];
        f3PartPi0DecayVertY = f3PartPi0DecayVert[1];
        f3PartPi0DecayVertZ = f3PartPi0DecayVert[2];
        f3PartPi0PhotPhotDCA = PhotElecDCA;
        f3PartPi0PhotConvx = f3PartPi0PhotConv[0];
        f3PartPi0PhotConvy = f3PartPi0PhotConv[1];
        f3PartPi0PhotConvz = f3PartPi0PhotConv[2];
        f3PartPrimVertX = primaryVtxPosX;
        f3PartPrimVertY = primaryVtxPosY;
        f3PartPrimVertZ = primaryVtxPosZ;
        f3PartPi0ElecPosx = f3PartPi0ElecPos[0];
        f3PartPi0ElecPosy = f3PartPi0ElecPos[1];
        f3PartPi0ElecPosz = f3PartPi0ElecPos[2];
        f3PartPi0CandTree->Fill();

      }

      if(pi0mass<0.05 || pi0mass>0.25) continue;  //Coarse mass cut to reduce combinatorics!

      /************************Sigma+ reconstruction************************************/

      for(Int_t k=0; k<nProton+nAntiProton; k++) {

        AliAODTrack *prot;
        if(k<nProton) prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray->at(k));
        else prot = (AliAODTrack*)aodEvent->GetTrack(fAntiProtonArray->at(k-nProton));
        if(!prot) continue;

        trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);
        trackSigmaplus = trackPi0 + trackProton;
        Double_t sigmaplusmass = trackSigmaplus.M();

        prot->GetXYZ(trackxyz);      
        prot->GetPxPyPz(trackpxpypz);
        for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
        for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
        prot->GetCovarianceXYZPxPyPz(covMatrix);
        if(k<nProton) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
        else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

        KFParticle KFSigmaPlus(KFProton,KFPi0);
        KFSigmaPlus.TransportToDecayVertex();

        Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));
        Float_t  DCAxy = -999., DCAz = -999.;
        prot->GetImpactParameters(DCAxy,DCAz);

        TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
        TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);

        Bool_t isReallySigma = kFALSE;
        if(isMonteCarlo && isReallyPi0fromSigma){
          AliAODMCParticle* ProtPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
          if(ProtPart){
            if(TMath::Abs(ProtPart->GetPdgCode())==2212&&ProtPart->GetMother()!=-1){
                if(ProtPart->GetMother()==Pi0MotherLabel) isReallySigma = kTRUE;    
            }//Proton has Mother
          }//Proton exists
        }//Is MonteCarlo and Pi0 is from Sigma        

        if(isReallySigma) FillHistogram("fHistMC4PartSigmaPA",TMath::Abs(sigmamomentum.Angle(sigmavertex)));

        // Fill the Sigma Candidate Trees
        f4PartIsMCSigma = kFALSE; if(isReallySigma) f4PartIsMCSigma = kTRUE;
        f4PartInvSigMass = sigmaplusmass; 
        f4PartSigCharge = prot->Charge(); 
        f4PartSigPx = trackSigmaplus.Px(); 
        f4PartSigPy = trackSigmaplus.Py(); 
        f4PartSigPz = trackSigmaplus.Pz(); 
        f4PartPrimVertX = primaryVtxPosX; 
        f4PartPrimVertY = primaryVtxPosY; 
        f4PartPrimVertZ = primaryVtxPosZ; 
        f4PartSigDecayVertX = KFSigmaPlus.GetX(); 
        f4PartSigDecayVertY = KFSigmaPlus.GetY(); 
        f4PartSigDecayVertZ = KFSigmaPlus.GetZ(); 
        f4PartInvPi0Mass = pi0mass; 
        f4PartPi0Px = trackPi0.Px(); 
        f4PartPi0Py = trackPi0.Py(); 
        f4PartPi0Pz = trackPi0.Pz(); 
        f4PartPi0DaughtPt1 = TMath::Sqrt(photon->MomPosX()*photon->MomPosX()+photon->MomPosY()*photon->MomPosY());
        f4PartPi0DaughtPt2 = TMath::Sqrt(photon->MomNegX()*photon->MomNegX()+photon->MomNegY()*photon->MomNegY());
        f4PartPi0DaughtPt3 = electron->Pt();
        f4PartPi0DecayVertX = KFPi0.GetX(); 
        f4PartPi0DecayVertY = KFPi0.GetY(); 
        f4PartPi0DecayVertZ = KFPi0.GetZ();
        f4PartPi0PhotPhotDCA = PhotElecDCA;
        f4PartProtonPx = prot->Px(); 
        f4PartProtonPy = prot->Py(); 
        f4PartProtonPz = prot->Pz(); 
        f4PartProtonDCAtoPVxy = DCAxy; 
        f4PartProtonDCAtoPVz = DCAz; 
        f4PartProtonPi0DCA = ProtPi0DCA; 
        if(f4PartInvSigMass<1.7) f4PartSigmaCandTree->Fill();

      }//End of Proton loop
    } //End of Electron Loop 
  } //End of Photon Loop

return;

} //End of PairPhotonandElectron()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticles() {

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  // Get Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  TLorentzVector trackPhoton1, trackPhoton2, trackPi0, trackProton, trackSigmaplus;
  KFParticle KFElectron1, KFElectron2, KFPositron1, KFPositron2, KFProton; 

  const Int_t nConvPhoton = fConvPhotonArray->size();
  const Int_t nProton = fProtonArray->size();
  const Int_t nAntiProton = fAntiProtonArray->size();
  Int_t countPi0 = 0;

  for(Int_t i=0; i<nConvPhoton-1; i++) {

    AliAODv0 *v0_1 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(i));
    if(!v0_1) continue;

    for(Int_t j=i+1; j<nConvPhoton; j++) {

      AliAODv0 *v0_2 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(j));
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
              if(V0Part1&&V0Part2){if(V0Part1->GetPdgCode()==22&&V0Part2->GetPdgCode()==22){
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
                      if(TMath::Abs(Pi0Mother->GetPdgCode())==3222) isReallyPi0fromSigma = kTRUE;
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
      FillHistogram("fHistpi0mass",pi0mass);
      if(isReallyPi0) FillHistogram("fHistMCpi0mass",pi0mass);

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

      //Reconstruct the Pi0 with the gammas
      KFParticle KFPi0(KFPhoton1,KFPhoton2);
      KFPi0.TransportToDecayVertex();
          
      //Get the invariant mass to compare with TLorentzvector!
      Float_t mass,masserr;
      KFPi0.GetMass(mass,masserr);

      FillHistogram("fHistpi0massKF",mass);
      //End of KF Pi0 calculation

      Double_t PhotPhotDCA = TMath::Abs(KFPhoton1.GetDistanceFromParticle(KFPhoton2));

      Double_t KFPi0DCAPV = TMath::Sqrt((KFPi0.GetX()-primaryVtxPosX)*(KFPi0.GetX()-primaryVtxPosX)+(KFPi0.GetY()-primaryVtxPosY)*(KFPi0.GetY()-primaryVtxPosY));  
      if(isReallyPi0){ 
        FillHistogram("fHistPi0VertexvsMC",KFPi0DCAPV-MCPi0DCAPV);
        FillHistogram("fHistPi0VertexMC",MCPi0DCAPV);
      }

      if(pi0mass>fMinPi0mass && pi0mass<fMaxPi0mass) countPi0++;

      if(pi0mass<0.05 || pi0mass>0.25) continue;  //Coarse mass cut to reduce combinatorics!

      /************************Sigma+ reconstruction************************************/

        for(Int_t k=0; k<nProton+nAntiProton; k++) {

          AliAODTrack *prot;
          if(k<nProton) prot = (AliAODTrack*)aodEvent->GetTrack(fProtonArray->at(k));
          else prot = (AliAODTrack*)aodEvent->GetTrack(fAntiProtonArray->at(k-nProton));
          if(!prot) continue;

          trackProton.SetXYZM(prot->Px(),prot->Py(),prot->Pz(),cProtonMass);
          trackSigmaplus = trackPi0 + trackProton;
          Double_t sigmaplusmass = trackSigmaplus.M();

          prot->GetXYZ(trackxyz);      
          prot->GetPxPyPz(trackpxpypz);
          for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
          for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
          prot->GetCovarianceXYZPxPyPz(covMatrix);
          if(k<nProton) KFProton.Create(trackparams,covMatrix,1,cProtonMass);
          else KFProton.Create(trackparams,covMatrix,-1,cProtonMass);

          KFParticle KFSigmaPlus(KFProton,KFPi0);
          KFSigmaPlus.TransportToDecayVertex();

          Double_t ProtPi0DCA = TMath::Abs(KFProton.GetDistanceFromParticle(KFPi0));
          Float_t  DCAxy = -999., DCAz = -999.;
          prot->GetImpactParameters(DCAxy,DCAz);

          TVector3 sigmamomentum(trackSigmaplus.Px(),trackSigmaplus.Py(),trackSigmaplus.Pz());
          TVector3 sigmavertex(KFSigmaPlus.GetX()-primaryVtxPosX,KFSigmaPlus.GetY()-primaryVtxPosY,KFSigmaPlus.GetZ()-primaryVtxPosZ);

          Bool_t isReallySigma = kFALSE;
          if(isMonteCarlo){
            AliAODMCParticle* ProtonPart = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(prot->GetLabel())));
            if(ProtonPart){if(TMath::Abs(ProtonPart->GetPdgCode())==2212){
              if(ProtonPart->GetMother()==Pi0MotherLabel){
                if(isReallyPi0fromSigma) {FillHistogram("fHistMCCounter",23); isReallySigma = kTRUE;}
                else if (isReallyPi0fromDelta) FillHistogram("fHistMCdeltamass",sigmaplusmass); 
              }//Proton and Pi0 have common Mother
            }}//MC Particle exists and is a Proton 
          }//End of isMonteCarlo

          if(isReallySigma) FillHistogram("fHistMCSigmaPA",TMath::Abs(sigmamomentum.Angle(sigmavertex)));

          // Fill the Sigma Candidate Trees
          fIsMCSigma = kFALSE; if(isReallySigma) fIsMCSigma = kTRUE;
          fInvSigMass = sigmaplusmass; 
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
          fInvPi0Mass = pi0mass; 
          fPi0Px = trackPi0.Px(); 
          fPi0Py = trackPi0.Py(); 
          fPi0Pz = trackPi0.Pz(); 
          fPi0DaughtPt1 = TMath::Sqrt(v0_1->MomPosX()*v0_1->MomPosX()+v0_1->MomPosY()*v0_1->MomPosY());
          fPi0DaughtPt2 = TMath::Sqrt(v0_1->MomNegX()*v0_1->MomNegX()+v0_1->MomNegY()*v0_1->MomNegY());
          fPi0DaughtPt3 = TMath::Sqrt(v0_2->MomPosX()*v0_2->MomPosX()+v0_2->MomPosY()*v0_2->MomPosY());
          fPi0DaughtPt4 = TMath::Sqrt(v0_2->MomNegX()*v0_2->MomNegX()+v0_2->MomNegY()*v0_2->MomNegY());
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
          if(fInvSigMass<1.7) fSigmaCandTree->Fill();

        }//End of Proton loop

      /************************End of Sigma+ reconstruction*****************************/        
        
      }//End of Photon 2 loop
    }//End of Photon 1 loop

  FillHistogram("fHistCandperEvent",countPi0,nProton); 

return;

} //End of ReconstructParticles()

//_____________________________________________________________________________

void AliAnalysisTaskSigmaPlus::ReconstructParticlesfromee() const{

  Double_t primaryVtxPos[3] = {primaryVtxPosX,primaryVtxPosY,primaryVtxPosZ};

  Double_t pi0mass, sigmaplusmass; //invariant masses

  // Track parameters
  Double_t trackxyz[3];
  Double_t trackpxpypz[3];
  Double_t trackparams[6];
  Double_t covMatrix[21];

  KFParticle KFElectron1, KFElectron2, KFPositron1, KFPositron2; 
  TLorentzVector trackElectron1, trackElectron2, trackPositron1, trackPositron2, trackPhoton1, trackPhoton2, trackPi0;

  const Int_t nElectrons = fElectronPairArray->size();
  const Int_t nProtons = fProtonArray->size();

  for(Int_t i=0; i<nElectrons-1; i++) {

    AliAODTrack *electrack1 = (AliAODTrack*)aodEvent->GetTrack(fElectronPairArray->at(i));
    AliAODTrack *positrack1 = (AliAODTrack*)aodEvent->GetTrack(fPositronPairArray->at(i));
    if(!electrack1||!positrack1) continue;

    for(Int_t j=i+1; j<nElectrons; j++) {

      AliAODTrack *electrack2 = (AliAODTrack*)aodEvent->GetTrack(fElectronPairArray->at(j));
      AliAODTrack *positrack2 = (AliAODTrack*)aodEvent->GetTrack(fPositronPairArray->at(j));
      if(!electrack2||!positrack2) continue;

      AliAODv0 *v0;
      Bool_t Phot1hasV0 = kFALSE; Int_t m_temp1 = 0;
      Bool_t Phot2hasV0 = kFALSE; Int_t m_temp2 = 0;
      Int_t mmax = fConvPhotonArray->size();
      for(Int_t m=0; m<mmax; m++){
        v0 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(m));
        if(!v0) continue;
        Int_t NID = v0->GetNegID();
        Int_t PID = v0->GetPosID();
        if(NID==electrack1->GetID()||PID==positrack1->GetID()){Phot1hasV0 = kTRUE; m_temp1 = m;}
        if(NID==electrack2->GetID()||PID==positrack2->GetID()){Phot2hasV0 = kTRUE; m_temp2 = m;}
      }

      if(Phot1hasV0&&Phot2hasV0) continue;

      if(Phot1hasV0){
        v0 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(m_temp1));        
        if(!v0) continue;
        AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
        AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
        if(!track1||!track2) continue;

        trackPhoton1.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),0);

        trackparams[0] = v0->DecayVertexV0X();
        trackparams[1] = v0->DecayVertexV0Y();
        trackparams[2] = v0->DecayVertexV0Z();
        trackparams[3] = v0->MomPosX();
        trackparams[4] = v0->MomPosY();
        trackparams[5] = v0->MomPosZ();
        if(track1->Charge()>0) track1->GetCovarianceXYZPxPyPz(covMatrix);
        else track2->GetCovarianceXYZPxPyPz(covMatrix);
        KFPositron1.Create(trackparams,covMatrix,1,cElectronMass);

        trackparams[0] = v0->DecayVertexV0X();
        trackparams[1] = v0->DecayVertexV0Y();
        trackparams[2] = v0->DecayVertexV0Z();
        trackparams[3] = v0->MomNegX();
        trackparams[4] = v0->MomNegY();
        trackparams[5] = v0->MomNegZ();
        if(track1->Charge()<0) track1->GetCovarianceXYZPxPyPz(covMatrix);
        else track2->GetCovarianceXYZPxPyPz(covMatrix);
        KFElectron1.Create(trackparams,covMatrix,-1,cElectronMass);
      }
      else{ //Else use KF
        // Set up KFParticle
        electrack1->GetXYZ(trackxyz);      
        electrack1->GetPxPyPz(trackpxpypz);
        for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
        for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
        electrack1->GetCovarianceXYZPxPyPz(covMatrix);
        KFElectron1.Create(trackparams,covMatrix,-1,cElectronMass);
        trackElectron1.SetXYZM(trackpxpypz[0],trackpxpypz[1],trackpxpypz[2],cElectronMass);

        // Repeat for all other particles
        positrack1->GetXYZ(trackxyz);      
        positrack1->GetPxPyPz(trackpxpypz);
        for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
        for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
        positrack1->GetCovarianceXYZPxPyPz(covMatrix);
        KFPositron1.Create(trackparams,covMatrix,1,cElectronMass);
        trackPositron1.SetXYZM(trackpxpypz[0],trackpxpypz[1],trackpxpypz[2],cElectronMass);

        //Reconstruct the Photon with KF
        KFParticle KFPhoton1(KFElectron1,KFPositron1);

        //Save KF Photon Parameters in TLorentzvectors
        trackPhoton1.SetXYZM(KFPhoton1.GetPx(),KFPhoton1.GetPy(),KFPhoton1.GetPz(),0);
      }

      if(Phot2hasV0){
        v0 = (AliAODv0*)aodEvent->GetV0(fConvPhotonArray->at(m_temp2));        
        if(!v0) continue;
        AliAODTrack* track3 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
        AliAODTrack* track4 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
        if(!track3||!track4) continue;

        trackPhoton1.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),0);

        trackparams[0] = v0->DecayVertexV0X();
        trackparams[1] = v0->DecayVertexV0Y();
        trackparams[2] = v0->DecayVertexV0Z();
        trackparams[3] = v0->MomPosX();
        trackparams[4] = v0->MomPosY();
        trackparams[5] = v0->MomPosZ();
        if(track3->Charge()>0) track3->GetCovarianceXYZPxPyPz(covMatrix);
        else track4->GetCovarianceXYZPxPyPz(covMatrix);
        KFPositron2.Create(trackparams,covMatrix,1,cElectronMass);

        trackparams[0] = v0->DecayVertexV0X();
        trackparams[1] = v0->DecayVertexV0Y();
        trackparams[2] = v0->DecayVertexV0Z();
        trackparams[3] = v0->MomNegX();
        trackparams[4] = v0->MomNegY();
        trackparams[5] = v0->MomNegZ();
        if(track3->Charge()<0) track3->GetCovarianceXYZPxPyPz(covMatrix);
        else track4->GetCovarianceXYZPxPyPz(covMatrix);
        KFElectron2.Create(trackparams,covMatrix,-1,cElectronMass);
      }
      else{ //Else use KF
        // Set up KFParticle
        electrack2->GetXYZ(trackxyz);      
        electrack2->GetPxPyPz(trackpxpypz);
        for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
        for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
        electrack2->GetCovarianceXYZPxPyPz(covMatrix);
        KFElectron2.Create(trackparams,covMatrix,-1,cElectronMass);
        trackElectron2.SetXYZM(trackpxpypz[0],trackpxpypz[1],trackpxpypz[2],cElectronMass);

        // Repeat for all other particles
        positrack2->GetXYZ(trackxyz);      
        positrack2->GetPxPyPz(trackpxpypz);
        for (Int_t q = 0; q<3;q++) {trackparams[q] = trackxyz[q];}
        for (Int_t q = 0; q<3;q++) {trackparams[q+3] = trackpxpypz[q];}
        positrack2->GetCovarianceXYZPxPyPz(covMatrix);
        KFPositron2.Create(trackparams,covMatrix,1,cElectronMass);
        trackPositron2.SetXYZM(trackpxpypz[0],trackpxpypz[1],trackpxpypz[2],cElectronMass);

        //Reconstruct the Photon with KF
        KFParticle KFPhoton2(KFElectron2,KFPositron2);

        //Save KF Photon Parameters in TLorentzvectors
        trackPhoton2.SetXYZM(KFPhoton2.GetPx(),KFPhoton2.GetPy(),KFPhoton2.GetPz(),0);
      }

      //Get the invariant mass with TLorentzVector!
      trackPi0 = trackPhoton1 + trackPhoton2;
      pi0mass = trackPi0.M();
      FillHistogram("fHisteepi0massTL",pi0mass);
        
    }//End of Photon 2 loop
  }//End of Photon 1 loop

return;

} //End of ReconstructParticlesfromee()

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
