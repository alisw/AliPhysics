/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
//#include "AliTriggerAnalysis.h"
#include "AliAODMCHeader.h"

// my headers
#include "AliAnalysisTaskUpcEtaC.h"

ClassImp(AliAnalysisTaskUpcEtaC);

using std::cout;
using std::endl;

//trees for UPC EtaC analysis,
// christopher.anson@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcEtaC::AliAnalysisTaskUpcEtaC() 
  : AliAnalysisTaskSE(),fType(0),fTracking(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fEtaCK0sChannelTree(0),fEtaCTree(0),fMeritCutChoice(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFmask(0),fIsPhysicsSelected(kFALSE),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fZDCAtime(0),fZDCCtime(0),fV0Adecision(0),fV0Cdecision(0),fADAdecision(0),fADCdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fEtaCAODTracks(0),fEtaCESDTracks(0),fGenPart(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fHistCTest58TriggersPerRun(0),fHistCTest59TriggersPerRun(0),fHistCTest60TriggersPerRun(0),fHistCTest61TriggersPerRun(0),fHistCcup8TriggersPerRun(0),fHistCcup9TriggersPerRun(0),fHistCcup10TriggersPerRun(0),fHistCcup11TriggersPerRun(0),fHistCcup12TriggersPerRun(0),fHistCtrueTriggersPerRun(0),
  fListHist(0),fHistNeventsEtaC(0),fMPiKvsMPiK(0),f2KstarPtPiPlus(0),f2KstarPtPiMinus(0),f2KstarPtKPlus(0),f2KstarPtKMinus(0),f2KstarTPCsignalPion(0),f2KstarTPCsignalKaon(0),f2KstarDedxVsPtPion(0),f2KstarDedxVsPtKaon(0),f2KstarTPCsignalVsQPtPion(0),f2KstarTPCsignalVsQPtKaon(0),f2KstarPtVsMinvFirstKstar(0),f2KstarPtVsMinvSecondKstar(0),f2KstarPtVsMinvEtaC(0),
  f1KstarPtPiPlus(0),f1KstarPtPiMinus(0),f1KstarPtKPlus(0),f1KstarPtKMinus(0),f1KstarTPCsignalPion(0),f1KstarTPCsignalKaon(0),f1KstarDedxVsPtPion(0),f1KstarDedxVsPtKaon(0),f1KstarTPCsignalVsQPtPion(0),f1KstarTPCsignalVsQPtKaon(0),f1KstarPtVsMinvKstar(0),f1KstarPtVsMinvOtherPiKcombo(0),f1KstarPtVsMinvEtaC(0),
  f0KstarPtPiPlus(0),f0KstarPtPiMinus(0),f0KstarPtKPlus(0),f0KstarPtKMinus(0),f0KstarTPCsignalPion(0),f0KstarTPCsignalKaon(0),f0KstarDedxVsPtPion(0),f0KstarDedxVsPtKaon(0),f0KstarTPCsignalVsQPtPion(0),f0KstarTPCsignalVsQPtKaon(0),f0KstarPtVsMinvFirstPiKcombo(0),f0KstarPtVsMinvSecondPiKcombo(0),f0KstarPtVsMinvEtaC(0),
  fHistK0sCandidatesPerEvent(0),fK0sPosDaughterPt(0),fK0sNegDaughterPt(0),fK0sPosVsNegDaughterPt(0),fPionK0sChannelPt(0),fKaonK0sChannelPt(0),fK0sPtVsMinv(0),fKPiPtVsMinvK0sChannel(0),fMK0sVsMKPiK0sChannel(0),fEtaCPtVsMinvK0sChannel(0),fK0sDecayLength(0),
  fHistNpion(0),fHistNK0sPion(0),fHistNkaon(0),fHistPiMinusK(0),
  fV0DaughterDca(0),fK0sDcaToPrimVertex(0),fK0sDaughterDcaToPrimVertex(0),fK0sMassDistribution(0),fV0DecayLength(0),fV0Eta(0),fCosPointingAngle(0),
    fHistNeventsEtaCK0sChannel(0),fHistEtaCMassVsPt(0),fHistEtaCMassCoherent(0),fHistZDCCuts(0),
  fNSigmaPionTPCvsNSigmaPionTOFLowPt(0),fNSigmaPionTPCvsNSigmaPionTOFMidPt(0),fNSigmaPionTPCvsNSigmaPionTOFHighPt(0),
  fNSigmaKaonTPCvsNSigmaKaonTOFLowPt(0),fNSigmaKaonTPCvsNSigmaKaonTOFMidPt(0),fNSigmaKaonTPCvsNSigmaKaonTOFHighPt(0),
  fTPCdEdxVsTOFbetaAll(0),fTPCdEdxVsTOFbetaPionsWithPID(0),fTPCdEdxVsTOFbetaKaonsWithPID(0),
  fTOFTimeVsTPCdEdxAll(0),fTOFTimeVsTPCdEdxPionsWithPID(0),fTOFTimeVsTPCdEdxKaonsWithPID(0),fNTracksWithTOFPIDPerEvent(0),fNTracksMissingDueToTOFPerEvent(0),
  fTOFbetaVsPtAll(0),fTOFbetaVsPtPionsWithPID(0),fTOFbetaVsPtKaonsWithPID(0),
    fListSystematics(0),fListJPsiLoose(0),fListJPsiTight(0),fListEtaCLoose(0),fListEtaCTight(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcEtaC


//_____________________________________________________________________________
AliAnalysisTaskUpcEtaC::AliAnalysisTaskUpcEtaC(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fTracking(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fEtaCK0sChannelTree(0),fEtaCTree(0),fMeritCutChoice(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFmask(0),fIsPhysicsSelected(kFALSE),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fZDCAtime(0),fZDCCtime(0),fV0Adecision(0),fV0Cdecision(0),fADAdecision(0),fADCdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fEtaCAODTracks(0),fEtaCESDTracks(0),fGenPart(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0), fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
  fHistCTest58TriggersPerRun(0),fHistCTest59TriggersPerRun(0),fHistCTest60TriggersPerRun(0),fHistCTest61TriggersPerRun(0),fHistCcup8TriggersPerRun(0),fHistCcup9TriggersPerRun(0),fHistCcup10TriggersPerRun(0),fHistCcup11TriggersPerRun(0),fHistCcup12TriggersPerRun(0),fHistCtrueTriggersPerRun(0),
  fListHist(0),fHistNeventsEtaC(0),fMPiKvsMPiK(0),f2KstarPtPiPlus(0),f2KstarPtPiMinus(0),f2KstarPtKPlus(0),f2KstarPtKMinus(0),f2KstarTPCsignalPion(0),f2KstarTPCsignalKaon(0),f2KstarTPCsignalVsQPtPion(0),f2KstarTPCsignalVsQPtKaon(0),f2KstarDedxVsPtPion(0),f2KstarDedxVsPtKaon(0),f2KstarPtVsMinvFirstKstar(0),f2KstarPtVsMinvSecondKstar(0),f2KstarPtVsMinvEtaC(0),
  f1KstarPtPiPlus(0),f1KstarPtPiMinus(0),f1KstarPtKPlus(0),f1KstarPtKMinus(0),f1KstarTPCsignalPion(0),f1KstarTPCsignalKaon(0),f1KstarDedxVsPtPion(0),f1KstarDedxVsPtKaon(0),f1KstarTPCsignalVsQPtPion(0),f1KstarTPCsignalVsQPtKaon(0),f1KstarPtVsMinvKstar(0),f1KstarPtVsMinvOtherPiKcombo(0),f1KstarPtVsMinvEtaC(0),
  f0KstarPtPiPlus(0),f0KstarPtPiMinus(0),f0KstarPtKPlus(0),f0KstarPtKMinus(0),f0KstarTPCsignalPion(0),f0KstarTPCsignalKaon(0),f0KstarDedxVsPtPion(0),f0KstarDedxVsPtKaon(0),f0KstarTPCsignalVsQPtPion(0),f0KstarTPCsignalVsQPtKaon(0),f0KstarPtVsMinvFirstPiKcombo(0),f0KstarPtVsMinvSecondPiKcombo(0),f0KstarPtVsMinvEtaC(0),
  fHistK0sCandidatesPerEvent(0),fK0sPosDaughterPt(0),fK0sNegDaughterPt(0),fK0sPosVsNegDaughterPt(0),fPionK0sChannelPt(0),fKaonK0sChannelPt(0),fK0sPtVsMinv(0),fKPiPtVsMinvK0sChannel(0),fMK0sVsMKPiK0sChannel(0),fEtaCPtVsMinvK0sChannel(0),fK0sDecayLength(0),
  fHistNpion(0),fHistNK0sPion(0),fHistNkaon(0),fHistPiMinusK(0),
  fV0DaughterDca(0),fK0sDcaToPrimVertex(0),fK0sDaughterDcaToPrimVertex(0),fK0sMassDistribution(0),fV0DecayLength(0),fV0Eta(0),fCosPointingAngle(0),
    fHistNeventsEtaCK0sChannel(0),fHistEtaCMassVsPt(0),fHistEtaCMassCoherent(0),fHistZDCCuts(0),
  fNSigmaPionTPCvsNSigmaPionTOFLowPt(0),fNSigmaPionTPCvsNSigmaPionTOFMidPt(0),fNSigmaPionTPCvsNSigmaPionTOFHighPt(0),
  fNSigmaKaonTPCvsNSigmaKaonTOFLowPt(0),fNSigmaKaonTPCvsNSigmaKaonTOFMidPt(0),fNSigmaKaonTPCvsNSigmaKaonTOFHighPt(0),
  fTPCdEdxVsTOFbetaAll(0),fTPCdEdxVsTOFbetaPionsWithPID(0),fTPCdEdxVsTOFbetaKaonsWithPID(0),
  fTOFTimeVsTPCdEdxAll(0),fTOFTimeVsTPCdEdxPionsWithPID(0),fTOFTimeVsTPCdEdxKaonsWithPID(0),fNTracksWithTOFPIDPerEvent(0),fNTracksMissingDueToTOFPerEvent(0),
  fTOFbetaVsPtAll(0),fTOFbetaVsPtPionsWithPID(0),fTOFbetaVsPtKaonsWithPID(0),
    fListSystematics(0),fListJPsiLoose(0),fListJPsiTight(0),fListEtaCLoose(0),fListEtaCTight(0)

{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;

  //cout << "##### fType = " << fType << " so AOD selected." << endl;

  fMeritCutChoice = 4;  //case for selecting best K0s candidate
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

}//AliAnalysisTaskUpcEtaC

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) {
  	fTrigger[i] = kFALSE;
	fTriggerInputsMC[i] = kFALSE;
	}
  for(Int_t i=0; i<4; i++) {
	fPIDTPCMuon[i] = -666;
	fPIDTPCElectron[i] = -666;
	fPIDTPCPion[i] = -666;
	fPIDTPCKaon[i] = -666;
	fPIDTPCProton[i] = -666;
	
	fPIDTOFMuon[i] = -666;
	fPIDTOFElectron[i] = -666;
	fPIDTOFPion[i] = -666;
	fPIDTOFKaon[i] = -666;
	fPIDTOFProton[i] = -666;
	
	fIsVtxContributor[i] = kFALSE;
	}
  for(Int_t i=0; i<4; i++) {
	fPIDTPCMuonPos[i] = -666;
	fPIDTPCElectronPos[i] = -666;
	fPIDTPCPionPos[i] = -666;
	fPIDTPCKaonPos[i] = -666;
	fPIDTPCProtonPos[i] = -666;
	
	fPIDTOFMuonPos[i] = -666;
	fPIDTOFElectronPos[i] = -666;
	fPIDTOFPionPos[i] = -666;
	fPIDTOFKaonPos[i] = -666;
	fPIDTOFProtonPos[i] = -666;

	fPIDTPCMuonNeg[i] = -666;
	fPIDTPCElectronNeg[i] = -666;
	fPIDTPCPionNeg[i] = -666;
	fPIDTPCKaonNeg[i] = -666;
	fPIDTPCProtonNeg[i] = -666;
	
	fPIDTOFMuonNeg[i] = -666;
	fPIDTOFElectronNeg[i] = -666;
	fPIDTOFPionNeg[i] = -666;
	fPIDTOFKaonNeg[i] = -666;
	fPIDTOFProtonNeg[i] = -666;
  }
  for(Int_t i=0; i<3; i++){
  	fVtxPos[i] = -666; 
	fMCVtxPos[i] = -666;
	fVtxErr[i] = -666;
	fKfVtxPos[i] = -666;
	fSpdVtxPos[i] = -666;
	}

  //cout << "##### end of Init()" << endl;
}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcEtaC::~AliAnalysisTaskUpcEtaC() 
{
  // Destructor
  if(fEtaCK0sChannelTree){
     delete fEtaCK0sChannelTree;
     fEtaCK0sChannelTree = 0x0;
  }
  if(fEtaCTree){
     delete fEtaCTree;
     fEtaCTree = 0x0;
  }
  if(fListTrig){
     delete fListTrig;
     fListTrig = 0x0;
  }
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
  }

}//~AliAnalysisTaskUpcEtaC


//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::UserCreateOutputObjects()
{
  //cout << "##### Start of UserCreateOutputObjects()" << endl;

  //PID response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  //input file
  fDataFilnam = new TObjString();
  fDataFilnam->SetString("");

    //tracks
  //DELETE fJPsiAODTracks = new TClonesArray("AliAODTrack", 1000);
  //DELETE  fJPsiESDTracks = new TClonesArray("AliESDtrack", 1000);
  fEtaCAODTracks = new TClonesArray("AliAODTrack", 1000);
  fEtaCESDTracks = new TClonesArray("AliESDtrack", 1000);
  fGenPart = new TClonesArray("TParticle", 1000);

  //output tree with JPsi candidate events
  //DELETE  fJPsiTree = new TTree("fJPsiTree", "fJPsiTree");
  //DELETE  fJPsiTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  //DELETE  fJPsiTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  //DELETE  fJPsiTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  //DELETE  fJPsiTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  //DELETE  fJPsiTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  //DELETE  fJPsiTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  //DELETE  fJPsiTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
   //DELETE fJPsiTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  //DELETE  fJPsiTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  //DELETE  fJPsiTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  //DELETE  fJPsiTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
  //DELETE  fJPsiTree ->Branch("fTOFmask", &fTOFmask);
  
  //DELETE  fJPsiTree ->Branch("fIsPhysicsSelected", &fIsPhysicsSelected, "fIsPhysicsSelected/O");
  
  //DELETE  fJPsiTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[2]/D");
  
  //DELETE  fJPsiTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[2]/D");
  //DELETE  fJPsiTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[2]/D");
  
  //DELETE  fJPsiTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[2]/O");
  
  //DELETE  fJPsiTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  //DELETE  fJPsiTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  //DELETE  fJPsiTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  //DELETE  fJPsiTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  //DELETE  fJPsiTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  //DELETE  fJPsiTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  //DELETE  fJPsiTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/D");
  //DELETE  fJPsiTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/D");
  //DELETE  fJPsiTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/D");
  //DELETE  fJPsiTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/D");
  //DELETE  fJPsiTree ->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/D");
  //DELETE  fJPsiTree ->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/D");
  //DELETE  fJPsiTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  //DELETE  fJPsiTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I"); 
  //DELETE  fJPsiTree ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  //DELETE  fJPsiTree ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");  
  //DELETE  fJPsiTree ->Branch("fDataFilnam", &fDataFilnam);
  //DELETE  fJPsiTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  //DELETE  fJPsiTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L"); 		       
  //DELETE  if( fType == 0 ) {
  //DELETE    fJPsiTree ->Branch("fJPsiESDTracks", &fJPsiESDTracks);
  //DELETE  }
  //DELETE  if( fType == 1 ) {
  //DELETE    fJPsiTree ->Branch("fJPsiAODTracks", &fJPsiAODTracks);
  //DELETE  }
  //DELETE  if(isMC) {
  //DELETE    fJPsiTree ->Branch("fGenPart", &fGenPart);
   //DELETE   fJPsiTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O", ntrg));
   //DELETE   fJPsiTree ->Branch("fMCVtxPos", &fMCVtxPos[0], "fMCVtxPos[3]/D");
   //DELETE }

 
 //output tree with EtaC candidate events
  fEtaCTree = new TTree("fEtaCTree", "fEtaCTree");
  fEtaCTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fEtaCTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fEtaCTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fEtaCTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fEtaCTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fEtaCTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fEtaCTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fEtaCTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fEtaCTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fEtaCTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fEtaCTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
  fEtaCTree ->Branch("fTOFmask", &fTOFmask);
  
  fEtaCTree ->Branch("fIsPhysicsSelected", &fIsPhysicsSelected, "fIsPhysicsSelected/O");
  
  fEtaCTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[4]/D");
  fEtaCTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[4]/D");
  fEtaCTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[4]/D");
  fEtaCTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[4]/D");
  fEtaCTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[4]/D");
  
  fEtaCTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[4]/D");
  fEtaCTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[4]/D");
  fEtaCTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[4]/D");
  fEtaCTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[4]/D");
  fEtaCTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[4]/D");
  
  fEtaCTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[4]/O");
  
  fEtaCTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fEtaCTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fEtaCTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fEtaCTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fEtaCTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  fEtaCTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fEtaCTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/D");
  fEtaCTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/D");
  fEtaCTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/D");
  fEtaCTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/D");
  fEtaCTree ->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/D");
  fEtaCTree ->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/D");
  fEtaCTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fEtaCTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I"); 
  fEtaCTree ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fEtaCTree ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");  
  fEtaCTree ->Branch("fDataFilnam", &fDataFilnam);
  fEtaCTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fEtaCTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fEtaCTree ->Branch("fEtaCESDTracks", &fEtaCESDTracks);
  }
  if( fType == 1 ) {
    fEtaCTree ->Branch("fEtaCAODTracks", &fEtaCAODTracks);
  }
  if(isMC) {
    fEtaCTree ->Branch("fGenPart", &fGenPart);
    fEtaCTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O", ntrg));
    fEtaCTree ->Branch("fMCVtxPos", &fMCVtxPos[0], "fMCVtxPos[3]/D");
  }

 //output tree with EtaC K0s Channel candidate events
  fEtaCK0sChannelTree = new TTree("fEtaCK0sChannelTree", "fEtaCK0sChannelTree");
  fEtaCK0sChannelTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fEtaCK0sChannelTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fEtaCK0sChannelTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fEtaCK0sChannelTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fEtaCK0sChannelTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fEtaCK0sChannelTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fEtaCK0sChannelTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fEtaCK0sChannelTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fEtaCK0sChannelTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fEtaCK0sChannelTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fEtaCK0sChannelTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
  fEtaCK0sChannelTree ->Branch("fTOFmask", &fTOFmask);
  
  fEtaCK0sChannelTree ->Branch("fIsPhysicsSelected", &fIsPhysicsSelected, "fIsPhysicsSelected/O");
  
  fEtaCK0sChannelTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[4]/D");
  
  fEtaCK0sChannelTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[4]/D");
  fEtaCK0sChannelTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[4]/D");
  
  fEtaCK0sChannelTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[4]/O");
  
  fEtaCK0sChannelTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fEtaCK0sChannelTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fEtaCK0sChannelTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fEtaCK0sChannelTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fEtaCK0sChannelTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  fEtaCK0sChannelTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fEtaCK0sChannelTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/D");
  fEtaCK0sChannelTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/D");
  fEtaCK0sChannelTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/D");
  fEtaCK0sChannelTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/D");
  fEtaCK0sChannelTree ->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/D");
  fEtaCK0sChannelTree ->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/D");
  fEtaCK0sChannelTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fEtaCK0sChannelTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I"); 
  fEtaCK0sChannelTree ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fEtaCK0sChannelTree ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");  
  fEtaCK0sChannelTree ->Branch("fDataFilnam", &fDataFilnam);
  fEtaCK0sChannelTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fEtaCK0sChannelTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fEtaCK0sChannelTree ->Branch("fEtaCESDTracks", &fEtaCESDTracks);
  }
  if( fType == 1 ) {
    fEtaCK0sChannelTree ->Branch("fEtaCAODTracks", &fEtaCAODTracks);
  }
  if(isMC) {
    fEtaCK0sChannelTree ->Branch("fGenPart", &fGenPart);
    fEtaCK0sChannelTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O", ntrg));
    fEtaCK0sChannelTree ->Branch("fMCVtxPos", &fMCVtxPos[0], "fMCVtxPos[3]/D");
  }
  //#####
  
  fListTrig = new TList();
  fListTrig ->SetOwner();
  
  fHistCcup4TriggersPerRun = new TH1D("fHistCcup4TriggersPerRun", "fHistCcup4TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup4TriggersPerRun);
  
  fHistCcup7TriggersPerRun = new TH1D("fHistCcup7TriggersPerRun", "fHistCcup7TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup7TriggersPerRun);
    
  fHistCcup2TriggersPerRun = new TH1D("fHistCcup2TriggersPerRun", "fHistCcup2TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup2TriggersPerRun);
  
  fHistCint1TriggersPerRun = new TH1D("fHistCint1TriggersPerRun", "fHistCint1TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCint1TriggersPerRun);
  
  fHistCint6TriggersPerRun = new TH1D("fHistCint6TriggersPerRun", "fHistCint6TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCint6TriggersPerRun);
  
  fHistC0tvxAndCint1TriggersPerRun = new TH1D("fHistC0tvxAndCint1TriggersPerRun", "fHistC0tvxAndCint1TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistC0tvxAndCint1TriggersPerRun);
  
  fHistZedTriggersPerRun = new TH1D("fHistZedTriggersPerRun", "fHistZedTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistZedTriggersPerRun);

  fHistCvlnTriggersPerRun = new TH1D("fHistCvlnTriggersPerRun", "fHistCvlnTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCvlnTriggersPerRun);
  
  fHistMBTriggersPerRun = new TH1D("fHistMBTriggersPerRun", "fHistMBTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistMBTriggersPerRun);
  
  fHistCentralTriggersPerRun = new TH1D("fHistCentralTriggersPerRun", "fHistCentralTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCentralTriggersPerRun);
  
  fHistSemiCentralTriggersPerRun = new TH1D("fHistSemiCentralTriggersPerRun", "fHistSemiCentralTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistSemiCentralTriggersPerRun);
  
  fHistCTest58TriggersPerRun = new TH1D("fHistCTest58TriggersPerRun", "fHistCTest58TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCTest58TriggersPerRun);
  
  fHistCTest59TriggersPerRun = new TH1D("fHistCTest59TriggersPerRun", "fHistCTest59TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCTest59TriggersPerRun);
  
  fHistCTest60TriggersPerRun = new TH1D("fHistCTest60TriggersPerRun", "fHistCTest60TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCTest60TriggersPerRun);
  
  fHistCTest61TriggersPerRun = new TH1D("fHistCTest61TriggersPerRun", "fHistCTest61TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCTest61TriggersPerRun);
  
  fHistCcup8TriggersPerRun = new TH1D("fHistCcup8TriggersPerRun", "fHistCcup8TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup8TriggersPerRun);
  
  fHistCcup9TriggersPerRun = new TH1D("fHistCcup9TriggersPerRun", "fHistCcup9TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup9TriggersPerRun);
  
  fHistCcup10TriggersPerRun = new TH1D("fHistCcup10TriggersPerRun", "fHistCcup10TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup10TriggersPerRun);
 
  fHistCcup11TriggersPerRun = new TH1D("fHistCcup11TriggersPerRun", "fHistCcup11TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup11TriggersPerRun);
  
  fHistCcup12TriggersPerRun = new TH1D("fHistCcup12TriggersPerRun", "fHistCcup12TriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCcup12TriggersPerRun);
  
  fHistCtrueTriggersPerRun = new TH1D("fHistCtrueTriggersPerRun", "fHistCtrueTriggersPerRun", 40000, 240000.5, 280000.5);
  fListTrig->Add(fHistCtrueTriggersPerRun);

  fListHist = new TList();
  fListHist ->SetOwner();

  cout << "##### Before fHistNevents... stuff is completed" << endl;
 
  TString CutNameEtaCK0sChannel[13] = {"Analyzed","Triggered","Vertex cut","V0 decision","Neutron ZDC cut","1 or more K0s candidates",
				       "Four good tracks","KPi and two daughter Pions","qK<0 && Sum(qPi)<0 (bad)","qK>0 && Sum(qPi)>0 (bad)",
				       "K-Pi+ && Opp sign daughters (good)","K+Pi- && Opp sign daughters","Total with correct charges"};

  fHistNeventsEtaCK0sChannel = new TH1D("fHistNeventsEtaCK0sChannel","fHistNeventsEtaCK0sChannel",13,0.5,13.5);
  for(Int_t i=0;i<13; i++) fHistNeventsEtaCK0sChannel->GetXaxis()->SetBinLabel(i+1,CutNameEtaCK0sChannel[i].Data());
  fListHist->Add(fHistNeventsEtaCK0sChannel);

  TString CutNameEtaC[15] = {"Analyzed","Triggered","Vertex cut","V0 decision","Neutron ZDC cut","Four good tracks",
			     "Two Kaons and Two Pions","Like sign kaons","Like sign pions","Like sign both","Opposite sign",
			     "Pion and Kaon mass check","2Kstar events","1Kstar events","0Kstar events"};

  fHistNeventsEtaC = new TH1D("fHistNeventsEtaC","fHistNeventsEtaC",15,0.5,15.5);
  for (Int_t i = 0; i<15; i++) fHistNeventsEtaC->GetXaxis()->SetBinLabel(i+1,CutNameEtaC[i].Data());
  fListHist->Add(fHistNeventsEtaC);

  //  cout << "##### After fHistNevents... stuff is completed" << endl;

  fMPiKvsMPiK = new TH2D("fMPiKvsMPiK","fMPiKvsMPiK",4000, 0.0, 40., 4000, 0.0, 40.);
  fListHist->Add(fMPiKvsMPiK);

  f2KstarPtPiPlus = new TH1D("f2KstarPtPiPlus","f2KstarPtPiPlus",300, 0., 3.);
  fListHist->Add(f2KstarPtPiPlus);
  f2KstarPtPiMinus = new TH1D("f2KstarPtPiMinus","f2KstarPtPiMinus",300, 0., 3.);
  fListHist->Add(f2KstarPtPiMinus);
  f2KstarPtKPlus = new TH1D("f2KstarPtKPlus","f2KstarPtKPlus",300, 0., 3.);
  fListHist->Add(f2KstarPtKPlus);
  f2KstarPtKMinus = new TH1D("2KstarPtKMinus","f2KstarPtKminus",300, 0., 3.);
  fListHist->Add(f2KstarPtKMinus);
  f2KstarTPCsignalPion = new TH2D("f2KstarTPCsignalPion","f2KstarTPCsignalPion",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f2KstarTPCsignalPion);
  f2KstarTPCsignalKaon = new TH2D("f2KstarTPCsignalKaon","f2KstarTPCsignalKaon",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f2KstarTPCsignalKaon);
  f2KstarDedxVsPtPion = new TH2D("f2KstarDedxVsPtPion","f2KstarDedxVsPtPion",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f2KstarDedxVsPtPion);
  f2KstarDedxVsPtKaon = new TH2D("f2KstarDedxVsPtKaon","f2KstarDedxVsPtKaon",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f2KstarDedxVsPtKaon);
  f2KstarPtVsMinvFirstKstar = new TH2D("f2KstarPtVsMinvFirstKstar","f2KstarPtVsMinvFirstKstar",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f2KstarPtVsMinvFirstKstar);
  f2KstarPtVsMinvSecondKstar = new TH2D("f2KstarPtVsMinvSecondKstar","f2KstarPtVsMinvSecondKstar",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f2KstarPtVsMinvSecondKstar);
  //f2KstarMinvFirstKstar = new TH1D("f2KstarMinvFirstKstar","f2KstarMinvFirstKstar",200, 0., 2.);
  //fListHist->Add(f2KstarMinvFirstKstar);
  //f2KstarMinvSecondKstar = new TH1D("f2KstarMinvSecondKstar","f2KstarMinvSecondKstar",200, 0., 2.);
  //fListHist->Add(f2KstarMinvSecondKstar);
  f2KstarPtVsMinvEtaC = new TH2D("f2KstarPtVsMinvEtaC","f2KstarPtVsMinvEtaC",1000, 0., 10.,1000, 0., 10.);
  fListHist->Add(f2KstarPtVsMinvEtaC);
  //f2KstarMinvEtaC = new TH1D("f2KstarMinvEtaC","f2KstarMinvEtaC",200, 2., 4.);
  //fListHist->Add(f2KstarMinvEtaC);

  f1KstarPtPiPlus = new TH1D("f1KstarPtPiPlus","f1KstarPtPiPlus",300,0.,3.);
  fListHist->Add(f1KstarPtPiPlus);
  f1KstarPtPiMinus = new TH1D("f1KstarPtPiMinus","f1KstarPtPiMinus",300,0.,3.);
  fListHist->Add(f1KstarPtPiMinus);
  f1KstarPtKPlus = new TH1D("f1KstarPtKPlus","f1KstarPtKPlus",300,0.,3.);
  fListHist->Add(f1KstarPtKPlus);
  f1KstarPtKMinus = new TH1D("f1KstarPtKMinus","f1KstarPtKMinus",300,0.,3.);
  fListHist->Add(f1KstarPtKMinus);
  f1KstarTPCsignalPion = new TH2D("f1KstarTPCsignalPion","f1KstarTPCsignalPion",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f1KstarTPCsignalPion);
  f1KstarTPCsignalKaon = new TH2D("f1KstarTPCsignalKaon","f1KstarTPCsignalKaon",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f1KstarTPCsignalKaon);
  f1KstarDedxVsPtPion = new TH2D("f1KstarDedxVsPtPion","f1KstarDedxVsPtPion",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f1KstarDedxVsPtPion);
  f1KstarDedxVsPtKaon = new TH2D("f1KstarDedxVsPtKaon","f1KstarDedxVsPtKaon",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f1KstarDedxVsPtKaon);
  f1KstarPtVsMinvKstar = new TH2D("f1KstarPtVsMinvKstar","f1KstarPtVsMinvKstar",1000, 0., 10.,1000, 0., 10.);
  fListHist->Add(f1KstarPtVsMinvKstar);
  f1KstarPtVsMinvOtherPiKcombo = new TH2D("f1KstarPtVsMinvOtherPiKcombo","f1KstarPtVsMinvOtherPiKcombo", 1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f1KstarPtVsMinvOtherPiKcombo);
  //f1KstarMinvKstar = new TH1D("f1KstarMinvKstar","f1KstarMinvKstar",200, 0., 2.);
  //fListHist->Add(f1KstarMinvKstar);
  //f1KstarMinvOtherPiKcombo = new TH1D("f1KstarMinvOtherPiKcombo","f1KstarMinvOtherPiKcombo",200, 0., 2.);
  //fListHist->Add(f1KstarMinvOtherPiKcombo);
  f1KstarPtVsMinvEtaC = new TH2D("f1KstarPtVsMinvEtaC","f1KstarPtVsMinvEtaC",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f1KstarPtVsMinvEtaC);
  //f1KstarMinvEtaC = new TH1D("f1KstarMinvEtaC","f1KstarMinvEtaC",200, 2., 4.);
  //fListHist->Add(f1KstarMinvEtaC);

  f0KstarPtPiPlus = new TH1D("f0KstarPtPiPlus","f0KstarPtPiPlus",300,0.,3.);  
  fListHist->Add(f0KstarPtPiPlus);
  f0KstarPtPiMinus = new TH1D("f0KstarPtPiMinus","fK0starPtPiMinus",300,0.,3.);
  fListHist->Add(f0KstarPtPiMinus);
  f0KstarPtKPlus = new TH1D("f0KstarPtKPlus","f0KstarPtKPlus",300,0.,3.);
  fListHist->Add(f0KstarPtKPlus);
  f0KstarPtKMinus = new TH1D("f0KstarPtKMinus","f0KstarPtKMinus",300,0.,3.);
  fListHist->Add(f0KstarPtKMinus);
  f0KstarTPCsignalPion = new TH2D("f0KstarTPCsignalPion","f0KstarTPCsignalPion",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f0KstarTPCsignalPion);
  f0KstarTPCsignalKaon = new TH2D("f0KstarTPCsignalKaon","f0KstarTPCsignalKaon",600, 0., 300., 600, 0., 300.);
  fListHist->Add(f0KstarTPCsignalKaon);
  f0KstarDedxVsPtPion = new TH2D("f0KstarDedxVsPtPion","f0KstarDedxVsPtPion",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f0KstarDedxVsPtPion);
  f0KstarDedxVsPtKaon = new TH2D("f0KstarDedxVsPtKaon","f0KstarDedxVsPtKaon",3000,0.,300., 3000,0.,300.);
  fListHist->Add(f0KstarDedxVsPtKaon);
  f0KstarPtVsMinvFirstPiKcombo = new TH2D("f0KstarPtVsMinvFirstPiKcombo","f0KstarPtVsMinvFirstPiKcombo",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f0KstarPtVsMinvFirstPiKcombo);
  f0KstarPtVsMinvSecondPiKcombo = new TH2D("f0KstarPtVsMinvSecondPiKcombo","f0KstarPtVsMinvSecondPiKcombo",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f0KstarPtVsMinvSecondPiKcombo);
  //f0KstarMinvFirstPiKcombo = new TH1D("f0KstarMinvFirstPiKcombo","f0KstarMinvFirstPiKcombo",200, 0., 2.);
  //fListHist->Add(f0KstarMinvFirstPiKcombo);
  //f0KstarPtVsMinvSecondPiKcombo = new TH1D("f0KstarMinvSecondPiKcombo","f0KstarMinvSecondPiKcombo",200, 0., 2.);
  //fListHist->Add(f0KstarMinvSecondPiKcombo);
  f0KstarPtVsMinvEtaC = new TH2D("f0KstarPtVsMinvEtaC","f0KstarPtVsMinvEtaC",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(f0KstarPtVsMinvEtaC);
  //f0KstarMinvEtaC = new TH1D("f0KstarMinvEtaC","f0KstarMinvEtaC",200, 2., 4.);
  //fListHist->Add(f0KstarMinvEtaC);

  //Some PID diagnostic plots
  fNSigmaPionTPCvsNSigmaPionTOFLowPt = new TH2D("fNSigmaPionTPCvsNSigmaPionTOFLowPt","fNSigmaPionTPCvsNSigmaPionTOFLowPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaPionTPCvsNSigmaPionTOFLowPt);
  fNSigmaPionTPCvsNSigmaPionTOFMidPt = new TH2D("fNSigmaPionTPCvsNSigmaPionTOFMidPt","fNSigmaPionTPCvsNSigmaPionTOFMidPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaPionTPCvsNSigmaPionTOFMidPt);
  fNSigmaPionTPCvsNSigmaPionTOFHighPt = new TH2D("fNSigmaPionTPCvsNSigmaPionTOFHighPt","fNSigmaPionTPCvsNSigmaPionTOFHighPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaPionTPCvsNSigmaPionTOFHighPt);
  fNSigmaKaonTPCvsNSigmaKaonTOFLowPt = new TH2D("fNSigmaKaonTPCvsNSigmaKaonTOFLowPt","fNSigmaKaonTPCvsNSigmaKaonTOFLowPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaKaonTPCvsNSigmaKaonTOFLowPt);
  fNSigmaKaonTPCvsNSigmaKaonTOFMidPt = new TH2D("fNSigmaKaonTPCvsNSigmaKaonTOFMidPt","fNSigmaKaonTPCvsNSigmaKaonTOFMidPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaKaonTPCvsNSigmaKaonTOFMidPt);
  fNSigmaKaonTPCvsNSigmaKaonTOFHighPt = new TH2D("fNSigmaKaonTPCvsNSigmaKaonTOFHighPt","fNSigmaKaonTPCvsNSigmaKaonTOFHighPt",900,-30,60,900,-30,60);
  fListHist->Add(fNSigmaKaonTPCvsNSigmaKaonTOFHighPt);
  fTPCdEdxVsTOFbetaAll = new TH2D("fTPCdEdxVsTOFbetaAll","fTPCdEdxVsTOFbetaAll",3000,0.,300.,200,0.,2.);
  fListHist->Add(fTPCdEdxVsTOFbetaAll);
  fTPCdEdxVsTOFbetaPionsWithPID = new TH2D("fTPCdEdxVsTOFbetaPionsWithPID","fTPCdEdxVsTOFbetaPionsWithPID",3000,0.,300.,200,0.,2.);
  fListHist->Add(fTPCdEdxVsTOFbetaPionsWithPID);
  fTPCdEdxVsTOFbetaKaonsWithPID = new TH2D("fTPCdEdxVsTOFbetaKaonsWithPID","fTPCdEdxVsTOFbetaKaonsWithPID",3000,0.,300.,200,0.,2.);
  fListHist->Add(fTPCdEdxVsTOFbetaKaonsWithPID);
  fTOFTimeVsTPCdEdxAll = new TH2D("fTOFTimeVsTPCdEdxAll","fTOFTimeVsTPCdEdxAll",3000,0.,300.,2000,-200.,200.);
  fListHist->Add(fTOFTimeVsTPCdEdxAll);
  fTOFTimeVsTPCdEdxPionsWithPID = new TH2D("fTOFTimeVsTPCdEdxPionsWithPID","fTOFTimeVsTPCdEdxPionsWithPID",3000,0.,300.,2000,-200.,200.);
  fListHist->Add(fTOFTimeVsTPCdEdxPionsWithPID);
  fTOFTimeVsTPCdEdxKaonsWithPID = new TH2D("fTOFTimeVsTPCdEdxKaonsWithPID","fTOFTimeVsTPCdEdxKaonsWithPID",3000,0.,300.,2000,-200.,200.);
  fListHist->Add(fTOFTimeVsTPCdEdxKaonsWithPID);
  fNTracksWithTOFPIDPerEvent = new TH1D("fNTracksWithTOFPIDPerEvent","fNTracksWithTOFPIDPerEvent",10,0.,10.);
  fListHist->Add(fNTracksWithTOFPIDPerEvent);
  fNTracksMissingDueToTOFPerEvent = new TH1D("fNTracksMissingDueToTOFPerEvent","fNTracksMissingDueToTOFPerEvent",10,0.,10.);
  fListHist->Add(fNTracksMissingDueToTOFPerEvent);
  fTOFbetaVsPtAll = new TH2D("fTOFbetaVsPtAll","fTOFbetaVsPtAll",1000,0.,10.,200,0.,2.);
  fListHist->Add(fTOFbetaVsPtAll);
  fTOFbetaVsPtPionsWithPID = new TH2D("fTOFbetaVsPtPionsWithPID","fTOFbetaVsPtPionsWithPID",1000,0.,10.,200,0.,2.);
  fListHist->Add(fTOFbetaVsPtPionsWithPID);
  fTOFbetaVsPtKaonsWithPID = new TH2D("fTOFbetaVsPtKaonsWithPID","fTOFbetaVsPtKaonsWithPID",1000,0.,10.,200,0.,2.);
  fListHist->Add(fTOFbetaVsPtKaonsWithPID);

  /*
  fHistEtaCMassVsPt = new TH2D("fHistEtaCMassVsPt","Mass vs p_{T} of #psi(2s) candidates",100,3,6,50,0,5);
  fHistEtaCMassVsPt->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fHistEtaCMassVsPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fListHist->Add(fHistEtaCMassVsPt);
  
  fHistEtaCMassCoherent = new TH1D("fHistEtaCMassAllCoherent","Invariant mass of coherent #psi(2s) candidates",50,2.5,5.5);
  fHistEtaCMassCoherent->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fListHist->Add(fHistEtaCMassCoherent);
  */

  fHistK0sCandidatesPerEvent = new TH1D("fHistK0sCandidatesPerEvent","fHistK0sCandidatesPerEvent",10,0.,10.);
  fListHist->Add(fHistK0sCandidatesPerEvent);
  fK0sPosDaughterPt = new TH1D("fK0sPosDaughterPt","fK0sPosDaughterPt",1000, 0., 10.);
  fListHist->Add(fK0sPosDaughterPt);
  fK0sNegDaughterPt = new TH1D("fK0sNegDaughterPt","fK0sNegDaughterPt",1000, 0., 10.);
  fListHist->Add(fK0sNegDaughterPt);
  fK0sPosVsNegDaughterPt = new TH2D("fK0sPosVsNegDaughterPt","fK0sPosVsNegDaughterPt",1000, 0., 10.,100, 0., 10.);
  fListHist->Add(fK0sPosVsNegDaughterPt);
  fPionK0sChannelPt = new TH1D("fPionK0sChannelPt","fPionK0sChannelPt",1000, 0., 10.);
  fListHist->Add(fPionK0sChannelPt);
  fKaonK0sChannelPt = new TH1D("fKaonK0sChannelPt","fKaonK0sChannelPt",1000, 0., 10.);
  fListHist->Add(fKaonK0sChannelPt);
  fK0sPtVsMinv = new TH2D("fK0sPtVsMinv","fK0sPtVsMinv",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(fK0sPtVsMinv);
  //fK0sMinv = new TH1D("fK0sMinv","fK0sMinv",200, 0.4, 0.6);
  //fListHist->Add(fK0sMinv);
  fKPiPtVsMinvK0sChannel = new TH2D("fKPiPtVsMinvK0sChannel","fKPiPtVsMinvK0sChannel",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(fKPiPtVsMinvK0sChannel);
  //fKPiMinvK0sChannel = new TH1D("fKPiMinvK0sChannel","fKPiMinvK0sChannel",500, 0.4, 0.9);
  //fListHist->Add(fKPiMinvK0sChannel);
  fMK0sVsMKPiK0sChannel = new TH2D("fMK0sVsMKPiK0sChannel","fMK0sVsMKPiK0sChannel",1000, 0., 1., 1000, 0., 1.);
  fListHist->Add(fMK0sVsMKPiK0sChannel);
  fEtaCPtVsMinvK0sChannel = new TH2D("fEtaCPtVsMinvK0sChannel","fEtaCPtVsMinvK0sChannel",1000, 0., 10., 1000, 0., 10.);
  fListHist->Add(fEtaCPtVsMinvK0sChannel);
  //fEtaCMinvK0sChannel = new TH1D("fEtaCMinvK0sChannel","fEtaCMinvK0sChannel",200, 2., 4.);
  //fListHist->Add(fEtaCMinvK0sChannel);
  fK0sDecayLength = new TH1D("fK0sDecayLength","fK0sDecayLength",300, 0., 30.);
  fListHist->Add(fK0sDecayLength);
  
  fHistNpion = new TH1D("fHistNpion","fHistNpion",10,0.,10.);
  fListHist->Add(fHistNpion);
  fHistNK0sPion = new TH1D("fHistNK0sPion","fHistNK0sPion",10,0.,10.);
  fListHist->Add(fHistNK0sPion);
  fHistNkaon = new TH1D("fHistNkaon","fHistNkaon",10,0.,10.);
  fListHist->Add(fHistNkaon);
  fHistPiMinusK = new TH1D("fHistPiMinusK","fHistPiMinusK",10,0.,10.);
  fListHist->Add(fHistPiMinusK);

  fV0DaughterDca = new TH1D("fV0DaughterDca","fV0DaughterDca",1000,0.,100.);
  fListHist->Add(fV0DaughterDca);
  fK0sDcaToPrimVertex = new TH1D("fK0sDcaToPrimVertex","fK0sDcaToPrimVertex",100,0.,10.);
  fListHist->Add(fK0sDcaToPrimVertex);
  fK0sDaughterDcaToPrimVertex = new TH1D("fK0sDaughterDcaToPrimVertex","fK0sDaughterDcaToPrimVertex",1000,0.,100.);
  fListHist->Add(fK0sDaughterDcaToPrimVertex);
  fK0sMassDistribution = new TH1D("fK0sMassDistribution","fK0sMassDistribution",600,0.2,0.8);
  fListHist->Add(fK0sMassDistribution);
  fV0DecayLength = new TH1D("fV0DecayLength","fV0DecayLength",200,0.,200.);
  fListHist->Add(fV0DecayLength);
  fV0Eta = new TH1D("fV0Eta","fV0Eta",160,-0.8,0.8);
  fListHist->Add(fV0Eta);
  fCosPointingAngle = new TH1D("fCosPointingAngle","fCosPointingAngle",500,0.,1.);
  fListHist->Add(fCosPointingAngle);
 
  TString CutNameZDC[4] = {"CCUP4","< 8 neutrons","0 netrons","No timing"};
  fHistZDCCuts = new TH1D("fHistZDCCuts","fHistZDCCuts",4,0.5,4.5);
  for (Int_t i = 0; i<4; i++) fHistZDCCuts->GetXaxis()->SetBinLabel(i+1,CutNameZDC[i].Data());
  fListHist->Add(fHistZDCCuts);
  
  
  fListSystematics = new TList();
  fListSystematics->SetOwner();
  fListSystematics->SetName("fListSystematics");
  fListHist->Add(fListSystematics);
  InitSystematics();

  cout << "##### End of UserCreateOutputObjects()" << endl;
  
  PostData(1, fEtaCK0sChannelTree);
  PostData(2, fEtaCTree);
  PostData(3, fListTrig);
  PostData(4, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::InitSystematics()
{ 

  //DELETEfListJPsiLoose = new TList();
  //DELETEfListJPsiLoose->SetOwner();
  //DELETEfListJPsiLoose->SetName("JPsiLoose");
  //DELETEfListSystematics->Add(fListJPsiLoose);

  //DELETETH1D *fHistJPsiNClusLoose = new TH1D("JPsiNClusLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiLoose->Add(fHistJPsiNClusLoose);

  //DELETETH1D *fHistJPsiChi2Loose = new TH1D("JPsiChi2Loose","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiLoose->Add(fHistJPsiChi2Loose);

  //DELETETH1D *fHistJPsiDCAzLoose = new TH1D("JPsiDCAzLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiLoose->Add(fHistJPsiDCAzLoose);

  //DELETETH1D *fHistJPsiDCAxyLoose = new TH1D("JPsiDCAxyLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiLoose->Add(fHistJPsiDCAxyLoose);

  //DELETETH1D *fHistJPsiITShitsLoose = new TH1D("JPsiITShitsLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiLoose->Add(fHistJPsiITShitsLoose);


  //DELETEfListJPsiTight = new TList();
  //DELETEfListJPsiTight->SetOwner();
  //DELETEfListJPsiTight->SetName("JPsiTight");
  //DELETEfListSystematics->Add(fListJPsiTight);

  //DELETETH1D *fHistJPsiNClusTight = new TH1D("JPsiNClusTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiTight->Add(fHistJPsiNClusTight);

  //DELETETH1D *fHistJPsiChi2Tight = new TH1D("JPsiChi2Tight","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiTight->Add(fHistJPsiChi2Tight);

  //DELETETH1D *fHistJPsiDCAzTight = new TH1D("JPsiDCAzTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiTight->Add(fHistJPsiDCAzTight);

  //DELETETH1D *fHistJPsiDCAxyTight = new TH1D("JPsiDCAxyTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
  //DELETEfListJPsiTight->Add(fHistJPsiDCAxyTight);


fListEtaCLoose = new TList();
fListEtaCLoose->SetOwner();
fListEtaCLoose->SetName("EtaCLoose");
fListSystematics->Add(fListEtaCLoose);

TH1D *fHistEtaCNClusLoose = new TH1D("EtaCNClusLoose","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCLoose->Add(fHistEtaCNClusLoose);

TH1D *fHistEtaCChi2Loose = new TH1D("EtaCChi2Loose","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCLoose->Add(fHistEtaCChi2Loose);

TH1D *fHistEtaCDCAzLoose = new TH1D("EtaCDCAzLoose","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCLoose->Add(fHistEtaCDCAzLoose);

TH1D *fHistEtaCDCAxyLoose = new TH1D("EtaCDCAxyLoose","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCLoose->Add(fHistEtaCDCAxyLoose);

TH1D *fHistEtaCITShitsLoose = new TH1D("EtaCITShitsLoose","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCLoose->Add(fHistEtaCITShitsLoose);


fListEtaCTight = new TList();
fListEtaCTight->SetOwner();
fListEtaCTight->SetName("EtaCTight");
fListSystematics->Add(fListEtaCTight);

TH1D *fHistEtaCNClusTight = new TH1D("EtaCNClusTight","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCTight->Add(fHistEtaCNClusTight);

TH1D *fHistEtaCChi2Tight = new TH1D("EtaCChi2Tight","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCTight->Add(fHistEtaCChi2Tight);

TH1D *fHistEtaCDCAzTight = new TH1D("EtaCDCAzTight","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCTight->Add(fHistEtaCDCAzTight);

TH1D *fHistEtaCDCAxyTight = new TH1D("EtaCDCAxyTight","Invariant mass of #psi(2S) candidates",50,2.5,5.5);
fListEtaCTight->Add(fHistEtaCDCAxyTight);


}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ){
    	RunESDtrig(); 
  	if(fRunHist) RunESDhist();
	if(fRunTree) RunESDtree();
	}

  if( fType == 1 ){
  	RunAODtrig(); 
  	if(fRunHist) RunAODhist();
	if(fRunTree) RunAODtree();
	}

}//UserExec
//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunAODtrig()
{
  //  cout << "########## Beginning of RunAODtrig()" << endl;

  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  fRunNum = aod ->GetRunNumber();
  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if(trigger.Contains("CCUP4-B")) fHistCcup4TriggersPerRun->Fill(fRunNum); //CCUP4 triggers
  if(trigger.Contains("CCUP7-B")) fHistCcup7TriggersPerRun->Fill(fRunNum); //CCUP7 triggers
  if(trigger.Contains("CCUP2-B")) fHistCcup2TriggersPerRun->Fill(fRunNum); //CCUP2 triggers
  
  if(trigger.Contains("CINT1-B")) fHistCint1TriggersPerRun->Fill(fRunNum); //CINT1 triggers
  
  if(trigger.Contains("CTEST58-B")) fHistCTest58TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST59-B")) fHistCTest59TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST60-B")) fHistCTest60TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST61-B")) fHistCTest61TriggersPerRun->Fill(fRunNum); //CTEST triggers
  
  if(trigger.Contains("CCUP8-B")) fHistCcup8TriggersPerRun->Fill(fRunNum); //CCUP8 triggers
  if(trigger.Contains("CCUP9-B")) fHistCcup9TriggersPerRun->Fill(fRunNum); //CCUP9 triggers
  if(trigger.Contains("CCUP10-B")) fHistCcup10TriggersPerRun->Fill(fRunNum); //CCUP10 triggers
  if(trigger.Contains("CCUP11-B")) fHistCcup11TriggersPerRun->Fill(fRunNum); //CCUP11 triggers
  if(trigger.Contains("CCUP12-B")) fHistCcup12TriggersPerRun->Fill(fRunNum); //CCUP12 triggers
  
  if(trigger.Contains("CTRUE-B")) fHistCtrueTriggersPerRun->Fill(fRunNum); //CTRUE triggers

  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  if(trigger.Contains("CINT1-B") && (fL0inputs & (1 << 3))) fHistC0tvxAndCint1TriggersPerRun->Fill(fRunNum); //0TVX triggers in CINT1 events
  
  if(trigger.Contains("CVLN_B2-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - synchronously downscaled
  if(trigger.Contains("CVLN_R1-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - randomly downscaled
  
  fL1inputs = aod->GetHeader()->GetL1TriggerInputs();
  if(fL1inputs & (1 << 18)) fHistZedTriggersPerRun->Fill(fRunNum); //1ZED trigger inputs
  
  //MB, Central and SemiCentral triggers
  AliCentrality *centrality = aod->GetCentrality();
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  Double_t percentile = centrality->GetCentralityPercentileUnchecked("V0M");
  //Double_t percentile = centrality->GetCentralityPercentile("V0M");
  
  if(((selectionMask & AliVEvent::kMB) == AliVEvent::kMB) && percentile<=80 && percentile>=0) fHistMBTriggersPerRun->Fill(fRunNum);
  
  if(((selectionMask & AliVEvent::kCentral) == AliVEvent::kCentral) && percentile<=6 && percentile>=0 && (trigger.Contains("CVHN_R2-B"))) fHistCentralTriggersPerRun->Fill(fRunNum);

  if(((selectionMask & AliVEvent::kSemiCentral) == AliVEvent::kSemiCentral) && percentile<=50 && percentile>=15) fHistSemiCentralTriggersPerRun->Fill(fRunNum);

  //  cout << "########## End of RunAODtrig()" << endl;

PostData(3, fListTrig);

}
//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunAODhist()
{

  //  cout << "##### Beginning of RunAODhist()" << endl;

  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Double_t kaonMass = partKaon->Mass();
    
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();

  TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
  Double_t kStarMass = partKstar->Mass();
  Double_t kStarWidth = partKstar->Width();

  TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
  Double_t k0ShortMass = partK0short->Mass();
  Double_t k0ShortWidth = partK0short->Width();

  //  cout << "mPi " << pionMass << ", mK " << kaonMass << ", mKstar " << kStarMass << ", wKstar " << kStarWidth << ", mk0s " << k0ShortMass << ", wk0s " << k0ShortWidth << endl;

  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;
  
 // cout<<"Event number: "<<((TTree*) GetInputData(0))->GetTree()->GetReadEntry()<<endl;

  fHistNeventsEtaCK0sChannel->Fill(1);
  fHistNeventsEtaC->Fill(1);

  //AliAODpidUtil
  //AliAODpidUtil *pidres = new AliAODpidUtil;
  //pidres->SetOADBPath("$ALICE_ROOT/OADB");

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if(!isMC && !trigger.Contains("CCUP") ) return;
  

  
  fHistNeventsEtaCK0sChannel->Fill(2);
  fHistNeventsEtaC->Fill(2);
  
  AliAODZDC *fZDCdata = aod->GetZDCData();
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  fZDCAtime = fZDCdata->GetZNATime();
  fZDCCtime = fZDCdata->GetZNCTime();
  
  if(trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1);
  if(fZNAenergy < 8200 && fZNCenergy < 8200) fHistZDCCuts->Fill(2);
  if(fZNAenergy < 683 && fZNCenergy < 683) fHistZDCCuts->Fill(3);
  if(fZDCAtime == 0 && fZDCCtime == 0) fHistZDCCuts->Fill(4);

  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsEtaCK0sChannel->Fill(3);
  fHistNeventsEtaC->Fill(3);

  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;
  
  fHistNeventsEtaCK0sChannel->Fill(4);
  fHistNeventsEtaC->Fill(4);

  if( fZNAenergy > 8200 || fZNCenergy > 8200) return;
  
  fHistNeventsEtaCK0sChannel->Fill(5);
  fHistNeventsEtaC->Fill(5); 
  
  //Systematics - cut variation
  if(fRunSystematics) RunAODsystematics(aod);

  //Two tracks loop
  Int_t nGoodTracks = 0;
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  Int_t missingTOFPID[5] = {-1,-1,-1,-1,-1};
  Int_t missingTOFPIDK0s[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vPion[4], vKaon[4], vK0sPion[4], vK0sKaon[4], vKPiK0sChannel, vK0s, vKstar[2], vCandidate;
  Short_t qKaon[4], qPion[4], qK0sPion[4], qK0sKaon[4];
  UInt_t nKaon=0, nPion=0, nK0sPion=0, nSpdHits=0;
  Double_t fRecTPCsignalPion[5], fRecTPCsignalKaon[5], fRecTPCsignalK0sPion[5];
  Int_t fChannel = 0;
  Double_t trackPt[5]={0,0,0,0,0};
  
   
  //Four track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));   
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = 4*(0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop

  //Diagnostic figures regarding nSigma TOF vs TPC, etc.  
  nKaon = 0; nPion = 0;
  //eventStartTime = pidres->GetTOFResponse().GetStartTime(trk->P());
  Float_t eventStartTime = 0.;
  Double_t beta = 0;
  Int_t nTracksWithTOFPID = 0;
  Int_t nLowPtTracks = 0;
  if(nGoodTracks == 4 && nSpdHits>1){
    for(Int_t i=0; i<4; i++){
      AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
      if(!trk) AliFatal("Not a standard AOD");

      //Get nsigma info for PID
      fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
      fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
      fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
      fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
      fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
      
      if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,trk) == AliPIDResponse::kDetPidOk) { //3 = kTOF
	fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
	fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
	fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
	fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
	fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
      } else {
	fPIDTOFMuon[i] = -999.;
	fPIDTOFElectron[i] = -999.;
	fPIDTOFPion[i] = -999.;
	fPIDTOFKaon[i] = -999.;
	fPIDTOFProton[i] = -999.;
      }

      //Check number of tracks with reasonable TOF PID available.
      if(trk->Pt() >= 0.6 && fPIDTOFPion[i] > -999.) nTracksWithTOFPID++;
      if(trk->Pt() < 0.6) nLowPtTracks++;

      //These are nSigmaTPC vs nSigmaTOF for all particles with Pion mass hypothesis and Kaon mass hypothesis.      
      if(trk->Pt() < 1.) {
            fNSigmaPionTPCvsNSigmaPionTOFLowPt->Fill(fPIDTOFPion[i],fPIDTPCPion[i]);
      	    fNSigmaKaonTPCvsNSigmaKaonTOFLowPt->Fill(fPIDTOFKaon[i],fPIDTPCKaon[i]);
      } else if(trk->Pt() >= 1. && trk->Pt() < 1.5) {
	    fNSigmaPionTPCvsNSigmaPionTOFMidPt->Fill(fPIDTOFPion[i],fPIDTPCPion[i]);
	    fNSigmaKaonTPCvsNSigmaKaonTOFMidPt->Fill(fPIDTOFKaon[i],fPIDTPCKaon[i]);
      } else if(trk->Pt() >= 1.5) {
      	    fNSigmaPionTPCvsNSigmaPionTOFHighPt->Fill(fPIDTOFPion[i],fPIDTPCPion[i]);
	    fNSigmaKaonTPCvsNSigmaKaonTOFHighPt->Fill(fPIDTOFKaon[i],fPIDTPCKaon[i]);
      }
      //These are TPC signal (dE/dx) vs TOF signal (time in ps) or TOF beta
      eventStartTime = aod->GetTOFHeader()->GetDefaultEventTimeVal(); //may have to define this differently. An array of start times is returned
      if(trk->GetIntegratedLength() >= 360. && trk->GetIntegratedLength() <= 800. && trk->GetTOFsignal() > 0. && eventStartTime < 999990.0) {

	beta = (trk->GetIntegratedLength()*0.01) / ( (trk->GetTOFsignal() - eventStartTime)*(1.e-12)*TMath::C() ); //0.01 converts cm to m, 1e-12 converts ps to sec

	fTPCdEdxVsTOFbetaAll->Fill(beta,trk->GetTPCsignal());
	fTOFbetaVsPtAll->Fill(trk->Pt(),beta);
	fTOFTimeVsTPCdEdxAll->Fill(trk->GetTPCsignal(),(trk->GetTOFsignal()-eventStartTime)); //time diff in ps. //(1e-12 converts ps to sec)
	if((trk->Pt() < 0.6 && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3.) || (trk->Pt() >=0.6 && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() > 3. && fabs(fPIDTPCPion[i]) < 2. && fabs(fPIDTPCKaon[i]) > 2.)) {
	  fTPCdEdxVsTOFbetaPionsWithPID->Fill(beta,trk->GetTPCsignal());
	  fTOFbetaVsPtPionsWithPID->Fill(trk->Pt(),beta);
	  fTOFTimeVsTPCdEdxPionsWithPID->Fill(trk->GetTPCsignal(),(trk->GetTOFsignal()-eventStartTime)); //time diff in ps. //(1e-12 converts ps to sec)
	}
	if((trk->Pt() < 0.6 && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTPCPion[i]) > 3.) || (trk->Pt() >=0.6 && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTOFPion[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTPCPion[i]) > 3.)) {
	  fTPCdEdxVsTOFbetaKaonsWithPID->Fill(beta,trk->GetTPCsignal());
	  fTOFbetaVsPtKaonsWithPID->Fill(trk->Pt(),beta);
	  fTOFTimeVsTPCdEdxKaonsWithPID->Fill(trk->GetTPCsignal(),(trk->GetTOFsignal()-eventStartTime)); //time diff in ps. //(1e-12 converts ps to sec)
	}
      }

    }
    fNTracksWithTOFPIDPerEvent->Fill(nTracksWithTOFPID);
    fNTracksMissingDueToTOFPerEvent->Fill(nGoodTracks-nLowPtTracks-nTracksWithTOFPID);
  }

  nKaon=0; nPion=0;
  Int_t nTracksWithoutTOFinfo = 0;

  if(nGoodTracks == 4 && nSpdHits>1){
  	  fHistNeventsEtaC->Fill(6);
  	  for(Int_t i=0; i<4; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
                if(!trk) AliFatal("Not a standard AOD");

		//		cout << "#################### Before PID block 1" << endl;

		//Get nsigma info for PID
		fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,trk) == AliPIDResponse::kDetPidOk) { //3 = kTOF
		  fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
		  fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
		  fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
		  fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
		  fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
		} else {
		  fPIDTOFMuon[i] = -999.;
		  fPIDTOFElectron[i] = -999.;
		  fPIDTOFPion[i] = -999.;
		  fPIDTOFKaon[i] = -999.;
		  fPIDTOFProton[i] = -999.;
		}

		//		cout << "#################### After PID block 1" << endl;

		//Here I need to identify Pions and Kaons. This block is for the 2pi2k final state
		//if(trk->GetMostProbablePID() == 2) { //Pions
		if((trk->Pt() < 0.6 && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >=3. && fabs(fPIDTPCPion[i]) < 2. && fabs(fPIDTPCKaon[i]) > 2.)) { 
		  fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
		  qPion[nPion] = trk->Charge();
		  vPion[nPion].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
		  nPion++;
		}
		//else if(trk->GetMostProbablePID() == 3) { //Kaons
		else if((trk->Pt() < 0.6 && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTPCPion[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTOFPion[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >= 3. && fabs(fPIDTPCKaon[i]) < 2. && fabs(fPIDTPCPion[i]) < 3.)) {
		  fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
		  qKaon[nKaon] = trk->Charge();
		  vKaon[nKaon].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),kaonMass);
		  nKaon++;
		} 
		else if(trk->Pt() >= 0.6 && trk->Pt() < 3. && fPIDTOFPion[i] == -999.) {
		  missingTOFPID[i] = trackIndex[i];
		  nTracksWithoutTOFinfo++; //Track index if TOF misses a track. Maybe due to missing TOF signal
		}

		if(nPion > 2 || nKaon > 2) break;
	  }
	    for(int i=0;i<4;i++) { //If one kaon or one pion is missing due to missing TOF PID info assume it is the fourth.
	      if(nPion == 2 && nKaon == 1 && missingTOFPID[i] > 0 && nTracksWithoutTOFinfo == 1) {
		AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(missingTOFPID[i]));
		qKaon[nKaon] = trk->Charge();
		vKaon[nKaon].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),kaonMass);
		nKaon++;
	      } else if(nKaon == 2 && nPion == 1 && missingTOFPID[i] > 0 && nTracksWithoutTOFinfo == 1) {
		AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(missingTOFPID[i]));
		qPion[nPion] = trk->Charge();
		vPion[nPion].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
		nPion++;
	      }
	    }
	    
	  //Analyze good events, fill the histos.
	  if( (nPion == 2) && (nKaon == 2) ) {
	    fHistNeventsEtaC->Fill(7);
	    if(qKaon[0]*qKaon[1] > 0) fHistNeventsEtaC->Fill(8);
	    if(qPion[0]*qPion[1] > 0) fHistNeventsEtaC->Fill(9);
	    if((qKaon[0]*qKaon[1] > 0) && (qPion[0]*qPion[1] > 0)) fHistNeventsEtaC->Fill(10);
	    if((qKaon[0]*qKaon[1] < 0) && (qPion[0]*qPion[1] < 0)) { 
	      fHistNeventsEtaC->Fill(11);
	      //TEMP if(vPion[0].M() == pionMass && vPion[1].M() == pionMass && vKaon[0].M() == kaonMass && vKaon[1].M() == kaonMass) {
		fHistNeventsEtaC->Fill(12);
		vCandidate = vPion[0] + vPion[1] + vKaon[0] + vKaon[1];
		//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

		  
		  //Get masses of potential intermediate Kstar(892)'s
		  if((qPion[0]*qKaon[0] < 0) && (qPion[1]*qKaon[1] < 0)) {
		    vKstar[0] = vPion[0] + vKaon[0];
		    vKstar[1] = vPion[1] + vKaon[1];
		  }
		  else if ((qPion[0]*qKaon[1] < 0) && (qPion[1]*qKaon[0] < 0)) {
		    vKstar[0] = vPion[0] + vKaon[1];
		    vKstar[1] = vPion[1] + vKaon[0];
		  }
		  
		  //Fill Dalitz plot with PiK masses Pi-K+ vs Pi+K-
		  if(qKaon[0] < 0) fMPiKvsMPiK->Fill(pow(vKstar[0].M(),2),pow(vKstar[1].M(),2));
		  else fMPiKvsMPiK->Fill(pow(vKstar[1].M(),2),pow(vKstar[0].M(),2));
		  
		  //Fill histos
		  //2 Kstar case
		  if((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
		     (vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth))) {
		    fHistNeventsEtaC->Fill(13);
		    if(qPion[0] > 0 && qPion[1] < 0) {
		      f2KstarPtPiPlus->Fill(vPion[0].Pt());
		      f2KstarPtPiMinus->Fill(vPion[1].Pt());
		    } else {
		      f2KstarPtPiPlus->Fill(vPion[1].Pt());
		      f2KstarPtPiMinus->Fill(vPion[0].Pt());
		    }
		    if(qKaon[0] > 0 && qKaon[1] < 0) {
		      f2KstarPtKPlus->Fill(vKaon[0].Pt());
		      f2KstarPtKMinus->Fill(vKaon[1].Pt());
		    } else {
		      f2KstarPtKPlus->Fill(vKaon[1].Pt());
		      f2KstarPtKMinus->Fill(vKaon[0].Pt());
		    }
		    //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		    f2KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		    f2KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		    f2KstarDedxVsPtPion->Fill(vPion[0].Pt(),fRecTPCsignalPion[0]);
		    f2KstarDedxVsPtPion->Fill(vPion[1].Pt(),fRecTPCsignalPion[1]);
		    f2KstarDedxVsPtKaon->Fill(vKaon[0].Pt(),fRecTPCsignalKaon[0]);
		    f2KstarDedxVsPtKaon->Fill(vKaon[1].Pt(),fRecTPCsignalKaon[1]);
		    //Fill intermediate Kstar histos, one for each candidate
		    f2KstarPtVsMinvFirstKstar->Fill(vKstar[0].M(),vKstar[0].Pt());
		    f2KstarPtVsMinvSecondKstar->Fill(vKstar[0].M(),vKstar[1].Pt());
		    //f2KstarMinvFirstKstar->Fill(vKstar[0].M());
		    //f2KstarMinvSecondKstar->Fill(vKstar[1].M());
		    //Fill EtaC histos
		    f2KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		    //f2KstarMinvEtaC->Fill(vCandidate.M());
		  }
		  //1 Kstar case
		  else if ((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
			   ((vKstar[1].M() > (kStarMass + kStarWidth)) || (vKstar[1].M() < (kStarMass - kStarWidth)))) {
		    //Fill using first Kstar candidate
		    fHistNeventsEtaC->Fill(14);
		    if(qPion[0] > 0 && qPion[1] < 0) {
		      f1KstarPtPiPlus->Fill(vPion[0].Pt());
		      f1KstarPtPiMinus->Fill(vPion[1].Pt());
		    } else {
		      f1KstarPtPiPlus->Fill(vPion[1].Pt());
		      f1KstarPtPiMinus->Fill(vPion[0].Pt());
		    }
		    if(qKaon[0] > 0 && qKaon[1] < 0) {
		      f1KstarPtKPlus->Fill(vKaon[0].Pt());
		      f1KstarPtKMinus->Fill(vKaon[1].Pt());
		    } else {
		      f1KstarPtKPlus->Fill(vKaon[1].Pt());
		      f1KstarPtKMinus->Fill(vKaon[0].Pt());
		    }
		    //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		    f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		    f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		    f1KstarDedxVsPtPion->Fill(vPion[0].Pt(),fRecTPCsignalPion[0]);
		    f1KstarDedxVsPtPion->Fill(vPion[1].Pt(),fRecTPCsignalPion[1]);
		    f1KstarDedxVsPtKaon->Fill(vKaon[0].Pt(),fRecTPCsignalKaon[0]);
		    f1KstarDedxVsPtKaon->Fill(vKaon[1].Pt(),fRecTPCsignalKaon[1]);
		    //Fill intermediate Kstar histos, one for each candidate
		    f1KstarPtVsMinvKstar->Fill(vKstar[0].M(),vKstar[0].Pt());
		    f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[1].M(),vKstar[1].Pt());
		    //f1KstarMinvKstar->Fill(vKstar[0].M());
		    //f1KstarMinvOtherPiKcombo->Fill(vKstar[1].M());
		    //Fill EtaC histos
		    f1KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		    //f1KstarMinvEtaC->Fill(vCandidate.M());
		  }
		  else if ((vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth)) &&
			   ((vKstar[0].M() > (kStarMass + kStarWidth)) || (vKstar[0].M() < (kStarMass - kStarWidth)))) {
		    //Fill using second Kstar candidate
		    fHistNeventsEtaC->Fill(14);
		    if(qPion[0] > 0 && qPion[1] < 0) {
		      f1KstarPtPiPlus->Fill(vPion[0].Pt());
		      f1KstarPtPiMinus->Fill(vPion[1].Pt());
		    } else {
		      f1KstarPtPiPlus->Fill(vPion[1].Pt());
		      f1KstarPtPiMinus->Fill(vPion[0].Pt());
		    }
		    if(qKaon[0] > 0 && qKaon[1] < 0) {
		      f1KstarPtKPlus->Fill(vKaon[0].Pt());
		      f1KstarPtKMinus->Fill(vKaon[1].Pt());
		    } else {
		      f1KstarPtKPlus->Fill(vKaon[1].Pt());
		      f1KstarPtKMinus->Fill(vKaon[0].Pt());
		    }
		    //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		    f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		    f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		    f1KstarDedxVsPtPion->Fill(vPion[0].Pt(),fRecTPCsignalPion[0]);
		    f1KstarDedxVsPtPion->Fill(vPion[1].Pt(),fRecTPCsignalPion[1]);
		    f1KstarDedxVsPtKaon->Fill(vKaon[0].Pt(),fRecTPCsignalKaon[0]);
		    f1KstarDedxVsPtKaon->Fill(vKaon[1].Pt(),fRecTPCsignalKaon[1]);
		    //Fill intermediate Kstar histos, one for each candidate
		    f1KstarPtVsMinvKstar->Fill(vKstar[1].M(),vKstar[1].Pt());
		    f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[0].M(),vKstar[0].Pt());
		    //f1KstarMinvKstar->Fill(vKstar[1].M());
		    //f1KstarMinvOtherPiKcombo->Fill(vKstar[0].M());
		    //Fill EtaC histos
		    f1KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		    //f1KstarMinvEtaC->Fill(vCandidate.M());
		  }
		  //0 Kstar case
		  else {
		    fHistNeventsEtaC->Fill(15);
		    if(qPion[0] > 0 && qPion[1] < 0) {
		      f0KstarPtPiPlus->Fill(vPion[0].Pt());
		      f0KstarPtPiMinus->Fill(vPion[1].Pt());
		    } else {
		      f0KstarPtPiPlus->Fill(vPion[1].Pt());
		      f0KstarPtPiMinus->Fill(vPion[0].Pt());
		    }
		    if(qKaon[0] > 0 && qKaon[1] < 0) {
		      f0KstarPtKPlus->Fill(vKaon[0].Pt());
		      f0KstarPtKMinus->Fill(vKaon[1].Pt());
		    } else {
		      f0KstarPtKPlus->Fill(vKaon[1].Pt());
		      f0KstarPtKMinus->Fill(vKaon[0].Pt());
		    }
		    //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		    f0KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		    f0KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		    f0KstarDedxVsPtPion->Fill(vPion[0].Pt(),fRecTPCsignalPion[0]);
		    f0KstarDedxVsPtPion->Fill(vPion[1].Pt(),fRecTPCsignalPion[1]);
		    f0KstarDedxVsPtKaon->Fill(vKaon[0].Pt(),fRecTPCsignalKaon[0]);
		    f0KstarDedxVsPtKaon->Fill(vKaon[1].Pt(),fRecTPCsignalKaon[1]);
		    //Fill intermediate Kstar histos, one for each candidate
		    f0KstarPtVsMinvFirstPiKcombo->Fill(vCandidate.M(),vKstar[0].Pt());
		    f0KstarPtVsMinvSecondPiKcombo->Fill(vCandidate.M(),vKstar[1].Pt());
		    //f0KstarMinvFirstPiKcombo->Fill(vKstar[0].M());
		    //f0KstarMinvSecondPiKcombo->Fill(vKstar[1].M());
		    //Fill EtaC histos
		    f0KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		    //f0KstarMinvEtaC->Fill(vCandidate.M());
		  }
		  //TEMP }
		}
	  }
  }


  //K0short case (using V0s)

  /*
  //Two track loop
  nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));   
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = 4*(0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
  */

  nGoodTracks = 0;
  
  //First identify best K0s candidate and it's daughter pions


  //Cuts and constants
  //const bool kMCCase = kFALSE;
  Bool_t fCutCheck = kFALSE; //change if cut checks are implimented.

  const float kMaxNumK0 = 300;
  const float kMinDCAPrimaryPion = 0.4; //Loosen or remove
  const float kMaxDCADaughtersK0 = 1.0; //originally was 0.3;
  const float kMaxDCAK0 = 0.3;
  const float kMaxDLK0 = 30.0;
  const float kMinDLK0 = 0.2; //fMinDecayLength;
  const float kEtaCut = 0.8;
  const float kMinCosAngle = 0.98; //originally was0.99;

  //const float kMinSeparation = fMinSep;
  const float kTOFLow = 0.8;
  const float kMaxTOFSigmaPion = 3.0;
  const float kMaxTPCSigmaPion = 3.0;

  //const float kMassPion = 0.13957;
  //const float kMassK0Short = 0.497614;

  //for cut checks
  double kCheckMassLow[3] = {0.49, 0.48, 0.45};
  double kCheckMassHigh[3] = {0.505, 0.515, 0.550};
  double kCheckDCAK0[3] = {0.1, 0.3, 1.0};
  double kCheckDCAPi[3] = {1.0, 0.4, 0.1};
  double kCheckDCAPiPi[3] = {0.1, 0.3, 1.0};
  double kCheckAvgSep[3] = {10.0, 5.0, 0.0};

  //v0 tester
  int v0Count = 0; //number of v0s (entries in array)
  int k0Count = 0; //number of good K0s

  ULong_t statusPos = 0;
  ULong_t statusNeg = 0;

  double newV0Pars[3] = { -666., -666., -666.};
  double oldV0Pars[3] = { -666., -666., -666.};

  int bestK0sCandidateIndex = 0;
  int k0ShortCount = 0;

  //skip AliFemtoK0Particle line
  //some MC stuff was commented out in femto code
  int UsedNMissingDaughters = 0;

  for(int i = 0; i < aod->GetNumberOfV0s(); i++) {
    bool goodK0 = kFALSE;
    bool goodPiPlus = kFALSE;
    bool goodPiMinus = kFALSE;
    int nMissingDaughters = 0;
    
    //PRINT    cout << "################################################## Num V0s " << aod->GetNumberOfV0s() << endl;

    if(aod->GetNumberOfV0s() > 10) break; 

    //load v0 track
    AliAODv0* v0 = aod->GetV0(i);
    if(!v0) continue;
    Bool_t fOnlineCase = kFALSE; //Change this code if implementing online case.
    if(fOnlineCase) { //need to define this somewhere (in input to analysis perhaps).
      if(!(v0->GetOnFlyStatus())) continue;
    } //for online
    else{
      if((v0->GetOnFlyStatus())) continue; //for offline
    }

    //for on-the-fly ordering
    AliAODTrack* tempTrack = (AliAODTrack*)v0->GetDaughter(0);
    short int pos0or1;
    short int neg0or1;
    bool orderswitch = kFALSE;
    if(tempTrack->Charge() > 0) { pos0or1 = 0; neg0or1 = 1;}
    else { pos0or1 = 1; neg0or1 = 0; orderswitch = kTRUE;}

    //load daughter tracks
    AliAODTrack* prongTrackPos = (AliAODTrack*)v0->GetDaughter(pos0or1);
    AliAODTrack* prongTrackNeg = (AliAODTrack*)v0->GetDaughter(neg0or1);
    if(!prongTrackPos) continue;
    if(!prongTrackNeg) continue;

    //daughter cuts
    if(v0->PtProng(pos0or1) < 0.15) continue;
    if(v0->PtProng(neg0or1) < 0.15) continue;
    if(fabs(v0->EtaProng(pos0or1)) > 0.8) continue;
    if(fabs(v0->EtaProng(neg0or1)) > 0.8) continue;

    //load status for PID
    statusPos=prongTrackPos->GetStatus();
    //TEMP COMMENT    if((statusPos&AliESDtrack::kTPCrefit)==0) continue;
    prongTrackPos->SetAODEvent(aod);
    statusNeg=prongTrackNeg->GetStatus();
    //TEMP COMMENT    if((statusNeg&AliESDtrack::kTPCrefit)==0) continue;
    prongTrackNeg->SetAODEvent(aod);

    //PRINT    cout << "#################### Before PID block 2" << endl;

    //TPC PID (fPIDAOD=fPIDResponse in this code) I may just use same selection code from above here.
    //if(fabs(fPidAOD->NumberOfSigmasTPC(prongTrackPos,AliPID::kPion)) < kMaxTPCSigmaPion) goodPiPlus = kTRUE;
    //if(fabs(fPidAOD->NumberOfSigmasTPC(prongTrackPos,AliPID::kPion)) < kMaxTPCSigmaPion) goodPiMinus = kTRUE;
    //Get nsigma info for PID
    fPIDTPCMuonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos,AliPID::kMuon);
    fPIDTPCElectronPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos,AliPID::kElectron);
    fPIDTPCPionPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos,AliPID::kPion);
    fPIDTPCKaonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos,AliPID::kKaon);
    fPIDTPCProtonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos,AliPID::kProton);
    
    if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackPos) == AliPIDResponse::kDetPidOk) { //3 = kTOF
      fPIDTOFMuonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos,AliPID::kMuon);
      fPIDTOFElectronPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos,AliPID::kElectron);
      fPIDTOFPionPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos,AliPID::kPion);
      fPIDTOFKaonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos,AliPID::kKaon);
      fPIDTOFProtonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos,AliPID::kProton);
    } else {
      fPIDTOFMuonPos[i] = -999.;
      fPIDTOFElectronPos[i] = -999.;
      fPIDTOFPionPos[i] = -999.;
      fPIDTOFKaonPos[i] = -999.;
      fPIDTOFProtonPos[i] = -999.;
    }

    fPIDTPCMuonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg,AliPID::kMuon);
    fPIDTPCElectronNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg,AliPID::kElectron);
    fPIDTPCPionNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg,AliPID::kPion);
    fPIDTPCKaonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg,AliPID::kKaon);
    fPIDTPCProtonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg,AliPID::kProton);
    
    if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackNeg) == AliPIDResponse::kDetPidOk) { //3 = kTOF
      fPIDTOFMuonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg,AliPID::kMuon);
      fPIDTOFElectronNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg,AliPID::kElectron);
      fPIDTOFPionNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg,AliPID::kPion);
      fPIDTOFKaonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg,AliPID::kKaon);
      fPIDTOFProtonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg,AliPID::kProton);
    } else {
      fPIDTOFMuonNeg[i] = -999.;
      fPIDTOFElectronNeg[i] = -999.;
      fPIDTOFPionNeg[i] = -999.;
      fPIDTOFKaonNeg[i] = -999.;
      fPIDTOFProtonNeg[i] = -999.;
    }

    //PRINT    cout << "#################### After PID block 2" << endl;

    //    if((fabs(fPIDTPCPionPos[i]) < 4.) && (fabs(fPIDTPCPionPos[i]) < fabs(fPIDTPCKaonPos[i]))) { 
      //More restrictive PID   fPIDTPCMuonPos[i] > 2 && fPIDTPCElectronPos[i] > 2 && fPIDTPCKaonPos[i] > 2 && fPIDTPCProtonPos[i] > 2) {
      //fRecTPCsignalPion[nPion] = prongTrackPos->GetTPCsignal();
      //qK0sPion[nPion] = prongTrackPos->Charge();
      //vK0sPion[nPion].SetPtEtaPhiM(prongTrackPos->Pt(),prongTrackPos->Eta(),prongTrackPos->Phi(),pionMass);
    if((prongTrackPos->Pt() < 0.6 && fabs(fPIDTPCPionPos[i]) < 3. && fabs(fPIDTPCKaonPos[i]) > 3.) || (prongTrackPos->Pt() >= 0.6 && prongTrackPos->Pt() < 3. && fabs(fPIDTOFPionPos[i]) < 3. && fabs(fPIDTOFKaonPos[i]) > 3. && fPIDTOFPionPos[i] > -999.) || (prongTrackPos->Pt() >= 3. && fabs(fPIDTPCPionPos[i]) < 2. && fabs(fPIDTPCKaonPos[i]) > 2.)) {
      goodPiPlus = kTRUE;
      //nPion++;
    } else if(prongTrackPos->Pt() >= 0.6 && prongTrackPos->Pt() < 3. && fPIDTOFPionPos[i] == -999. && nMissingDaughters == 0) {
      goodPiPlus = kTRUE; //This may be redundant. Maybe just assume the daughter is a pion. Or if prongTrackNeg has PID as a pion assume this is a positive pion.
      nMissingDaughters++;
    }
    //    if((fabs(fPIDTPCPionNeg[i]) < 4.) && (fabs(fPIDTPCPionNeg[i]) < fabs(fPIDTPCKaonNeg[i]))) { 
      //More restrictive PID   fPIDTPCMuonNeg[i] > 2 && fPIDTPCElectronNeg[i] > 2 && fPIDTPCKaonNeg[i] > 2 && fPIDTPCProtonNeg[i] > 2) {
      //fRecTPCsignalPion[nPion] = prongTrackNeg->GetTPCsignal();
      //qK0sPion[nPion] = prongTrackNeg->Charge();
      //vK0sPion[nPion].SetPtEtaPhiM(prongTrackNeg->Pt(),prongTrackNeg->Eta(),prongTrackNeg->Phi(),pionMass);
    if((prongTrackNeg->Pt() < 0.6 && fabs(fPIDTPCPionNeg[i]) < 3. && fabs(fPIDTPCKaonNeg[i]) > 3.) || (prongTrackNeg->Pt() >= 0.6 && prongTrackNeg->Pt() < 3. && fabs(fPIDTOFPionNeg[i]) < 3. && fabs(fPIDTOFKaonNeg[i]) > 3. && fPIDTOFKaonNeg[i] > -999.) || (prongTrackNeg->Pt() >= 3. && fabs(fPIDTPCPionNeg[i]) < 2. && fabs(fPIDTPCKaonNeg[i]) > 2.)) {
      goodPiMinus = kTRUE;
      //nPion++;
      } else if(prongTrackNeg->Pt() >= 0.6 && prongTrackNeg->Pt() < 3. && fPIDTOFPionNeg[i] == -999. && nMissingDaughters == 0) {
	goodPiMinus = kTRUE;
	nMissingDaughters++;
      }
    //We now have selected a V0 that has decayed to two pions

    //PRINT    cout << "A" << endl;

    //Skip PID using TOF for now

    //K0s cuts
    //PRINT    cout << "V0 mass " << v0->MassK0Short() << ", Decay Length " << v0->DecayLength(fAODVertex) << ", Eta " << v0->Eta() << ", CosPA " << v0->CosPointingAngle(fAODVertex) << endl;
    fV0DecayLength->Fill(v0->DecayLength(fAODVertex));
    fV0Eta->Fill(v0->Eta());
    fCosPointingAngle->Fill(v0->CosPointingAngle(fAODVertex));
    if(!goodPiMinus || !goodPiPlus) continue;
    if(v0->Eta() > kEtaCut) continue;
    if(v0->CosPointingAngle(fAODVertex) < kMinCosAngle) continue;
    if(v0->MassK0Short() < 0.2 || v0->MassK0Short() > 0.8) continue;
    if(v0->DecayLength(fAODVertex) > kMaxDLK0) continue;
    if(v0->DecayLength(fAODVertex) < kMinDLK0) continue;

    //PRINT    cout << "A 2" << endl;

    double v0Dca = v0->DcaV0ToPrimVertex();
    fK0sDcaToPrimVertex->Fill(v0Dca);
    fK0sDaughterDcaToPrimVertex->Fill(v0->DcaNegToPrimVertex());
    fK0sDaughterDcaToPrimVertex->Fill(v0->DcaPosToPrimVertex());
    fV0DaughterDca->Fill(v0->DcaV0Daughters());
    fK0sMassDistribution->Fill(v0->MassK0Short());
    if(!fCutCheck) {
      //PRINT      cout << "A 3" << endl;
      if(v0->DcaNegToPrimVertex() < kMinDCAPrimaryPion) continue;
      if(v0->DcaPosToPrimVertex() < kMinDCAPrimaryPion) continue;
      if(v0->DcaV0Daughters() > kMaxDCADaughtersK0) continue;
      if(v0Dca > kMaxDCAK0) continue;
      //PRINT      cout << "A 4" << endl;
    } else {
      //PRINT      cout << "A 5" << endl;
      if(v0->DcaNegToPrimVertex() < kCheckDCAPi[2]) continue;
      if(v0->DcaPosToPrimVertex() < kCheckDCAPi[2]) continue;
      if(v0->DcaV0Daughters() > kCheckDCAPiPi[2]) continue;
      if(v0Dca > kCheckDCAK0[2]) continue;
      //PRINT      cout << "A 6" << endl;
    }

    //EVERYTHING BELOW HERE PASSES SINGLE PARTICLE CUTS, PION PID, and LOOSE MASS CUT

    //Skipping some MC stuff that is commented out

    if(!fCutCheck) {
      //PRINT      cout << "A 7" << endl;
      if(v0->MassK0Short() > 0.48 && v0->MassK0Short() < 0.515) goodK0 = kTRUE;
    } else {
      //PRINT      cout << "A 8" << endl;
      if(v0->MassK0Short() > kCheckMassLow[2] && v0->MassK0Short() < kCheckMassHigh[2]) goodK0 = kTRUE;
    }

    //PRINT    cout << "B" << endl;

    //Now we need to check for the best K0s. The femto code eliminates cases with shared daughters because
    //that would add fake correlation signal. We will have one or two candidates only (unless there are v0s
    //that might be from pile up. We shall require (1) closest mass to K0s mass, (2) DCA K0s closest to 
    //primary vertex, and (3) DCA daughters closest to each other (at point of decay).
    //This loop is over the number of V0s in the event. I will add a histo fNumberV0sPerEvent in events 
    //expecting it to be low and another with fNumber K0sPerEvent.
    //We will only select the best K0s according to the criteria above which should be very efficient.
    //Will skip shared daughters check. I just want the best one.

    //Check for shared daughters, using v0 DCA, etc.
    bool v0JudgeNew; //true if new v0 beats old
    //tempK0[v0Count].fSkipShared = kFALSE;
    newV0Pars[0] = fabs(v0->MassK0Short()-k0ShortMass);
    newV0Pars[1] = v0Dca;
    newV0Pars[2] = v0->DcaV0Daughters(); //parameters used in merit cut

    //PRINT    cout << "B 1" << endl;

    if(i == 0) {
      //PRINT      cout << "B 2" << endl;
      oldV0Pars[0] = fabs(v0->MassK0Short()-k0ShortMass);
      oldV0Pars[1] = v0Dca;
      oldV0Pars[2] = v0->DcaV0Daughters(); //first time through, skip comparison
      bestK0sCandidateIndex = i;
      UsedNMissingDaughters = nMissingDaughters;
      //PRINT      cout << "B 3" << endl;
    } else {
      //PRINT      cout << "B 4" << endl;
      v0JudgeNew = CheckMeritCutWinner(fMeritCutChoice, oldV0Pars, newV0Pars); //true if new wins
      if(v0JudgeNew) {
	//PRINT	cout << "B 5" << endl;
	oldV0Pars[0] = fabs(v0->MassK0Short()-k0ShortMass);
	oldV0Pars[1] = v0Dca;
	oldV0Pars[2] = v0->DcaV0Daughters(); //new beats old, set oldV0Par values for next pass
       	bestK0sCandidateIndex = i; //update index of best K0s candidate
	UsedNMissingDaughters = nMissingDaughters; //Update number of missing daughters for the current, best K0s candidate
	//PRINT	cout << "B 6" << endl;
      }
    }

    //PRINT    cout << "c " << k0ShortCount << endl;

    k0ShortCount++;
  }	//Just found the index for the best k0s candidate
  if(k0ShortCount > 0) fHistK0sCandidatesPerEvent->Fill(k0ShortCount); //Diagnostic histo.
  
  if(k0ShortCount > 0) {
    fHistNeventsEtaCK0sChannel->Fill(6);

    //Now create data structures for Pions from best K0s candidate.
    AliAODv0* v0 = aod->GetV0(bestK0sCandidateIndex);
    AliAODTrack* tempTrack = (AliAODTrack*)v0->GetDaughter(0);
    short int pos0or1;
    short int neg0or1;
    bool orderswitch = kFALSE;
    if(tempTrack->Charge() > 0) { pos0or1 = 0; neg0or1 = 1;}
    else { pos0or1 = 1; neg0or1 = 0; orderswitch = kTRUE;}
    AliAODTrack* prongTrackPos = (AliAODTrack*)v0->GetDaughter(pos0or1);
    AliAODTrack* prongTrackNeg = (AliAODTrack*)v0->GetDaughter(neg0or1);
    
    //PRINT    cout << "D" << endl;
    
    //Four track loop
    nPion = 0; nK0sPion = 0; nKaon = 0;
    Int_t currentID = -666;
    Int_t posTrackIndex = -666;
    Int_t negTrackIndex = -666;
    Bool_t skipTrack = kFALSE;
    Int_t nProngFound = 0;

    for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
      AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));   
      if( !trk ) continue;
      if(!(trk->TestFilterBit(1<<0))) continue;
      
      currentID = trk->GetID();
      skipTrack=kFALSE;
      if(currentID == prongTrackPos->GetID()) { posTrackIndex = itr; skipTrack = kTRUE; nProngFound++; }
      if(currentID == prongTrackNeg->GetID()) { negTrackIndex = itr; skipTrack = kTRUE; nProngFound++; }

      if(!skipTrack) {      
	if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
	//if(!skipTrack) { 
	if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      }
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!skipTrack) {
	if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
	delete trk_clone;
      }
      if(!skipTrack) { //Only apply DCA and ITS hit requirements on primary tracks, not on K0s daughters
        if(TMath::Abs(dca[1]) > 2) continue;
        Double_t cut_DCAxy = 4*(0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
        if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
        if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
      }
      
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
      
      if(nGoodTracks > 4) break;  
    }//Track loop
    
    //PRINT    cout << "E" << endl;

    //now identify the additional kaon and pion, excluding the pions from K0s daughters by fID.
    nKaon=0; nPion=0; nK0sPion=0;
    qK0sPion[0] = 0;
    qK0sPion[1] = 0;
    qPion[0] = 0;
    qKaon[0] = 0;
    Int_t nTracksWithoutTOFinfoK0s = 0;

    if(nGoodTracks == 4 && nSpdHits>1 && nProngFound == 2){
      fHistNeventsEtaCK0sChannel->Fill(7);
      for(Int_t i=0; i<4; i++){
	AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
	if(!trk) AliFatal("Not a standard AOD");
	
	//PRINT	cout << "#################### Before PID block 3" << endl;

	//Get nsigma info for PID
	fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
	fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
	fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
	fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
	fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
	
    if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,trk) == AliPIDResponse::kDetPidOk) { //3 = kTOF
	fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
	fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
	fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
	fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
	fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
    } else {
      fPIDTOFMuonPos[i] = -999.;
      fPIDTOFElectronPos[i] = -999.;
      fPIDTOFPionPos[i] = -999.;
      fPIDTOFKaonPos[i] = -999.;
      fPIDTOFProtonPos[i] = -999.;
    }
	
	//PRINT	cout << "#################### After PID block 3" << endl;

	//Here I need to identify Pions and Kaons. This block is for the 2pi2k final state
	//if(trk->GetMostProbablePID() == 2) { //Pions
	if(i == posTrackIndex) {
	  //	  if(fabs(fPIDTPCPion[i]) < 4 && fabs(fPIDTPCPion[i]) < fabs(fPIDTPCKaon[i])) { 
	    //More restrictive PID   && fPIDTPCElectron[i] > 2 && fPIDTPCKaon[i] > 2 && fPIDTPCProton[i] > 2) {
	  if((trk->Pt() < 0.6 && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >= 3. && fabs(fPIDTPCPion[i]) < 1. && fabs(fPIDTPCKaon[i]) > 1.)) {
	    fRecTPCsignalK0sPion[0] = trk->GetTPCsignal();
	    qK0sPion[0] = trk->Charge();
	    vK0sPion[0].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
	    nK0sPion++;
	  }
	}
	else if(i == negTrackIndex) {
	  //	  if(fabs(fPIDTPCPion[i]) < 4 && fabs(fPIDTPCPion[i]) < fabs(fPIDTPCKaon[i])) { 
	    //More restrictive PID   && fPIDTPCElectron[i] > 2 && fPIDTPCKaon[i] > 2 && fPIDTPCProton[i] > 2) {
	  if((trk->Pt() < 0.6 && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >= 3. && fabs(fPIDTPCPion[i]) < 1. && fabs(fPIDTPCKaon[i]) > 1.)) {
	    fRecTPCsignalK0sPion[1] = trk->GetTPCsignal();
	    qK0sPion[1] = trk->Charge();
	    vK0sPion[1].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
	    nK0sPion++;
	  }
	}
	//	else if(fabs(fPIDTPCPion[i]) < 4 && fabs(fPIDTPCPion[i]) < fabs(fPIDTPCKaon[i])) { 
	  //More restrictive PID   && fPIDTPCElectron[i] > 2 && fPIDTPCKaon[i] > 2 && fPIDTPCProton[i] > 2) {
	else if((trk->Pt() < 0.6 && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >= 3. && fabs(fPIDTPCPion[i]) < 1. && fabs(fPIDTPCKaon[i]) > 1.)) {
	  fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
	  qPion[nPion] = trk->Charge();
	  vPion[nPion].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
	  nPion++;
	}
	//	else if(fabs(fPIDTPCKaon[i]) < 4 && fabs(fPIDTPCKaon[i]) < fabs(fPIDTPCPion[i])) { 
	  //More restrictive PID   && fPIDTPCElectron[i] > 2 && fPIDTPCPion[i] > 2 && fPIDTPCProton[i] > 2) {
	else if((trk->Pt() < 0.6 && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTPCPion[i]) > 3.) || (trk->Pt() >= 0.6 && trk->Pt() < 3. && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTOFPion[i]) > 3. && fPIDTOFPion[i] > -999.) || (trk->Pt() >= 3. && fabs(fPIDTPCKaon[i]) < 1. && fabs(fPIDTPCPion[i]) > 1.)) {
	  fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
	  qKaon[nKaon] = trk->Charge();
	  vKaon[nKaon].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),kaonMass);
	  nKaon++;
	}
	//***********************************
	//Must add code to cound number of missing tracks, accounting for those that might already be missing from the K0s daughters. Then add code outside the loop to assign the appropriate mass and increment counter (using which of nKaon, nPion, and UsedNMissingDaughter to determing what mass to assume for the missing track. If UsedNMissing Daughter is 1 then no more tracks can be missing but if it is zero then either the kaon or pion may be assumed depending on nPion and nKaon.
	//***********************************
	else if(trk->Pt() >= 0.6 && trk->Pt() < 3. && fPIDTOFPion[i] == -999.) { //assume missing track is pion.
	  missingTOFPIDK0s[i] = trackIndex[i];
	  nTracksWithoutTOFinfoK0s++; //Track index if TOF misses a track. Maybe due to missing TOF signal
	}

	if(nPion > 1 || nK0sPion > 2 || nKaon > 1) break;
      }
      for(int i=0;i<4;i++) { //If one kaon or one pion is missing due to missing TOF PID info assume it is the fourth.
	if(nK0sPion == 2 && nPion == 1 && nKaon == 0 && missingTOFPIDK0s[i] > 0 && nTracksWithoutTOFinfoK0s == 1) {
	  AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(missingTOFPIDK0s[i]));
	  qKaon[nKaon] = trk->Charge();
	  vKaon[nKaon].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),kaonMass);
	  nKaon++;
	} else if(((nKaon == 1 && nK0sPion == 1 && nPion == 1) || (nKaon == 1 && nK0sPion == 2 && nPion == 0)) && missingTOFPIDK0s[i] > 0 && nTracksWithoutTOFinfoK0s == 1) {
	  AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(missingTOFPIDK0s[i]));
	  qPion[nPion] = trk->Charge();
	  vPion[nPion].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
	  nPion++;
	}
      }

      //Histos to understand the events the code is seeing.
      if(nProngFound == 2) {
	fHistNpion->Fill(nPion);
	fHistNK0sPion->Fill(nK0sPion);
	fHistNkaon->Fill(nKaon);      
	fHistPiMinusK->Fill((nPion+nK0sPion-nKaon));
      }

      if( (nPion == 1) && (nK0sPion == 2) && (nKaon == 1) ) {
	fHistNeventsEtaCK0sChannel->Fill(8);
	if(qKaon[0] < 0 && ((qPion[0]+qK0sPion[0]+qK0sPion[1]) < 1)) fHistNeventsEtaCK0sChannel->Fill(9); //K- and sumQPions < 0 K-pi+pi-pi- or K-pi-pi-pi-
	if(qKaon[0] > 0 && ((qPion[0]+qK0sPion[0]+qK0sPion[1]) > 0)) fHistNeventsEtaCK0sChannel->Fill(10); //K+ and sumQPions > 0 K+pi+pi+pi- or K+pi+pi+pi+
	if(qKaon[0] < 0 && (qPion[0] > 0) && ((qK0sPion[0]+qK0sPion[1]) == 0)) fHistNeventsEtaCK0sChannel->Fill(11); //Requires K-Pi+Pi+Pi- (K- events)
	if(qKaon[0] > 0 && (qPion[0] < 0) && ((qK0sPion[0]+qK0sPion[1]) == 0)) fHistNeventsEtaCK0sChannel->Fill(12); //Requires K+Pi-Pi+Pi- (K+ events)
	//Continue for good cases
	if( ((qKaon[0] < 0) && (qPion[0] > 0) && ((qK0sPion[0]+qK0sPion[1]) == 0)) || ((qKaon[0] > 0) && (qPion[0] < 0) && ((qK0sPion[0]+qK0sPion[1]) == 0)) ) {

	  //PRINT	  cout << "Filling K0s histos" << endl;

	  fHistNeventsEtaCK0sChannel->Fill(13);
	  //Now we need to use the V0s to identify a K0s and select best K0s candidate.
	  //Individual Track Ptś
	  fK0sPosDaughterPt->Fill(vK0sPion[0].Pt());
	  fK0sNegDaughterPt->Fill(vK0sPion[1].Pt());
	  fK0sPosVsNegDaughterPt->Fill(vK0sPion[1].Pt(),vK0sPion[0].Pt());
	  fPionK0sChannelPt->Fill(vPion[0].Pt());
	  fKaonK0sChannelPt->Fill(vKaon[0].Pt());
	  //Compute K0s info
	  vK0s = vK0sPion[0] + vK0sKaon[1];
	  fK0sPtVsMinv->Fill(vK0s.M(),vK0s.Pt());
	  //fK0sMinv->Fill(vK0s.M());
	  //Compute PiK info
	  vKPiK0sChannel = vPion[0] + vKaon[0];
	  fKPiPtVsMinvK0sChannel->Fill(vKPiK0sChannel.M(),vKPiK0sChannel.Pt());
	  //fKPiMinvK0sChannel->Fill(vKPiK0sChannel.M());
	  //Dalitz plot K0s vs PiK combo
	  fMK0sVsMKPiK0sChannel->Fill(vKPiK0sChannel.M(),vK0s.M());
	  //Compute EtaC info
	  vCandidate = vK0sPion[0] + vK0sKaon[1] + vPion[0] + vKaon[0];
	  fEtaCPtVsMinvK0sChannel->Fill(vCandidate.M(),vCandidate.Pt());
	  //fEtaCMinvK0sChannel->Fill(vCandidate.M());
	  fK0sDecayLength->Fill(v0->DecayLength(fAODVertex));
	}
      }
    }
  }
  

  //  cout << "##### End of RunAODHist()" << endl;

  PostData(4, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunAODtree()
{

  //  cout << "############### Beginning of RunAODtree()" << endl;
  /*
  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  if(isMC) RunAODMC(aod);

  //input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  fDataFilnam->Clear();
  fDataFilnam->SetString(filnam);
  fEvtNum = ((TTree*) GetInputData(0))->GetTree()->GetReadEntry();
  fRunNum = aod ->GetRunNumber();

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  if(fTracking == 0) {
    fTrigger[0]  = trigger.Contains("CCUP4-B"); // Central UPC Pb-Pb 2011
    fTrigger[1]  = trigger.Contains("CCUP2-B"); // Double gap
    fTrigger[2]  = trigger.Contains("CCUP7-B"); // Central UPC p-Pb 2013
    fTrigger[3]  = trigger.Contains("CINT1-B"); // MB trigger
    fTrigger[4]  = trigger.Contains("CTEST58-B"); // *0VBA *0VBC *0UBA *0UBC 0SH1
    fTrigger[5]  = trigger.Contains("CTEST59-B"); // *0VBA *0VBC *0UBA *0UBC 0STP
    fTrigger[6]  = trigger.Contains("CTEST60-B"); // *0VBA *0VBC *0UBA *0UBC 0OM2
    fTrigger[7]  = trigger.Contains("CTEST61-B"); // *0VBA *0VBC *0UBA *0UBC 0OMU
    fTrigger[8]  = trigger.Contains("CCUP8-B"); //*0VBA *0VBC *0UBA *0UBC 0STP 0OMU
    fTrigger[9]  = trigger.Contains("CCUP9-B"); //*0VBA *0VBC *0UBA *0UBC 0STP
    fTrigger[10]  = trigger.Contains("CCUP10-B"); //*0VBA *0VBC *0UBA *0UBC 0SH1
    fTrigger[11]  = trigger.Contains("CCUP11-B"); //*0UBA *0UBC 0STP 0OMU
    fTrigger[12]  = trigger.Contains("CCUP12-B"); //*0UBA *0UBC 0STP
    fTrigger[13]  = trigger.Contains("CTRUE-B"); //Unbiased trigger
  }
  if(fTracking == 8) {
    fTrigger[0]  = trigger.Contains("CMUP10-B");    // *0VBA *0UBA *0UBC 0MSL                       
    fTrigger[1]  = trigger.Contains("CMUP11-B");    // !0VBA & !0UBA & !0UBC & 0MUL                 
    fTrigger[2]  = trigger.Contains("CMUP12-B");    // !0VBA & !0UBA & !0UBC & 0MSL & 0SMB
    fTrigger[3]  = trigger.Contains("CMUP14-B");   // 0MSL & !0VBA & !0UBA
    fTrigger[4]  = trigger.Contains("CMUP15-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL
    fTrigger[5]  = trigger.Contains("CMUP16-B");   // 0MSL *0VBA *0UBA *0UGC *0VGA
    fTrigger[6]  = trigger.Contains("CMUP17-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL *0UGC *0VGA
    fTrigger[7]  = trigger.Contains("CMUP21-B");   // *0VBA *0UBA *0VBC 0SH1 *0SH2 *0UGC *0VGA
    fTrigger[8]  = trigger.Contains("CMUP22-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MSL 0SMB
    fTrigger[9]  = trigger.Contains("CMUP23-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MUL
  }

  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if(!isMC && !isTriggered ) return;
  
  //Physics selection
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if((fTrigger[3] || fTrigger[4]) &&((selectionMask & AliVEvent::kMB) == AliVEvent::kMB)) fIsPhysicsSelected = kTRUE;
  else fIsPhysicsSelected = kFALSE;

  //trigger inputs
  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  fL1inputs = aod->GetHeader()->GetL1TriggerInputs();  

  //TOF trigger mask
  const AliTOFHeader *tofH = aod->GetTOFHeader();
  fTOFmask = tofH->GetTriggerMask();

  //Event identification
  fPerNum = aod ->GetPeriodNumber();
  fOrbNum = aod ->GetOrbitNumber();
  fBCrossNum = aod ->GetBunchCrossNumber();

  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();
  fVtxPos[0] = fAODVertex->GetX();
  fVtxPos[1] = fAODVertex->GetY();
  fVtxPos[2] = fAODVertex->GetZ();
  Double_t CovMatx[6];
  fAODVertex->GetCovarianceMatrix(CovMatx); 
  fVtxErr[0] = CovMatx[0];
  fVtxErr[1] = CovMatx[1];
  fVtxErr[2] = CovMatx[2];
  fVtxChi2 = fAODVertex->GetChi2();
  fVtxNDF = fAODVertex->GetNDF();
  
  //SPD primary vertex
  AliAODVertex *fSPDVertex = aod->GetPrimaryVertexSPD();
  fSpdVtxContrib = fSPDVertex->GetNContributors();
  fSpdVtxPos[0] = fSPDVertex->GetX();
  fSpdVtxPos[1] = fSPDVertex->GetY();
  fSpdVtxPos[2] = fSPDVertex->GetZ();

  //Tracklets
  fNtracklets = aod->GetTracklets()->GetNumberOfTracklets();

  //VZERO, ZDC, AD
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  AliAODZDC *fZDCdata = aod->GetZDCData();
  AliAODAD *fADdata = aod ->GetADData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  
  if(fADdata){
  	fADAdecision = fADdata->GetADADecision();
  	fADCdecision = fADdata->GetADCDecision();
  }
  
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  fZPAenergy = fZDCdata->GetZPATowerEnergy()[0];
  fZPCenergy = fZDCdata->GetZPCTowerEnergy()[0];  
  fZDCAtime = fZDCdata->GetZNATime();
  fZDCCtime = fZDCdata->GetZNCTime();
  
  fNLooseTracks = 0;
  
  //Track loop - loose cuts
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(fTracking == 0){
      if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliAODTrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 20)continue;
      fNLooseTracks++;
      }
    if(fTracking == 1){
      if(!(trk->TestFilterBit(1<<1))) continue;
      
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      fNLooseTracks++;
      }  
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  
  //EtaC->Pi+Pi-K+K- channel
   nGoodTracks = 0;
   //Four track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    
    if(fTracking == 0){
      if(!(trk->TestFilterBit(1<<0))) continue;
      if(!(trk->GetStatus() & AliAODTrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = 4*(0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;

      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!(trk->TestFilterBit(1<<1))) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->Chi2perNDF() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
      
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
      
  fEtaCAODTracks->Clear("C");  
  if(nGoodTracks == 4){

  	  TDatabasePDG *pdgdat = TDatabasePDG::Instance();

	  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
	  Double_t kaonMass = partKaon->Mass();
	  
	  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
	  Double_t pionMass = partPion->Mass();
	  
	  TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
	  Double_t kStarMass = partKstar->Mass();
	  Double_t kStarWidth = partKstar->Width();
	  
	  TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
	  Double_t k0ShortMass = partK0short->Mass();
	  Double_t k0ShortWidth = partK0short->Width();
 
  	  Double_t KFcov[21];
  	  Double_t KFpar[6];
	  Double_t KFmass = pionMass;
	  Double_t fRecTPCsignal;
  	  AliKFParticle *KFpart[4];
  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(aod->GetMagneticField()); 
	  	  
  	  for(Int_t i=0; i<4; i++){
	    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
	    if(!trk) AliFatal("Not a standard AOD");
	    
	    if(fAODVertex->HasDaughter(trk) && trk->GetUsedForVtxFit())fIsVtxContributor[i] = kTRUE;
	    else fIsVtxContributor[i] = kFALSE;
	    
	    Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
	    AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
	    if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
	    delete trk_clone;
	    
	    new((*fEtaCAODTracks)[i]) AliAODTrack(*trk);
	    ((AliAODTrack*)((*fEtaCAODTracks)[i]))->SetDCA(dca[0],dca[1]);//to get DCAxy trk->DCA(); to get DCAz trk->ZAtDCA();
	    
	    
	    fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
	    fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
	    fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
	    fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
	    fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
	    
	    fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
	    fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
	    fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
	    fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
	    fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
	    
	    trk->GetPosition(KFpar);
	    trk->PxPyPz(KFpar+3);
	    trk->GetCovarianceXYZPxPyPz(KFcov);
	    
	    //Do the PID, set masses and create KFparticles		
	    if(fPIDTPCPion[i] < 2 && fPIDTPCMuon[i] > 2 && fPIDTPCElectron[i] > 2 && fPIDTPCKaon[i] > 2 && fPIDTPCProton[i] > 2) {
	      KFmass = pionMass;
	    }
	    //else if(trk->GetMostProbablePID() == 3) { //Kaons
	    else if(fPIDTPCKaon[i] < 2 && fPIDTPCMuon[i] > 2 && fPIDTPCElectron[i] > 2 && fPIDTPCPion[i] > 2 && fPIDTPCProton[i] > 2) {
	      KFmass = kaonMass;
	    }
	    
	    
	    KFpart[i] = new AliKFParticle();
	    KFpart[i]->SetField(aod->GetMagneticField());
	    KFpart[i]->AliKFParticleBase::Initialize(KFpar,KFcov,(Int_t) trk->Charge(), KFmass);
	    KFvtx->AddDaughter(*KFpart[i]); 
	    
	  }
	  fKfVtxPos[0]= KFvtx->GetX();
	  fKfVtxPos[1]= KFvtx->GetY();
	  fKfVtxPos[2]= KFvtx->GetZ();
	  for(UInt_t i=0; i<4; i++)delete KFpart[i];
	  delete KFvtx; 
	  if(!isMC){
	    fEtaCK0sChannelTree->Fill(); 
	    fEtaCTree ->Fill();
	  }
  }
  
  if(isMC){
  	fEtaCK0sChannelTree ->Fill();
	fEtaCTree ->Fill();
  }
  */  
  //  cout << "############### End of RunAODtree()" << endl;

  //  PostData(1, fEtaCK0sChannelTree);
  //PostData(2, fEtaCTree);

}//RunAOD


//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunAODMC(AliAODEvent *aod)
{

  for(Int_t i=0; i<ntrg; i++) fTriggerInputsMC[i] = kFALSE;
  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  fTriggerInputsMC[0] = fL0inputs & (1 << 0);   //0VBA VZERO A
  fTriggerInputsMC[1] = fL0inputs & (1 << 1);   //0VBC VZERO C
  fTriggerInputsMC[2] = fL0inputs & (1 << 11);  //0OMU TOF two hits with topology
  fTriggerInputsMC[3] = fL0inputs & (1 << 23);	//0OM2 TOF two hits
						
  //SPD inputs
  const AliAODTracklets *mult = aod->GetMultiplicity();
  Int_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=0;
  Int_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=0;

  Int_t nInner(0), nOuter(0);
  for (Int_t i(0); i<1200; ++i) {
    Bool_t isFired(mult->TestFastOrFiredChips(i));
    if (i<400) {
      vPhiInner[i/20] += isFired;
      nInner += isFired;
    } else {
      vPhiOuter[(i-400)/20] += isFired;
      nOuter += isFired;
    }
  }
 
  Int_t fired(0);
  for (Int_t i(0); i<10; ++i) {
    for (Int_t j(0); j<2; ++j) {
      const Int_t k(2*i+j);
      fired += ((   vPhiOuter[k]    || vPhiOuter[k+1]       ||
                    vPhiOuter[k+2]      )
                && (vPhiOuter[k+20] || vPhiOuter[(k+21)%40] ||
                    vPhiOuter[(k+22)%40])
                && (vPhiInner[i]    || vPhiInner[i+1]       )
                && (vPhiInner[i+10] || vPhiInner[(i+11)%20]));
    }
  }
  //0SMB - At least one hit in SPD
  if (nOuter > 0 || nInner > 0) fTriggerInputsMC[4] = kTRUE;
  //0SM2 - Two hits on outer layer
  if (nOuter > 1) fTriggerInputsMC[5] = kTRUE;
  //0STP - Topological SPD trigger (two pairs)
  if (fired != 0) fTriggerInputsMC[6] = kTRUE;
  //0SH1 - More then 6 hits on outer layer
  if (nOuter >= 7) fTriggerInputsMC[7] = kTRUE;
  

  fGenPart->Clear("C");

  TClonesArray *arrayMC = (TClonesArray*) aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC) return;

  Int_t nmc=0;
  //loop over mc particles
  for(Int_t imc=0; imc<arrayMC->GetEntriesFast(); imc++) {
    AliAODMCParticle *mcPart = (AliAODMCParticle*) arrayMC->At(imc);
    if(!mcPart) continue;

    if(mcPart->GetMother() >= 0) continue;

    TParticle *part = (TParticle*) fGenPart->ConstructedAt(nmc++);
    part->SetMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->E());
    part->SetPdgCode(mcPart->GetPdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles
  
  AliAODMCHeader *mcHeader = (AliAODMCHeader*) aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) return;
  
  fMCVtxPos[0] = mcHeader->GetVtxX();
  fMCVtxPos[1] = mcHeader->GetVtxY();
  fMCVtxPos[2] = mcHeader->GetVtxZ(); 

}//RunAODMC


//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunESDtrig()
{

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  fRunNum = esd ->GetRunNumber();
  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  if(trigger.Contains("CCUP4-B")) fHistCcup4TriggersPerRun->Fill(fRunNum); //CCUP4 triggers
  if(trigger.Contains("CCUP7-B")) fHistCcup7TriggersPerRun->Fill(fRunNum); //CCUP7 triggers
  if(trigger.Contains("CCUP2-B")) fHistCcup2TriggersPerRun->Fill(fRunNum); //CCUP2 triggers
  
  if(trigger.Contains("CINT1-B")) fHistCint1TriggersPerRun->Fill(fRunNum); //CINT1 triggers
  
  if(trigger.Contains("CTEST58-B")) fHistCTest58TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST59-B")) fHistCTest59TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST60-B")) fHistCTest60TriggersPerRun->Fill(fRunNum); //CTEST triggers
  if(trigger.Contains("CTEST61-B")) fHistCTest61TriggersPerRun->Fill(fRunNum); //CTEST triggers
  
  if(trigger.Contains("CCUP8-B")) fHistCcup8TriggersPerRun->Fill(fRunNum); //CCUP8 triggers
  if(trigger.Contains("CCUP9-B")) fHistCcup9TriggersPerRun->Fill(fRunNum); //CCUP9 triggers
  
  fL0inputs = esd->GetHeader()->GetL0TriggerInputs();
  if(trigger.Contains("CINT1-B") && (fL0inputs & (1 << 3))) fHistC0tvxAndCint1TriggersPerRun->Fill(fRunNum); //0TVX triggers in CINT1 events
  
  if(trigger.Contains("CVLN_B2-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - synchronously downscaled
  if(trigger.Contains("CVLN_R1-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - randomly downscaled
  
  if(esd->GetHeader()->IsTriggerInputFired("1ZED")) fHistZedTriggersPerRun->Fill(fRunNum); //1ZED trigger inputs
  
   //MB, Central and SemiCentral triggers
  AliCentrality *centrality = esd->GetCentrality();
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  //Double_t percentile = centrality->GetCentralityPercentile("V0M");
  Double_t percentile = centrality->GetCentralityPercentileUnchecked("V0M");
  
  if(((selectionMask & AliVEvent::kMB) == AliVEvent::kMB) && percentile<=80 && percentile>=0) fHistMBTriggersPerRun->Fill(fRunNum);
  
  if(((selectionMask & AliVEvent::kCentral) == AliVEvent::kCentral) && percentile<=6 && percentile>=0 && (trigger.Contains("CVHN_R2-B"))) fHistCentralTriggersPerRun->Fill(fRunNum);

  if(((selectionMask & AliVEvent::kSemiCentral) == AliVEvent::kSemiCentral) && percentile<=50 && percentile>=15) fHistSemiCentralTriggersPerRun->Fill(fRunNum);

  
PostData(3, fListTrig);

}
//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunESDhist()
{

  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Double_t kaonMass = partKaon->Mass();
    
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();

  TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
  Double_t kStarMass = partKstar->Mass();
  Double_t kStarWidth = partKstar->Width();

  TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
  Double_t k0ShortMass = partK0short->Mass();
  Double_t k0ShortWidth = partK0short->Width();

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  fHistNeventsEtaCK0sChannel->Fill(1);
  fHistNeventsEtaC->Fill(1);

  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  if(!isMC && !trigger.Contains("CCUP") ) return;
  
  fHistNeventsEtaCK0sChannel->Fill(2);
  fHistNeventsEtaC->Fill(2);
  
  
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  if(fZDCdata->IsZNAhit()) fZDCAtime= fZDCdata->GetZDCTDCCorrected(12,0);
  else fZDCAtime=-666;
  if(fZDCdata->IsZNChit()) fZDCCtime= fZDCdata->GetZDCTDCCorrected(10,0);
  else fZDCCtime=-666;
  
  if(trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1);
  if(fZNAenergy < 8200 && fZNCenergy < 8200) fHistZDCCuts->Fill(2);
  if(fZNAenergy < 683 && fZNCenergy < 683) fHistZDCCuts->Fill(3);
  if(fZDCAtime == -666 && fZDCCtime == -666) fHistZDCCuts->Fill(4);

  //primary vertex
  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
  fVtxContrib = fESDVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsEtaCK0sChannel->Fill(3);
  fHistNeventsEtaC->Fill(3);

  //VZERO, ZDC
  AliESDVZERO *fV0data = esd->GetVZEROData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliESDVZERO::kV0Empty || fV0Cdecision != AliESDVZERO::kV0Empty) return;

  if( fZNAenergy > 8200 || fZNCenergy > 8200) return;
  
  fHistNeventsEtaCK0sChannel->Fill(4);
  fHistNeventsEtaC->Fill(4);

   Int_t nGoodTracks=0;
  //Two tracks loop
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vPion[4], vKaon[4], vKstar[4], vCandidate;
  Short_t qKaon[4], qPion[4], qKstar[4];
  UInt_t nKaon=0, nPion=0;
  Double_t fRecTPCsignalPion[5], fRecTPCsignalKaon[5];



  nGoodTracks = 0; nPion=0; nKaon=0;
  
    //Four Track loop
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
      Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
      if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
      trk->GetImpactParameters(dca[0],dca[1]);
      if(TMath::Abs(dca[1]) > 2) continue;
      
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 4) break;   
  }//Track loop
  
  if(nGoodTracks == 4){ //The nSpdHits>1 requirement is not used for ESD analysis apparently.
  	  fHistNeventsEtaC->Fill(5);
  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(trackIndex[i]);
		
		//Start of my code
		//Get nsigma info for PID
		fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
		fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
		fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
		fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
		fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);

		//Here I need to identify Pions and Kaons. This block is for the 2pi2k final state
		//if(trk->GetMostProbablePID() == 2) { //Pions
		if(fPIDTPCPion[i] < 2 && fPIDTPCMuon[i] > 2 && fPIDTPCElectron[i] > 2 && fPIDTPCKaon[i] > 2 && fPIDTPCProton[i] > 2) {
		  fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
		  qPion[nPion] = trk->Charge();
		  vPion[nPion].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),pionMass);
		  nPion++;
		}
		//else if(trk->GetMostProbablePID() == 3) { //Kaons
		else if(fPIDTPCKaon[i] < 2 && fPIDTPCMuon[i] > 2 && fPIDTPCElectron[i] > 2 && fPIDTPCPion[i] > 2 && fPIDTPCProton[i] > 2) {
		  fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
		  qKaon[nKaon] = trk->Charge();
		  vKaon[nKaon].SetPtEtaPhiM(trk->Pt(),trk->Eta(),trk->Phi(),kaonMass);
		  nKaon++;
		}

		if(nPion > 2 || nKaon > 2) break;
	  }

	  //Analyze good events, fill the histos.
	  if( (nPion == 2) && (nKaon == 2) ) {
	    fHistNeventsEtaC->Fill(7);
	    if(qKaon[0]*qKaon[1] > 0) fHistNeventsEtaC->Fill(8);
	    if(qPion[0]*qPion[1] > 0) fHistNeventsEtaC->Fill(9);
	    if((qKaon[0]*qKaon[1] > 0) && (qPion[0]*qPion[1] > 0)) fHistNeventsEtaC->Fill(10);
	    if((qKaon[0]*qKaon[1] < 0) && (qPion[0]*qPion[1] < 0)) { 
	      fHistNeventsEtaC->Fill(11);
	      if(vPion[0].M() == pionMass && vPion[1].M() == pionMass && vKaon[0].M() == kaonMass && vKaon[1].M() == kaonMass) {
		fHistNeventsEtaC->Fill(12);
		vCandidate = vPion[0] + vPion[1] + vKaon[0] + vKaon[1];

		//Get masses of potential intermediate Kstar(892)'s
		if((qPion[0]*qKaon[0] < 0) && (qPion[1]*qKaon[1] < 0)) {
		  vKstar[0] = vPion[0] + vKaon[0];
		  vKstar[1] = vPion[1] + vKaon[1];
		}
		else if ((qPion[0]*qKaon[1] < 0) && (qPion[1]*qKaon[0] < 0)) {
		  vKstar[0] = vPion[0] + vKaon[1];
		  vKstar[1] = vPion[1] + vKaon[0];
		}

		//Fill Dalitz plot with PiK masses Pi-K+ vs Pi+K-
		if(qKaon[0] < 0) fMPiKvsMPiK->Fill(vKstar[0].M(),vKstar[1].M());
		else fMPiKvsMPiK->Fill(vKstar[1].M(),vKstar[0].M());

		//Fill histos
		//2 Kstar case
		if((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
		   (vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth))) {
		  fHistNeventsEtaC->Fill(13);
		  if(qPion[0] > 0 && qPion[1] < 0) {
		    f2KstarPtPiPlus->Fill(vPion[0].Pt());
		    f2KstarPtPiMinus->Fill(vPion[1].Pt());
		  } else {
		    f2KstarPtPiPlus->Fill(vPion[1].Pt());
		    f2KstarPtPiMinus->Fill(vPion[0].Pt());
		  }
		  if(qKaon[0] > 0 && qKaon[1] < 0) {
		    f2KstarPtKPlus->Fill(vKaon[0].Pt());
		    f2KstarPtKMinus->Fill(vKaon[1].Pt());
		  } else {
		    f2KstarPtKPlus->Fill(vKaon[1].Pt());
		    f2KstarPtKMinus->Fill(vKaon[0].Pt());
		  }
		  //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		  f2KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		  f2KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		  //Fill intermediate Kstar histos, one for each candidate
		  f2KstarPtVsMinvFirstKstar->Fill(vKstar[0].M(),vKstar[0].Pt());
		  f2KstarPtVsMinvSecondKstar->Fill(vKstar[1].M(),vKstar[1].Pt());
		  //    f2KstarMinvFirstKstar->Fill(vKstar[0].M());
		  //    f2KstarMinvSecondKstar->Fill(vKstar[1].M());
		  //Fill EtaC histos
		  f2KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		  //    f2KstarMinvEtaC->Fill(vCandidate.M());
		}
		//1 Kstar case
		else if ((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
			 ((vKstar[1].M() > (kStarMass + kStarWidth)) || (vKstar[1].M() < (kStarMass - kStarWidth)))) {
		  //Fill using first Kstar candidate
		  fHistNeventsEtaC->Fill(14);
		  if(qPion[0] > 0 && qPion[1] < 0) {
		    f1KstarPtPiPlus->Fill(vPion[0].Pt());
		    f1KstarPtPiMinus->Fill(vPion[1].Pt());
		  } else {
		    f1KstarPtPiPlus->Fill(vPion[1].Pt());
		    f1KstarPtPiMinus->Fill(vPion[0].Pt());
		  }
		  if(qKaon[0] > 0 && qKaon[1] < 0) {
		    f1KstarPtKPlus->Fill(vKaon[0].Pt());
		    f1KstarPtKMinus->Fill(vKaon[1].Pt());
		  } else {
		    f1KstarPtKPlus->Fill(vKaon[1].Pt());
		    f1KstarPtKMinus->Fill(vKaon[0].Pt());
		  }
		  //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		  f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		  f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		  //Fill intermediate Kstar histos, one for each candidate
		  f1KstarPtVsMinvKstar->Fill(vKstar[0].M(),vKstar[0].Pt());
		  f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[1].M(),vKstar[1].Pt());
		  //    f1KstarMinvKstar->Fill(vKstar[0].M());
		  //    f1KstarMinvOtherPiKcombo->Fill(vKstar[1].M());
		  //Fill EtaC histos
		  f1KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		  //    f1KstarMinvEtaC->Fill(vCandidate.M());
		}
		else if ((vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth)) &&
			 ((vKstar[0].M() > (kStarMass + kStarWidth)) || (vKstar[0].M() < (kStarMass - kStarWidth)))) {
		  //Fill using second Kstar candidate
		  fHistNeventsEtaC->Fill(14);
		  fHistNeventsEtaC->Fill(14);
		  if(qPion[0] > 0 && qPion[1] < 0) {
		    f1KstarPtPiPlus->Fill(vPion[0].Pt());
		    f1KstarPtPiMinus->Fill(vPion[1].Pt());
		  } else {
		    f1KstarPtPiPlus->Fill(vPion[1].Pt());
		    f1KstarPtPiMinus->Fill(vPion[0].Pt());
		  }
		  if(qKaon[0] > 0 && qKaon[1] < 0) {
		    f1KstarPtKPlus->Fill(vKaon[0].Pt());
		    f1KstarPtKMinus->Fill(vKaon[1].Pt());
		  } else {
		    f1KstarPtKPlus->Fill(vKaon[1].Pt());
		    f1KstarPtKMinus->Fill(vKaon[0].Pt());
		  }
		  //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		  f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		  f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		  //Fill intermediate Kstar histos, one for each candidate
		  f1KstarPtVsMinvKstar->Fill(vKstar[1].M(),vKstar[1].Pt());
		  f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[0].M(),vKstar[0].Pt());
		  //    f1KstarMinvKstar->Fill(vKstar[1].M());
		  //    f1KstarMinvOtherPiKcombo->Fill(vKstar[0].M());
		  //Fill EtaC histos
		  f1KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		  //    f1KstarMinvEtaC->Fill(vCandidate.M());
		}
		  //0 Kstar case
		else {
		  fHistNeventsEtaC->Fill(15);
		  if(qPion[0] > 0 && qPion[1] < 0) {
		    f0KstarPtPiPlus->Fill(vPion[0].Pt());
		    f0KstarPtPiMinus->Fill(vPion[1].Pt());
		  } else {
		    f0KstarPtPiPlus->Fill(vPion[1].Pt());
		    f0KstarPtPiMinus->Fill(vPion[0].Pt());
		  }
		  if(qKaon[0] > 0 && qKaon[1] < 0) {
		    f0KstarPtKPlus->Fill(vKaon[0].Pt());
		    f0KstarPtKMinus->Fill(vKaon[1].Pt());
		  } else {
		    f0KstarPtKPlus->Fill(vKaon[1].Pt());
		    f0KstarPtKMinus->Fill(vKaon[0].Pt());
		  }
		  //Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
		  f0KstarTPCsignalPion->Fill(fRecTPCsignalPion[0],fRecTPCsignalPion[1]);
		  f0KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0],fRecTPCsignalKaon[1]);
		  //Fill intermediate Kstar histos, one for each candidate
		  f0KstarPtVsMinvFirstPiKcombo->Fill(vKstar[0].M(),vKstar[0].Pt());
		  f0KstarPtVsMinvSecondPiKcombo->Fill(vKstar[1].M(),vKstar[1].Pt());
		  //    f0KstarMinvFirstPiKcombo->Fill(vKstar[0].M());
		  //    f0KstarMinvSecondPiKcombo->Fill(vKstar[1].M());
		  //Fill EtaC histos
		  f0KstarPtVsMinvEtaC->Fill(vCandidate.M(),vCandidate.Pt());
		  //    f0KstarMinvEtaC->Fill(vCandidate.M());
		}
	      }
	    }
	  }
  }

  //End of my EtaC->PiPiKK code


  
  PostData(4, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunESDtree()
{
  /*
  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;
  
  if(isMC) RunESDMC(esd);

  //input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  fDataFilnam->Clear();
  fDataFilnam->SetString(filnam);
  fEvtNum = ((TTree*) GetInputData(0))->GetTree()->GetReadEntry();
  fRunNum = esd->GetRunNumber();

   //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  fTrigger[0]  = trigger.Contains("CCUP4-B"); // Central UPC Pb-Pb 2011
  fTrigger[1]  = trigger.Contains("CCUP2-B"); // Double gap
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // Central UPC p-Pb 2013
  fTrigger[3]  = trigger.Contains("CINT1-B"); // MB trigger
  fTrigger[4]  = trigger.Contains("CTEST58-B"); // *0VBA *0VBC *0UBA *0UBC 0SH1
  fTrigger[5]  = trigger.Contains("CTEST59-B"); // *0VBA *0VBC *0UBA *0UBC 0STP
  fTrigger[6]  = trigger.Contains("CTEST60-B"); // *0VBA *0VBC *0UBA *0UBC 0OM2
  fTrigger[7]  = trigger.Contains("CTEST61-B"); // *0VBA *0VBC *0UBA *0UBC 0OMU
  fTrigger[8]  = trigger.Contains("CCUP8-B"); //*0VBA *0VBC *0UBA *0UBC 0STP 0OMU
  fTrigger[9]  = trigger.Contains("CCUP9-B"); //*0VBA *0VBC *0UBA *0UBC 0STP
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if(!isMC && !isTriggered ) return;
  
  //Physics selection
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if((fTrigger[3] || fTrigger[4]) &&((selectionMask & AliVEvent::kMB) == AliVEvent::kMB)) fIsPhysicsSelected = kTRUE;
  else fIsPhysicsSelected = kFALSE;
  
  //trigger inputs
  fL0inputs = esd->GetHeader()->GetL0TriggerInputs();
  fL1inputs = esd->GetHeader()->GetL1TriggerInputs();
  
  //Event identification
  fPerNum = esd->GetPeriodNumber();
  fOrbNum = esd->GetOrbitNumber();
  fBCrossNum = esd->GetBunchCrossNumber();

  //TOF trigger mask
  const AliTOFHeader *tofH = esd->GetTOFHeader();
  fTOFmask = tofH->GetTriggerMask();

  //primary vertex
  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
  fVtxContrib = fESDVertex->GetNContributors();
  fVtxPos[0] = fESDVertex->GetX();
  fVtxPos[1] = fESDVertex->GetY();
  fVtxPos[2] = fESDVertex->GetZ();
  Double_t CovMatx[6];
  fESDVertex->GetCovarianceMatrix(CovMatx); 
  fVtxErr[0] = CovMatx[0];
  fVtxErr[1] = CovMatx[1];
  fVtxErr[2] = CovMatx[2];
  fVtxChi2 = fESDVertex->GetChi2();
  fVtxNDF = fESDVertex->GetNDF();
    
  //SPD primary vertex
  AliESDVertex *fSPDVertex = (AliESDVertex*) esd->GetPrimaryVertexSPD();
  fSpdVtxContrib = fSPDVertex->GetNContributors();
  fSpdVtxPos[0] = fSPDVertex->GetX();
  fSpdVtxPos[1] = fSPDVertex->GetY();
  fSpdVtxPos[2] = fSPDVertex->GetZ();

  //Tracklets
  fNtracklets = esd->GetMultiplicity()->GetNumberOfTracklets();

  //VZERO, ZDC, AD
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  AliESDAD *fADdata = esd->GetADData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fADdata){
  	fADAdecision = fADdata->GetADADecision();
  	fADCdecision = fADdata->GetADCDecision();
	}
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  fZPAenergy = fZDCdata->GetZPATowerEnergy()[0];
  fZPCenergy = fZDCdata->GetZPCTowerEnergy()[0];
  if(fZDCdata->IsZNAhit()) fZDCAtime= fZDCdata->GetZDCTDCCorrected(12,0);
  else fZDCAtime=-666;
  if(fZDCdata->IsZNChit()) fZDCCtime= fZDCdata->GetZDCTDCCorrected(10,0);
  else fZDCCtime=-666;
  
  fNLooseTracks = 0;
  
  //Track loop - loose cuts
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;
    
    if(fTracking == 0){
      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 20)continue;
      fNLooseTracks++;
      }
    if(fTracking == 1){
      if(!(trk->GetStatus() & AliESDtrack::kITSpureSA) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      fNLooseTracks++;
      }  
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  
  
  nGoodTracks = 0;
  //Four track loop
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;
          
    if(fTracking == 0){
      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
      Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
      if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
      trk->GetImpactParameters(dca[0],dca[1]);
      if(!isMC){
      	if(TMath::Abs(dca[1]) > 2) continue;
      	Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
        if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
	}
      
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!(trk->GetStatus() & AliESDtrack::kITSpureSA) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->GetTPCchi2()/trk->GetITSNcls() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
      
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
      
      if(nGoodTracks > 4) break;   
  }//Track loop
  
  fEtaCESDTracks->Clear("C");
  if(nGoodTracks == 4){
  
    TDatabasePDG *pdgdat = TDatabasePDG::Instance();
    
    TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
    Double_t kaonMass = partKaon->Mass();
    
    TParticlePDG *partPion = pdgdat->GetParticle( 211 );
    Double_t pionMass = partPion->Mass();
    
    TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
    Double_t kStarMass = partKstar->Mass();
    Double_t kStarWidth = partKstar->Width();
    
    TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
    Double_t k0ShortMass = partK0short->Mass();
    Double_t k0ShortWidth = partK0short->Width();
  
  	  Double_t KFcov[21];
  	  Double_t KFpar[6];
	  Double_t KFmass = pionMass;
	  Double_t fRecTPCsignal;
  	  AliKFParticle *KFpart[2];
  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(esd->GetMagneticField()); 

  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(trackIndex[i]);
		
		if(fESDVertex->UsesTrack(trackIndex[i]))fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fEtaCESDTracks)[i]) AliESDtrack(*trk);
		
		fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
		fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
		fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
		fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
		fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
		
		trk->GetXYZ(KFpar);
    		trk->PxPyPz(KFpar+3);
    		trk->GetCovarianceXYZPxPyPz(KFcov);
		
		if(trk->Pt() > 1){   
      			fRecTPCsignal = trk->GetTPCsignal();      
      			if(fRecTPCsignal > 40 && fRecTPCsignal < 70) KFmass = muonMass;
      			if(fRecTPCsignal > 70 && fRecTPCsignal < 100)KFmass = electronMass;
			}
		else KFmass = pionMass;
		
		KFpart[i] = new AliKFParticle();
    		KFpart[i]->SetField(esd->GetMagneticField());
    		KFpart[i]->AliKFParticleBase::Initialize(KFpar,KFcov,(Int_t) trk->Charge(), KFmass);
		KFvtx->AddDaughter(*KFpart[i]); 		
  		}
		
  fKfVtxPos[0]= KFvtx->GetX();
  fKfVtxPos[1]= KFvtx->GetY();
  fKfVtxPos[2]= KFvtx->GetZ();
  for(UInt_t i=0; i<2; i++)delete KFpart[i];
  delete KFvtx; 
  
  if(!isMC) {
  fEtaCK0sChannelTree->Fill();
  fEtaCTree ->Fill();
  }
  }
  
  if(isMC){
  	fJPsiTree ->Fill();
	fEtaCTree ->Fill();
  }
  */ 
  PostData(1, fEtaCK0sChannelTree);
  PostData(2, fEtaCTree);

}//RunESD


//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunESDMC(AliESDEvent* esd)
{
  for(Int_t i=0; i<ntrg; i++) fTriggerInputsMC[i] = kFALSE;
  fTriggerInputsMC[0] = esd->GetHeader()->IsTriggerInputFired("0VBA"); //VZERO A
  fTriggerInputsMC[1] = esd->GetHeader()->IsTriggerInputFired("0VBC"); //VZERO C
  fTriggerInputsMC[2] = esd->GetHeader()->IsTriggerInputFired("0OMU"); //TOF two hits with topology
  fTriggerInputsMC[3] = esd->GetHeader()->IsTriggerInputFired("0OM2"); //TOF two hits
  //SPD inputs
  const AliMultiplicity *mult = esd->GetMultiplicity();
  Int_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=0;
  Int_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=0;

  Int_t nInner(0), nOuter(0);
  for (Int_t i(0); i<1200; ++i) {
    Bool_t isFired(mult->TestFastOrFiredChips(i));
    if (i<400) {
      vPhiInner[i/20] += isFired;
      nInner += isFired;
    } else {
      vPhiOuter[(i-400)/20] += isFired;
      nOuter += isFired;
    }
  }
 
  Int_t fired(0);
  for (Int_t i(0); i<10; ++i) {
    for (Int_t j(0); j<2; ++j) {
      const Int_t k(2*i+j);
      fired += ((   vPhiOuter[k]    || vPhiOuter[k+1]       ||
                    vPhiOuter[k+2]      )
                && (vPhiOuter[k+20] || vPhiOuter[(k+21)%40] ||
                    vPhiOuter[(k+22)%40])
                && (vPhiInner[i]    || vPhiInner[i+1]       )
                && (vPhiInner[i+10] || vPhiInner[(i+11)%20]));
    }
  }
  //0SMB - At least one hit in SPD
  if (nOuter > 0 || nInner > 0) fTriggerInputsMC[4] = kTRUE;
  //0SM2 - Two hits on outer layer
  if (nOuter > 1) fTriggerInputsMC[5] = kTRUE;
  //0STP - Topological SPD trigger (two pairs)
  if (fired != 0) fTriggerInputsMC[6] = kTRUE;
  //0SH1 - More then 6 hits on outer layer
  if (nOuter >= 7) fTriggerInputsMC[7] = kTRUE;
  

  fGenPart->Clear("C");

  AliMCEvent *mc = MCEvent();
  if(!mc) return;

  Int_t nmc = 0;
  //loop over mc particles
  for(Int_t imc=0; imc<mc->GetNumberOfTracks(); imc++) {
    AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
    if(!mcPart) continue;

    if(mcPart->GetMother() >= 0) continue;

    TParticle *part = (TParticle*) fGenPart->ConstructedAt(nmc++);
    part->SetMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->E());
    part->SetPdgCode(mcPart->PdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles
  
  AliESDVertex *fMCVertex = (AliESDVertex*) mc->GetPrimaryVertex();
  fMCVtxPos[0] = fMCVertex->GetX();
  fMCVtxPos[1] = fMCVertex->GetY();
  fMCVtxPos[2] = fMCVertex->GetZ();

}//RunESDMC



//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcEtaC::GetMedian(Double_t *daArray) {
    // Allocate an array of the same size and sort it.
    Double_t dpSorted[4];
    for (Int_t i = 0; i < 4; ++i) {
        dpSorted[i] = daArray[i];
    }
    for (Int_t i = 3; i > 0; --i) {
        for (Int_t j = 0; j < i; ++j) {
            if (dpSorted[j] > dpSorted[j+1]) {
                Double_t dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    Double_t dMedian = 0.0;
    dMedian = (dpSorted[2] + dpSorted[1])/2.0;
    
    return dMedian;
}

bool AliAnalysisTaskUpcEtaC::CheckMeritCutWinner(int cutChoice, double oldPars[3], double newPars[3]){
 //performs "merit cut" judgement check on v0s with shared daughters, using one of three criteria.
 //if cutChoice = 4, it uses all three criteria, needed 2 of 3 'points'

  cout << "#################### Beginning of CheckMeritCutWinner()" << endl;

 bool newV0Wins = kFALSE;
 double pardiff[3] = {newPars[0]-oldPars[0],
                      newPars[1]-oldPars[1],
                      newPars[2]-oldPars[2]};
 if(cutChoice > 0 && cutChoice < 4){
  if(pardiff[cutChoice] <= 0.) newV0Wins = kTRUE;
 }
 else if(cutChoice == 4){
  int newWinCount = 0;
  for(int i=0;i<3;i++){if(pardiff[i] <= 0) newWinCount++;} //I think that original pardiff[i+1] was an error.
  if(newWinCount > 1) newV0Wins = kTRUE;  
 }
 else{};

  cout << "#################### End of CheckMeritCutWinner()" << endl;

 return newV0Wins;
}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaC::RunAODsystematics(AliAODEvent* aod)
{

  /*
  Double_t fJPsiSels[4];

  fJPsiSels[0] =   70; //min number of TPC clusters
  fJPsiSels[1] =   4; //chi2
  fJPsiSels[2] =   2; //DCAz
  fJPsiSels[3] =   1; // DCAxy 1x 

  Double_t fJPsiSelsMid[4];

  fJPsiSelsMid[0] =   70; //min number of TPC clusters
  fJPsiSelsMid[1] =   4; //chi2
  fJPsiSelsMid[2] =   2; //DCAz
  fJPsiSelsMid[3] =   1; // DCAxy 1x 
  
  Double_t fJPsiSelsLoose[4];

  fJPsiSelsLoose[0] =   60; //min number of TPC clusters
  fJPsiSelsLoose[1] =   5; //chi2
  fJPsiSelsLoose[2] =   3; //DCAz
  fJPsiSelsLoose[3] =   2; // DCAxy 2x 

  Double_t fJPsiSelsTight[4];

  fJPsiSelsTight[0] =   80; //min number of TPC clusters
  fJPsiSelsTight[1] =   3.5; //chi2
  fJPsiSelsTight[2] =   1; //DCAz
  fJPsiSelsTight[3] =   0.5; // DCAxy 0.5x 

  Int_t nGoodTracks = 0;
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4],qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0;
  Double_t fRecTPCsignal[5], fRecTPCsignalDist;
  Int_t fChannel = 0;

  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Double_t kaonMass = partKaon->Mass();
    
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();

  TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
  Double_t kStarMass = partKstar->Mass();
  Double_t kStarWidth = partKstar->Width();

  TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
  Double_t k0ShortMass = partK0short->Mass();
  Double_t k0ShortWidth = partK0short->Width();

  
for(Int_t i=0; i<5; i++){
	  //cout<<"Loose sytematics, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsLoose[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
  //Two track loop
  nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(i!=4){ if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;}
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
     
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > 1) nHighPt++;     
      		fRecTPCsignal[nLepton] = trk->GetTPCsignal();     
      		qLepton[nLepton] = trk->Charge();
      		if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				mass[nLepton] = 0;
				}
      		if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				mass[nLepton] = 1;
				}
       		nLepton++;		
  		}		
  	if(nLepton == 2){
		if(qLepton[0]*qLepton[1] < 0 && nHighPt > 0 && (mass[0]!=-1 || mass[1]!=-1)){
			vCandidate = vLepton[0]+vLepton[1];		  
  			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  			if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  			else { 
				fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  				if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
				}
			if(fChannel == -1 && vCandidate.Pt()<0.15) ((TH1D*)(fListJPsiLoose->At(i)))->Fill(vCandidate.M()); 
  			if(fChannel == 1 && vCandidate.Pt()<0.3) ((TH1D*)(fListJPsiLoose->At(i)))->Fill(vCandidate.M()); 	
			}
		}
  }
}//loose cuts

for(Int_t i=0; i<4; i++){
	  //cout<<"Tight sytematics, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsTight[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
  //Two track loop
  nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
     
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");
    
      		if(trk->Pt() > 1) nHighPt++;     
      		fRecTPCsignal[nLepton] = trk->GetTPCsignal();     
      		qLepton[nLepton] = trk->Charge();
      		if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				mass[nLepton] = 0;
				}
      		if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				mass[nLepton] = 1;
				}
       		nLepton++;		
  		}		
  	if(nLepton == 2){
		if(qLepton[0]*qLepton[1] < 0 && nHighPt > 0 && (mass[0]!=-1 || mass[1]!=-1)){
			vCandidate = vLepton[0]+vLepton[1];		  
  			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  			if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  			else { 
				fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  				if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
				}
			if(fChannel == -1 && vCandidate.Pt()<0.15) ((TH1D*)(fListJPsiTight->At(i)))->Fill(vCandidate.M()); 
  			if(fChannel == 1 && vCandidate.Pt()<0.3) ((TH1D*)(fListJPsiTight->At(i)))->Fill(vCandidate.M()); 	
			}
		}
  }
}//tight cuts

//---------------------------------------------EtaC------------------------------------------------------------------------

  Double_t fEtaCSels[4];

  fEtaCSels[0] =   50; //min number of TPC clusters
  fEtaCSels[1] =   4; //chi2
  fEtaCSels[2] =   2; //DCAz
  fEtaCSels[3] =   4; // DCAxy 1x 

  Double_t fEtaCSelsMid[4];

  fEtaCSelsMid[0] =   50; //min number of TPC clusters
  fEtaCSelsMid[1] =   4; //chi2
  fEtaCSelsMid[2] =   2; //DCAz
  fEtaCSelsMid[3] =   4; // DCAxy 1x 
  
  Double_t fEtaCSelsLoose[4];

  fEtaCSelsLoose[0] =   60; //min number of TPC clusters
  fEtaCSelsLoose[1] =   5; //chi2
  fEtaCSelsLoose[2] =   3; //DCAz
  fEtaCSelsLoose[3] =   6; // DCAxy 2x 

  Double_t fEtaCSelsTight[4];

  fEtaCSelsTight[0] =   70; //min number of TPC clusters
  fEtaCSelsTight[1] =   3.5; //chi2
  fEtaCSelsTight[2] =   1; //DCAz
  fEtaCSelsTight[3] =   2; // DCAxy 0.5x 

  nGoodTracks = 0; nLepton=0; nHighPt=0; fChannel = 0;
  Int_t nSpdHits = 0;
  Double_t trackPt[5]={0,0,0,0,0};

for(Int_t i=0; i<5; i++){
	  //cout<<"Loose systematics psi2s, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsLoose[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
 
  //Four track loop
  nGoodTracks = 0; nSpdHits = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nPion=0; nHighPt=0;
  
  if(nGoodTracks == 4){
  	  if(i!=4){ if(nSpdHits<2) continue;} 
    	  MeanPt = GetMedian(trackPt);
  	  for(Int_t k=0; k<4; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > MeanPt){   
      			fRecTPCsignal[nLepton] = trk->GetTPCsignal();      
      			qLepton[nLepton] = trk->Charge();
      			if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
					mass[nLepton] = 0;
					}
      			if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
					mass[nLepton] = 1;
					}
			nLepton++;
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}	      
    		}
	if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0) && mass[0] != -1 && mass[1] != -1){
  		vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  		if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  		else { 
			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  			if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
			}			
		if(fChannel == -1) if(vDilepton.M() > 3.0 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.15) ((TH1D*)(fListEtaCLoose->At(i)))->Fill(vCandidate.M());		
  		if(fChannel == 1) if(vDilepton.M() > 2.6 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.3) ((TH1D*)(fListEtaCLoose->At(i)))->Fill(vCandidate.M());
	}
  }   
}//loose cuts

for(Int_t i=0; i<4; i++){
	  //cout<<"Tight systematics psi2s, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsTight[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
 
  //Four track loop
  nGoodTracks = 0; nSpdHits = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
    nLepton=0; nPion=0; nHighPt=0;
  
  if(nGoodTracks == 4){
  	  if(nSpdHits<2) continue; 
    	  MeanPt = GetMedian(trackPt);
  	  for(Int_t k=0; k<4; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > MeanPt){   
      			fRecTPCsignal[nLepton] = trk->GetTPCsignal();      
      			qLepton[nLepton] = trk->Charge();
      			if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
					mass[nLepton] = 0;
					}
      			if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
					mass[nLepton] = 1;
					}
			nLepton++;
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}	      
    		}
	if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0) && mass[0] != -1 && mass[1] != -1){
  		vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  		if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  		else { 
			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  			if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
			}			
		if(fChannel == -1) if(vDilepton.M() > 3.0 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.15) ((TH1D*)(fListEtaCTight->At(i)))->Fill(vCandidate.M());		
  		if(fChannel == 1) if(vDilepton.M() > 2.6 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.3) ((TH1D*)(fListEtaCTight->At(i)))->Fill(vCandidate.M());
	}
  }   
}//Tight cuts
  */
}
