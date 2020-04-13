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
#include "TRandom.h"

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
#include "AliTriggerAnalysis.h"
#include "AliAODMCHeader.h"
#include "AliDataFile.h"

// my headers
#include "AliAnalysisTaskUpcPsi2s.h"

ClassImp(AliAnalysisTaskUpcPsi2s);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s() 
  : AliAnalysisTaskSE(),fType(0),fTracking(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFmask(0),fIsPhysicsSelected(kFALSE),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fV0Adecision(0),fV0Cdecision(0),fADAdecision(0),fADCdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),fFOFiredChips(0),fIR1Map(0),fIR2Map(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),fGenPart(0),
    fEveTree(0),fPt(0), fY(0), fM(0), fDiLeptonM(0), fDiLeptonPt(0), fPIDsigma(0), fChannel(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fHistCcup8TriggersPerRun(0),fHistCcup9TriggersPerRun(0),fHistCcup10TriggersPerRun(0),fHistCcup11TriggersPerRun(0),fHistCcup12TriggersPerRun(0),
    fHistCcup25TriggersPerRun(0),fHistCcup26TriggersPerRun(0),fHistCcup27TriggersPerRun(0),
    fHistCcup29TriggersPerRun(0),fHistCcup30TriggersPerRun(0),fHistCcup31TriggersPerRun(0),
    fHistCtrueTriggersPerRun(0),
    fListHist(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),fHistDiLeptonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0),
    fListSystematics(0),fListJPsiLoose(0),fListJPsiTight(0),
    fSPDfile(0),hBCmod4(0),hSPDeff(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fTracking(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFmask(0),fIsPhysicsSelected(kFALSE),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fV0Adecision(0),fV0Cdecision(0),fADAdecision(0),fADCdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),fFOFiredChips(0),fIR1Map(0),fIR2Map(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),fGenPart(0),
    fEveTree(0),fPt(0), fY(0), fM(0), fDiLeptonM(0), fDiLeptonPt(0), fPIDsigma(0), fChannel(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0), fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fHistCcup8TriggersPerRun(0),fHistCcup9TriggersPerRun(0),fHistCcup10TriggersPerRun(0),fHistCcup11TriggersPerRun(0),fHistCcup12TriggersPerRun(0),
    fHistCcup25TriggersPerRun(0),fHistCcup26TriggersPerRun(0),fHistCcup27TriggersPerRun(0),
    fHistCcup29TriggersPerRun(0),fHistCcup30TriggersPerRun(0),fHistCcup31TriggersPerRun(0),
    fHistCtrueTriggersPerRun(0),
    fListHist(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),fHistDiLeptonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0),
    fListSystematics(0),fListJPsiLoose(0),fListJPsiTight(0),
    fSPDfile(0),hBCmod4(0),hSPDeff(0)

{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

}//AliAnalysisTaskUpcPsi2s

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::Init()
{
  
  for(Int_t i=0; i<ntrgMB; i++) {
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
	
	fZNATDCm[i] = -666;
        fZNCTDCm[i] = -666;
	fZPATDCm[i] = -666;
        fZPCTDCm[i] = -666;
	
	}
  for(Int_t i=0; i<3; i++){
  	fVtxPos[i] = -666; 
	fMCVtxPos[i] = -666;
	fVtxErr[i] = -666;
	fKfVtxPos[i] = -666;
	fSpdVtxPos[i] = -666;
	}

}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::~AliAnalysisTaskUpcPsi2s() 
{
  // Destructor
  if(fJPsiTree){
     delete fJPsiTree;
     fJPsiTree = 0x0;
  }
  if(fPsi2sTree){
     delete fPsi2sTree;
     fPsi2sTree = 0x0;
  }
  if(fListTrig){
     delete fListTrig;
     fListTrig = 0x0;
  }
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
  }

}//~AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserCreateOutputObjects()
{
  //PID response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  //input file
  fDataFilnam = new TObjString();
  fDataFilnam->SetString("");

    //tracks
  fJPsiAODTracks = new TClonesArray("AliAODTrack", 1000);
  fJPsiESDTracks = new TClonesArray("AliESDtrack", 1000);
  fPsi2sAODTracks = new TClonesArray("AliAODTrack", 1000);
  fPsi2sESDTracks = new TClonesArray("AliESDtrack", 1000);
  fGenPart = new TClonesArray("TParticle", 1000);

  //output tree with JPsi candidate events
  fJPsiTree = new TTree("fJPsiTree", "fJPsiTree");
  fJPsiTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
//  fJPsiTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
//  fJPsiTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
//  fJPsiTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fJPsiTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrgMB));
  fJPsiTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
//  fJPsiTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
//  fJPsiTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
//  fJPsiTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fJPsiTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
//  fJPsiTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
//  fJPsiTree ->Branch("fTOFmask", &fTOFmask);
  
//  fJPsiTree ->Branch("fIsPhysicsSelected", &fIsPhysicsSelected, "fIsPhysicsSelected/O");
  
  fJPsiTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[2]/F");
  fJPsiTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[2]/F");
  fJPsiTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[2]/F");
//  fJPsiTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[2]/F");
  fJPsiTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[2]/F");
  
//  fJPsiTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[2]/F");
//  fJPsiTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[2]/F");
//  fJPsiTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[2]/F");
//  fJPsiTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[2]/F");
  fJPsiTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[2]/F");
  
//  fJPsiTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[2]/O");
  
  fJPsiTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/F");
//  fJPsiTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/F");
//  fJPsiTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/F");
//  fJPsiTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/F");
  
//  fJPsiTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/F");
//  fJPsiTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/F");
  
  fJPsiTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/F");
  fJPsiTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/F");
  fJPsiTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/F");
  fJPsiTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/F");
  fJPsiTree ->Branch("fZNATDCm", &fZNATDCm[0], "fZNATDCm[4]/F");
  fJPsiTree ->Branch("fZNCTDCm", &fZNCTDCm[0], "fZNCTDCm[4]/F");
//  fJPsiTree ->Branch("fZPATDCm", &fZPATDCm[0], "fZPATDCm[4]/F");
//  fJPsiTree ->Branch("fZPCTDCm", &fZPCTDCm[0], "fZPCTDCm[4]/F");
  fJPsiTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fJPsiTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I"); 
  fJPsiTree ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fJPsiTree ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");  
//  fJPsiTree ->Branch("fDataFilnam", &fDataFilnam);
//  fJPsiTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
//  fJPsiTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");		       
  if( fType == 0 ) {
    fJPsiTree ->Branch("fJPsiESDTracks", &fJPsiESDTracks);
  }
  if( fType == 1 ) {
    fJPsiTree ->Branch("fJPsiAODTracks", &fJPsiAODTracks);
  }
  if(isMC) {
    fJPsiTree ->Branch("fGenPart", &fGenPart);
    fJPsiTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O", ntrgMB));
//    fJPsiTree ->Branch("fMCVtxPos", &fMCVtxPos[0], "fMCVtxPos[3]/F");
//    fJPsiTree ->Branch("fFOFiredChips", &fFOFiredChips);
  }

 
 //output tree with Psi2s candidate events
  fPsi2sTree = new TTree("fPsi2sTree", "fPsi2sTree");
  fPsi2sTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
//  fPsi2sTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
//  fPsi2sTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
//  fPsi2sTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fPsi2sTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrgMB));
  fPsi2sTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
//  fPsi2sTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
//  fPsi2sTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
//  fPsi2sTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fPsi2sTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
//  fPsi2sTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
//  fPsi2sTree ->Branch("fTOFmask", &fTOFmask);
  
//  fPsi2sTree ->Branch("fIsPhysicsSelected", &fIsPhysicsSelected, "fIsPhysicsSelected/O");
  
  fPsi2sTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[4]/F");
  fPsi2sTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[4]/F");
  fPsi2sTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[4]/F");
//  fPsi2sTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[4]/F");
  fPsi2sTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[4]/F");
  
//  fPsi2sTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[4]/F");
//  fPsi2sTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[4]/F");
//  fPsi2sTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[4]/F");
//  fPsi2sTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[4]/F");
  fPsi2sTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[4]/F");
  
//  fPsi2sTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[4]/O");
  
  fPsi2sTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/F");
//  fPsi2sTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/F");
//  fPsi2sTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/F");
//  fPsi2sTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/F");
  
//  fPsi2sTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/F");
//  fPsi2sTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/F");
  
  fPsi2sTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/F");
  fPsi2sTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/F");
  fPsi2sTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/F");
  fPsi2sTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/F");
  fPsi2sTree ->Branch("fZNATDCm", &fZNATDCm[0], "fZNATDCm[4]/F");
  fPsi2sTree ->Branch("fZNCTDCm", &fZNCTDCm[0], "fZNCTDCm[4]/F");
  fPsi2sTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fPsi2sTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I"); 
  fPsi2sTree ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fPsi2sTree ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");  
//  fPsi2sTree ->Branch("fDataFilnam", &fDataFilnam);
//  fPsi2sTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
//  fPsi2sTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");	       
  if( fType == 0 ) {
    fPsi2sTree ->Branch("fPsi2sESDTracks", &fPsi2sESDTracks);
  }
  if( fType == 1 ) {
    fPsi2sTree ->Branch("fPsi2sAODTracks", &fPsi2sAODTracks);
  }
  if(isMC) {
    fPsi2sTree ->Branch("fGenPart", &fGenPart);
    fPsi2sTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O", ntrgMB));
//    fPsi2sTree ->Branch("fMCVtxPos", &fMCVtxPos[0], "fMCVtxPos[3]/F");
//    fPsi2sTree ->Branch("fFOFiredChips", &fFOFiredChips);
  }
  
  fListTrig = new TList();
  fListTrig ->SetOwner();
  
  fHistCcup4TriggersPerRun = new TH1D("fHistCcup4TriggersPerRun", "fHistCcup4TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup4TriggersPerRun);
  
  fHistCcup7TriggersPerRun = new TH1D("fHistCcup7TriggersPerRun", "fHistCcup7TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup7TriggersPerRun);
    
  fHistCcup2TriggersPerRun = new TH1D("fHistCcup2TriggersPerRun", "fHistCcup2TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup2TriggersPerRun);
  
  fHistCint1TriggersPerRun = new TH1D("fHistCint1TriggersPerRun", "fHistCint1TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCint1TriggersPerRun);
  
  fHistCint6TriggersPerRun = new TH1D("fHistCint6TriggersPerRun", "fHistCint6TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCint6TriggersPerRun);
  
  fHistC0tvxAndCint1TriggersPerRun = new TH1D("fHistC0tvxAndCint1TriggersPerRun", "fHistC0tvxAndCint1TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistC0tvxAndCint1TriggersPerRun);
  
  fHistZedTriggersPerRun = new TH1D("fHistZedTriggersPerRun", "fHistZedTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistZedTriggersPerRun);

  fHistCvlnTriggersPerRun = new TH1D("fHistCvlnTriggersPerRun", "fHistCvlnTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCvlnTriggersPerRun);
  
  fHistMBTriggersPerRun = new TH1D("fHistMBTriggersPerRun", "fHistMBTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistMBTriggersPerRun);
  
  fHistCentralTriggersPerRun = new TH1D("fHistCentralTriggersPerRun", "fHistCentralTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCentralTriggersPerRun);
  
  fHistSemiCentralTriggersPerRun = new TH1D("fHistSemiCentralTriggersPerRun", "fHistSemiCentralTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistSemiCentralTriggersPerRun);
  
  fHistCcup8TriggersPerRun = new TH1D("fHistCcup8TriggersPerRun", "fHistCcup8TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup8TriggersPerRun);
  
  fHistCcup9TriggersPerRun = new TH1D("fHistCcup9TriggersPerRun", "fHistCcup9TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup9TriggersPerRun);
  
  fHistCcup10TriggersPerRun = new TH1D("fHistCcup10TriggersPerRun", "fHistCcup10TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup10TriggersPerRun);
  
  fHistCcup11TriggersPerRun = new TH1D("fHistCcup11TriggersPerRun", "fHistCcup11TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup11TriggersPerRun);
  
  fHistCcup12TriggersPerRun = new TH1D("fHistCcup12TriggersPerRun", "fHistCcup12TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup12TriggersPerRun);
  
  fHistCcup25TriggersPerRun = new TH1D("fHistCcup25TriggersPerRun", "fHistCcup25TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup25TriggersPerRun);
  
  fHistCcup26TriggersPerRun = new TH1D("fHistCcup26TriggersPerRun", "fHistCcup26TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup26TriggersPerRun);
  
  fHistCcup27TriggersPerRun = new TH1D("fHistCcup27TriggersPerRun", "fHistCcup27TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup27TriggersPerRun);
  
  fHistCcup29TriggersPerRun = new TH1D("fHistCcup29TriggersPerRun", "fHistCcup29TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup29TriggersPerRun);
  
  fHistCcup30TriggersPerRun = new TH1D("fHistCcup30TriggersPerRun", "fHistCcup30TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup30TriggersPerRun);
  
  fHistCcup31TriggersPerRun = new TH1D("fHistCcup31TriggersPerRun", "fHistCcup31TriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCcup31TriggersPerRun);

  fHistCtrueTriggersPerRun = new TH1D("fHistCtrueTriggersPerRun", "fHistCtrueTriggersPerRun", 60000, 240000.5, 300000.5);
  fListTrig->Add(fHistCtrueTriggersPerRun);
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
  TString CutNameJPsi[13] = {"Analyzed","Triggered","Vertex cut","V0 decision","Neutron ZDC cut","Two good tracks",
  				"Like sign","Oposite sign","One p_{T}>1", "Both p_{T}>1","PID","Dimuom","Dielectron"};
  fHistNeventsJPsi = new TH1D("fHistNeventsJPsi","fHistNeventsPsi2s",13,0.5,13.5);
  for (Int_t i = 0; i<13; i++) fHistNeventsJPsi->GetXaxis()->SetBinLabel(i+1,CutNameJPsi[i].Data());
  fListHist->Add(fHistNeventsJPsi);
  
  fHistTPCsignalJPsi = new TH2D("fHistTPCsignalJPsi","fHistTPCsignalJPsi",240,0,120,240,0,120);
  fListHist->Add(fHistTPCsignalJPsi);
  
  fHistDiLeptonPtJPsi = new TH2D("fHistDiLeptonPtJPsi","fHistDiLeptonPtJPsi",350,0,3.5,350,0,3.5);
  fListHist->Add(fHistDiLeptonPtJPsi);

  fHistDiElectronMass = new TH1D("fHistDiElectronMass","Invariant mass of J/#psi candidates",100,2,5);
  fHistDiElectronMass->GetXaxis()->SetTitle("Invariant mass(e^{+}e^{-}) (GeV/c)");
  fListHist->Add(fHistDiElectronMass);
  
  fHistDiMuonMass = new TH1D("fHistDiMuonMass","Invariant mass of J/#psi candidates",100,2,5);
  fHistDiMuonMass->GetXaxis()->SetTitle("Invariant mass(#mu^{+}#mu^{-}) (GeV/c)");
  fListHist->Add(fHistDiMuonMass);
  
  fHistDiLeptonMass = new TH1D("fHistDiLeptonMass","Invariant mass of J/#psi candidates",130,2.1,6.0);
  fHistDiLeptonMass->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}) (GeV/c)");
  fListHist->Add(fHistDiLeptonMass);

  TString CutNamePsi2s[14] = {"Analyzed","Triggered","Vertex cut","V0 decision","Neutron ZDC cut","Four good tracks",
  				"DiLepton - DiPion","Like sign leptons","Like sign pions","Like sign both","Oposite sign","PID","Dimuom","Dielectron"};

  fHistNeventsPsi2s = new TH1D("fHistNeventsPsi2s","fHistNeventsPsi2s",14,0.5,14.5);
  for (Int_t i = 0; i<14; i++) fHistNeventsPsi2s->GetXaxis()->SetBinLabel(i+1,CutNamePsi2s[i].Data());
  fListHist->Add(fHistNeventsPsi2s);

  fHistPsi2sMassVsPt = new TH2D("fHistPsi2sMassVsPt","Mass vs p_{T} of #psi(2s) candidates",100,3,6,50,0,5);
  fHistPsi2sMassVsPt->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fHistPsi2sMassVsPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fListHist->Add(fHistPsi2sMassVsPt);
  
  fHistPsi2sMassCoherent = new TH1D("fHistPsi2sMassAllCoherent","Invariant mass of coherent #psi(2s) candidates",50,2.5,5.5);
  fHistPsi2sMassCoherent->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fListHist->Add(fHistPsi2sMassCoherent);
  
  fEveTree = new TTree("fEveTree", "fEveTree");
  fEveTree ->Branch("fPt", &fPt, "fPt/F");
  fEveTree ->Branch("fY", &fY, "fY/F");
  fEveTree ->Branch("fM", &fM, "fM/F");
  fEveTree ->Branch("fChannel", &fChannel, "fChannel/I");
  fEveTree ->Branch("fDiLeptonM", &fDiLeptonM, "fDiLeptonM/F");
  fEveTree ->Branch("fDiLeptonPt", &fDiLeptonPt, "fDiLeptonPt/F");
  fEveTree ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fEveTree ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fEveTree ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fEveTree ->Branch("fDataFilnam", &fDataFilnam);
  fEveTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L"); 
  fListHist->Add(fEveTree);
  
  
  fListSystematics = new TList();
  fListSystematics->SetOwner();
  fListSystematics->SetName("fListSystematics");
  fListHist->Add(fListSystematics);
  InitSystematics();
  
  fSPDfile = AliDataFile::OpenOADB("PWGUD/UPC/SPDFOEfficiency_run245067.root");
  fSPDfile->Print();
  fSPDfile->Map();
  hSPDeff = (TH2D*) fSPDfile->Get("hEff");
  hSPDeff->SetDirectory(0);
  TH2D *hBCmod4_2D = (TH2D*) fSPDfile->Get("hCounts");
  hBCmod4 = hBCmod4_2D->ProjectionY();
  fSPDfile->Close();
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);
  PostData(3, fListTrig);
  PostData(4, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::InitSystematics()
{ 

fListJPsiLoose = new TList();
fListJPsiLoose->SetOwner();
fListJPsiLoose->SetName("JPsiLoose");
fListSystematics->Add(fListJPsiLoose);

TH1D *fHistJPsiNClusLoose = new TH1D("JPsiNClusLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiLoose->Add(fHistJPsiNClusLoose);

TH1D *fHistJPsiChi2Loose = new TH1D("JPsiChi2Loose","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiLoose->Add(fHistJPsiChi2Loose);

TH1D *fHistJPsiDCAzLoose = new TH1D("JPsiDCAzLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiLoose->Add(fHistJPsiDCAzLoose);

TH1D *fHistJPsiDCAxyLoose = new TH1D("JPsiDCAxyLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiLoose->Add(fHistJPsiDCAxyLoose);

TH1D *fHistJPsiITShitsLoose = new TH1D("JPsiITShitsLoose","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiLoose->Add(fHistJPsiITShitsLoose);


fListJPsiTight = new TList();
fListJPsiTight->SetOwner();
fListJPsiTight->SetName("JPsiTight");
fListSystematics->Add(fListJPsiTight);

TH1D *fHistJPsiNClusTight = new TH1D("JPsiNClusTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiTight->Add(fHistJPsiNClusTight);

TH1D *fHistJPsiChi2Tight = new TH1D("JPsiChi2Tight","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiTight->Add(fHistJPsiChi2Tight);

TH1D *fHistJPsiDCAzTight = new TH1D("JPsiDCAzTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiTight->Add(fHistJPsiDCAzTight);

TH1D *fHistJPsiDCAxyTight = new TH1D("JPsiDCAxyTight","Invariant mass of J/#psi candidates",130,2.1,6.0);
fListJPsiTight->Add(fHistJPsiDCAxyTight);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ){
    	//RunESDtrig(); 
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
void AliAnalysisTaskUpcPsi2s::RunAODtrig()
{

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
  
  if(trigger.Contains("CCUP8-B")) fHistCcup8TriggersPerRun->Fill(fRunNum); //CCUP8 triggers
  if(trigger.Contains("CCUP9-B")) fHistCcup9TriggersPerRun->Fill(fRunNum); //CCUP9 triggers
  if(trigger.Contains("CCUP10-B")) fHistCcup10TriggersPerRun->Fill(fRunNum); //CCUP10 triggers
  if(trigger.Contains("CCUP11-B")) fHistCcup11TriggersPerRun->Fill(fRunNum); //CCUP11 triggers
  if(trigger.Contains("CCUP12-B")) fHistCcup12TriggersPerRun->Fill(fRunNum); //CCUP12 triggers
  
  if(trigger.Contains("CCUP25-B")) fHistCcup25TriggersPerRun->Fill(fRunNum); //CCUP25 triggers
  if(trigger.Contains("CCUP26-B")) fHistCcup26TriggersPerRun->Fill(fRunNum); //CCUP26 triggers
  if(trigger.Contains("CCUP27-B")) fHistCcup27TriggersPerRun->Fill(fRunNum); //CCUP27 triggers
  
  if(fRunNum>=295881){
    if(trigger.Contains("CCUP30-B-SPD2-CENTNOTRD")) fHistCcup30TriggersPerRun->Fill(fRunNum); //CCUP30 triggers
    if(trigger.Contains("CCUP31-B-SPD2-CENTNOTRD")) fHistCcup31TriggersPerRun->Fill(fRunNum); //CCUP31 triggers
    if(fRunNum>=296594 && trigger.Contains("CCUP29-U-SPD2-CENTNOTRD")) fHistCcup29TriggersPerRun->Fill(fRunNum); //CCUP29 triggers
    if(fRunNum<296594  && trigger.Contains("CCUP29-B-SPD2-CENTNOTRD")) fHistCcup29TriggersPerRun->Fill(fRunNum); //CCUP29 triggers
	}
  else{ 
	if(trigger.Contains("CCUP29-B-NOPF-CENTNOTRD")) fHistCcup29TriggersPerRun->Fill(fRunNum); //CCUP29 triggers
	if(trigger.Contains("CCUP30-B-NOPF-CENTNOTRD")) fHistCcup30TriggersPerRun->Fill(fRunNum); //CCUP30 triggers
	if(trigger.Contains("CCUP31-B-NOPF-CENTNOTRD")) fHistCcup31TriggersPerRun->Fill(fRunNum); //CCUP31 triggers
	}

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
    
PostData(3, fListTrig);

}
//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunAODhist()
{

  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  Double_t muonMass = partMuon->Mass();
  
  TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
  Double_t electronMass = partElectron->Mass();
  
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();

  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;
  
 // cout<<"Event number: "<<((TTree*) GetInputData(0))->GetTree()->GetReadEntry()<<endl;

  fHistNeventsJPsi->Fill(1);
  fHistNeventsPsi2s->Fill(1);

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if(!isMC && !trigger.Contains("CCUP") ) return;
  
  fHistNeventsJPsi->Fill(2);
  fHistNeventsPsi2s->Fill(2);
  
  AliAODZDC *fZDCdata = aod->GetZDCData();
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  
  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsJPsi->Fill(3);
  fHistNeventsPsi2s->Fill(3);

  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;
  
  fHistNeventsJPsi->Fill(4);
  fHistNeventsPsi2s->Fill(4);

  if( fZNAenergy > 8200 || fZNCenergy > 8200) return;
  
  fHistNeventsJPsi->Fill(5);
  fHistNeventsPsi2s->Fill(5); 
  
  //Systematics - cut variation
  if(fRunSystematics) RunAODsystematics(aod);

  //Two tracks loop
  Int_t nGoodTracks = 0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4], qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0, nSpdHits=0;
  Double_t fRecTPCsignal[5], fRecTPCsignalDist;
  Int_t fChannel = 0;
  Int_t mass[3]={-1,-1,-1};
  Double_t TrackPt[5]={0,0,0,0,0};
  Double_t MeanPt = -1;
  
   
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
     
      TrackIndex[nGoodTracks] = itr;
      TrackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
  
  nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;
  
  if(nGoodTracks == 4 && nSpdHits>1){
    	  MeanPt = GetMedian(TrackPt);
  	  fHistNeventsPsi2s->Fill(6);
  	  for(Int_t i=0; i<4; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));
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
				      
      		if(nLepton > 2 || nPion > 2) break;
    		}
	if((nLepton == 2) && (nPion == 2)){
		fHistNeventsPsi2s->Fill(7);
		if(qLepton[0]*qLepton[1] > 0) fHistNeventsPsi2s->Fill(8);
		if(qPion[0]*qPion[1] > 0) fHistNeventsPsi2s->Fill(9);
		if((qLepton[0]*qLepton[1] > 0) && (qPion[0]*qPion[1] > 0)) fHistNeventsPsi2s->Fill(10);
		if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0)){
			fHistNeventsPsi2s->Fill(11);
	 		if(mass[0] != -1 && mass[1] != -1) {
				fHistNeventsPsi2s->Fill(12); 
  				vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  				vDilepton = vLepton[0]+vLepton[1];
				fHistPsi2sMassVsPt->Fill(vCandidate.M(),vCandidate.Pt());
				fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  				if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  				else { 
					fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  					if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
					}
				
				if(fChannel == -1) {
					fHistNeventsPsi2s->Fill(13);
					if(vDilepton.M() > 3.0 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.15) fHistPsi2sMassCoherent->Fill(vCandidate.M());
					}	
  				if(fChannel == 1){ 
					fHistNeventsPsi2s->Fill(14);
					if(vDilepton.M() > 2.6 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.3) fHistPsi2sMassCoherent->Fill(vCandidate.M());
					}
				}
			}
		}
  }
  
  nGoodTracks = 0;
  //Two track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 70)continue;
      if(trk->Chi2perNDF() > 4)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
  
   nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;

  if(nGoodTracks == 2){
  	  fHistNeventsJPsi->Fill(6);
  	  for(Int_t i=0; i<2; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));		
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
	 	if(qLepton[0]*qLepton[1] > 0) fHistNeventsJPsi->Fill(7);
		if(qLepton[0]*qLepton[1] < 0){
			fHistNeventsJPsi->Fill(8);
			if(nHighPt > 0){
				fHistNeventsJPsi->Fill(9);
				fHistTPCsignalJPsi->Fill(fRecTPCsignal[0],fRecTPCsignal[1]);
				if(nHighPt == 2) fHistNeventsJPsi->Fill(10);
				if(mass[0] != -1 && mass[1] != -1) {
					fHistNeventsJPsi->Fill(11);
					vCandidate = vLepton[0]+vLepton[1];
					fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  					if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  					else { 
						fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  						if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
						}
					if( vCandidate.M() > 2.8 && vCandidate.M() < 3.2) fHistDiLeptonPtJPsi->Fill(vLepton[0].Pt(),vLepton[1].Pt());
					if(fChannel == -1) {
						fHistDiMuonMass->Fill(vCandidate.M());
						if(vCandidate.Pt()<0.15)fHistDiLeptonMass->Fill(vCandidate.M());
						fHistNeventsJPsi->Fill(12);
						}
  					if(fChannel == 1) {
						fHistDiElectronMass->Fill(vCandidate.M());
						if(vCandidate.Pt()<0.3)fHistDiLeptonMass->Fill(vCandidate.M());
						fHistNeventsJPsi->Fill(13);
   						}
					}
				}
			}
		}
  }
  
 
  PostData(4, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunAODtree()
{
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
  if(fTracking == 0){ 
  	fTrigger[0]  = trigger.Contains("CCUP4-B"); // *0VBA *0VBC 0SM2 0OMU
  	fTrigger[1]  = trigger.Contains("CCUP2-B"); // *0VBA *0VBC 0SM2 0OM2
  	fTrigger[2]  = trigger.Contains("CCUP7-B"); // *0VBA *0VBC 0STP 0OMU
  	fTrigger[3]  = trigger.Contains("CINT1-B"); //  0VBA || 0VBC || 0SMB
  	fTrigger[4]  = trigger.Contains("CCUP8-B"); //*0VBA *0VBC *0UBA *0UBC 0STP 0OMU
  	fTrigger[5]  = trigger.Contains("CCUP9-B"); //*0VBA *0VBC *0UBA *0UBC 0STP
  	fTrigger[6]  = trigger.Contains("CCUP10-B"); //*0VBA *0VBC *0UBA *0UBC 0SH1
  	fTrigger[7]  = trigger.Contains("CCUP11-B"); //*0UBA *0UBC 0STP 0OMU
  	fTrigger[8]  = trigger.Contains("CCUP12-B"); //*0UBA *0UBC 0STP
  	fTrigger[9]  = trigger.Contains("CTRUE-B"); //Unbiased trigger
	fTrigger[10]  = trigger.Contains("CCUP25-B");//*0VBA *0VBC 0STG 0OM2
	fTrigger[11]  = trigger.Contains("CCUP26-B");//*0VBA *0VBC 0SH1 
	fTrigger[12]  = trigger.Contains("CCUP27-B");//*0VBA *0VBC 0STG 
	fTrigger[13]  = trigger.Contains("CCUP29-U-SPD2-CENTNOTRD");//*0VBA *0VBC *0UBA *0UBC 0STG
	fTrigger[14]  = trigger.Contains("CCUP29-B-SPD2-CENTNOTRD");//*0VBA *0VBC *0UBA *0UBC 0STG
	fTrigger[15]  = trigger.Contains("CCUP29-B-NOPF-CENTNOTRD");//*0VBA *0VBC *0UBA *0UBC 0STG
	fTrigger[16]  = trigger.Contains("CCUP30-B-NOPF-CENTNOTRD");//*0VBA *0VBC *0UBA *0UBC 0STG 0OM2
	fTrigger[17]  = trigger.Contains("CCUP30-B-SPD2-CENTNOTRD");//*0VBA *0VBC *0UBA *0UBC 0STG 0OM2
	fTrigger[18]  = trigger.Contains("CCUP31-B-NOPF-CENTNOTRD");// *0VBA *0VBC *0UBA *0UBC 0STG 0OMU
	fTrigger[19]  = trigger.Contains("CCUP31-B-SPD2-CENTNOTRD");// *0VBA *0VBC *0UBA *0UBC 0STG 0OMU
	}
  if(fTracking == 1){ 
  	fTrigger[0] = trigger.Contains("CCUP14-B"); 
  	fTrigger[1] = trigger.Contains("CCUP15-B"); 
  	fTrigger[2] = trigger.Contains("CCUP20-B"); 
  	fTrigger[3] = trigger.Contains("CCUP21-B"); 
  	fTrigger[4] = trigger.Contains("CCUP22-B"); 
  	fTrigger[5] = trigger.Contains("CCUP23-B"); 
  	fTrigger[6] = trigger.Contains("CCUP24-B");
	fTrigger[7] = trigger.Contains("CCUP26-B");//*0VBA *0VBC 0SH1 
	fTrigger[8] = trigger.Contains("CCUP27-B");//*0VBA *0VBC 0STG
	fTrigger[9]  = trigger.Contains("CCUP29-B");//*0VBA *0VBC *0UBA *0UBC 0STG
	fTrigger[10]  = trigger.Contains("CCUP30-B");//*0VBA *0VBC *0UBA *0UBC 0STG 0OM2
	fTrigger[11]  = trigger.Contains("CCUP31-B");// *0VBA *0VBC *0UBA *0UBC 0STG 0OMU
	} 
   if(fTracking == 8){ 
  	fTrigger[0]  = trigger.Contains("CMUP10-B");	// *0VBA *0UBA *0UBC 0MSL			
  	fTrigger[1]  = trigger.Contains("CMUP11-B");	// !0VBA & !0UBA & !0UBC & 0MUL 		
  	fTrigger[2]  = trigger.Contains("CMUP12-B");	// !0VBA & !0UBA & !0UBC & 0MSL & 0SMB
  	fTrigger[3]  = trigger.Contains("CMUP14-B");   // 0MSL & !0VBA & !0UBA
  	fTrigger[4]  = trigger.Contains("CMUP15-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL
  	fTrigger[5]  = trigger.Contains("CMUP16-B");   // 0MSL *0VBA *0UBA *0UGC *0VGA
  	fTrigger[6]  = trigger.Contains("CMUP17-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL *0UGC *0VGA
  	fTrigger[7]  = trigger.Contains("CMUP21-B");   // *0VBA *0UBA *0VBC 0SH1 *0SH2 *0UGC *0VGA
  	fTrigger[8]  = trigger.Contains("CMUP22-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MSL 0SMB
  	fTrigger[9]  = trigger.Contains("CMUP23-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MUL
	}
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrgMB; i++) {
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
  
  //Past-future protection maps
  fIR1Map = aod->GetHeader()->GetIRInt1InteractionMap();
  fIR2Map = aod->GetHeader()->GetIRInt2InteractionMap();

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
  for (Int_t i=0;i<4;i++){
  	fZNATDCm[i] = fZDCdata->GetZNATDCm(i);
  	fZNCTDCm[i] = fZDCdata->GetZNCTDCm(i);
	fZPATDCm[i] = fZDCdata->GetZPATDCm(i);
  	fZPCTDCm[i] = fZDCdata->GetZPCTDCm(i);
	}
  
  fNLooseTracks = 0;
  
  //Track loop - loose cuts
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(fTracking == 0){
      if(!(trk->TestFilterBit(1<<0))) continue;
      //if(!trk->HasPointOnITSLayer(0)||!trk->HasPointOnITSLayer(1)) continue;
      
      fNLooseTracks++;
      }
    if(fTracking == 1){
      if(!(trk->TestFilterBit(1<<1))) continue;
      
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      //fNLooseTracks++;
      }  
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  //Two track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    
    if(fTracking == 0){
      if(!trk->TestFilterBit(1<<5)) continue;
    
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!trk->TestFilterBit(1<<1)) continue;
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 8){  
      if(!trk->IsMuonTrack())continue;
      if( trk->GetRAtAbsorberEnd() < 17.5 || trk->GetRAtAbsorberEnd() > 89.5 ) continue;
      if( trk->Eta() < -4.0 || trk->Eta() > -2.5 ) continue;
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
				  
      if(nGoodTracks > 2) break;
  }//Track loop
     	
  fJPsiAODTracks->Clear("C");
  if(nGoodTracks == 2){

   	  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
	  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  	  Double_t muonMass = partMuon->Mass();  
          TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
          Double_t electronMass = partElectron->Mass();  
  	  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  	  Double_t pionMass = partPion->Mass();
  
  	  Double_t KFcov[21];
  	  Double_t KFpar[6];
	  Double_t KFmass = pionMass;
	  Double_t fRecTPCsignal;
  	  AliKFParticle *KFpart[2];
  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(aod->GetMagneticField()); 
  
  	  for(Int_t i=0; i<2; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));
                if(!trk) AliFatal("Not a standard AOD");
		
		if(fAODVertex->HasDaughter(trk) && trk->GetUsedForVtxFit())fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;

		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      		if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      		delete trk_clone;
				
		new((*fJPsiAODTracks)[i]) AliAODTrack(*trk); 
		((AliAODTrack*)((*fJPsiAODTracks)[i]))->SetDCA(dca[0],dca[1]);//to get DCAxy trk->DCA(); to get DCAz trk->ZAtDCA();
		
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
		
		if(trk->Pt() > 1){   
      			fRecTPCsignal = trk->GetTPCsignal();      
      			if(fRecTPCsignal > 40 && fRecTPCsignal < 70) KFmass = muonMass;
      			if(fRecTPCsignal > 70 && fRecTPCsignal < 100)KFmass = electronMass;
			}
		else KFmass = pionMass;
		
		KFpart[i] = new AliKFParticle();
    		KFpart[i]->SetField(aod->GetMagneticField());
    		KFpart[i]->AliKFParticleBase::Initialize(KFpar,KFcov,(Int_t) trk->Charge(), KFmass);
		KFvtx->AddDaughter(*KFpart[i]); 		
  		}
  fKfVtxPos[0]= KFvtx->GetX();
  fKfVtxPos[1]= KFvtx->GetY();
  fKfVtxPos[2]= KFvtx->GetZ();
  for(UInt_t i=0; i<2; i++)delete KFpart[i];
  delete KFvtx; 

  if(!isMC) fJPsiTree ->Fill();
  }
  
   nGoodTracks = 0;
   //Four track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    
    if(fTracking == 0){
      if(!(trk->TestFilterBit(1<<4))) continue;
  
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!(trk->TestFilterBit(1<<1))) continue;
       
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 8){  
      if(!trk->IsMuonTrack())continue;
      if( trk->GetRAtAbsorberEnd() < 17.5 || trk->GetRAtAbsorberEnd() > 89.5 ) continue;
      if( trk->Eta() < -4.0 || trk->Eta() > -2.5 ) continue;
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
      
  fPsi2sAODTracks->Clear("C");  
  if(nGoodTracks == 4){

  	  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
	  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  	  Double_t muonMass = partMuon->Mass();  
          TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
          Double_t electronMass = partElectron->Mass();  
  	  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  	  Double_t pionMass = partPion->Mass();
  
  	  Double_t KFcov[21];
  	  Double_t KFpar[6];
	  Double_t KFmass = pionMass;
	  Double_t fRecTPCsignal;
  	  AliKFParticle *KFpart[4];
  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(aod->GetMagneticField()); 
	  	  
  	  for(Int_t i=0; i<4; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));
                if(!trk) AliFatal("Not a standard AOD");
		
		if(fAODVertex->HasDaughter(trk) && trk->GetUsedForVtxFit())fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;

		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      		if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      		delete trk_clone;
		
		new((*fPsi2sAODTracks)[i]) AliAODTrack(*trk);
		((AliAODTrack*)((*fPsi2sAODTracks)[i]))->SetDCA(dca[0],dca[1]);//to get DCAxy trk->DCA(); to get DCAz trk->ZAtDCA();
		
		
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
		
		if(trk->Pt() > 1){   
      			fRecTPCsignal = trk->GetTPCsignal();      
      			if(fRecTPCsignal > 40 && fRecTPCsignal < 70) KFmass = muonMass;
      			if(fRecTPCsignal > 70 && fRecTPCsignal < 100)KFmass = electronMass;
			}
		else KFmass = pionMass;
		
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
  if(!isMC) fPsi2sTree ->Fill();
  }
  
  if(isMC){
  	fJPsiTree ->Fill();
	fPsi2sTree ->Fill();
  }
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);

}//RunAOD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunAODMC(AliAODEvent *aod)
{

  for(Int_t i=0; i<11; i++) fTriggerInputsMC[i] = kFALSE;
  
  UShort_t fTriggerAD = aod->GetADData()->GetTriggerBits();
  UShort_t fTriggerVZERO = aod->GetVZEROData()->GetTriggerBits();
  UInt_t fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  
  fTriggerInputsMC[0] = fTriggerVZERO & (1 << 12); //0VBA VZERO A
  fTriggerInputsMC[1] = fTriggerVZERO & (1 << 13); //0VBC VZERO C
  fTriggerInputsMC[2] = fTriggerAD & (1 << 12);   //0UBA ADA
  fTriggerInputsMC[3] = fTriggerAD & (1 << 13);   //0UBC ADC
  fTriggerInputsMC[4] = fL0inputs & (1 << 22);  //0OMU TOF two hits with topology
  fTriggerInputsMC[5] = fL0inputs & (1 << 19);	//0OM2 TOF two hits
  					
  //SPD inputs
  const Int_t bcMod4 = TMath::Nint(hBCmod4->GetRandom());
  const AliAODTracklets *mult = aod->GetMultiplicity();
  fFOFiredChips = mult->GetFastOrFiredChips();
  Int_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=0;
  Int_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=0;

  Int_t nInner(0), nOuter(0);
  for (Int_t i(0); i<1200; ++i) {
    const Double_t eff = hSPDeff->GetBinContent(1+i, 1+bcMod4);
    Bool_t isFired = (mult->TestFastOrFiredChips(i)) && (gRandom->Uniform(0,1) < eff);
    if (i<400) {
      vPhiInner[i/20] += isFired;
      nInner += isFired;
    } else {
      vPhiOuter[(i-400)/20] += isFired;
      nOuter += isFired;
    }
  }
 
  Int_t fired = 0;
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
  if (nOuter > 0 || nInner > 0) fTriggerInputsMC[6] = kTRUE;
  //0SM2 - Two hits on outer layer
  if (nOuter > 1) fTriggerInputsMC[7] = kTRUE;
  //0STP - Topological SPD trigger (two pairs)
  if (fired != 0) fTriggerInputsMC[8] = kTRUE;
  //0SH1 - More then 6 hits on outer layer
  if (nOuter >= 7) fTriggerInputsMC[9] = kTRUE;
  //0SH1 2017 - Two hits on inner and outer layer
  if (nInner >= 2 && nOuter >= 2) fTriggerInputsMC[10] = kTRUE;
  

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
void AliAnalysisTaskUpcPsi2s::RunESDtrig()
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
void AliAnalysisTaskUpcPsi2s::RunESDhist()
{

  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  Double_t muonMass = partMuon->Mass();
  
  TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
  Double_t electronMass = partElectron->Mass();
  
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();
  
  TParticlePDG *partProton = pdgdat->GetParticle( 2212 );
  Float_t protonMass = partProton->Mass();

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;
  
    //input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  fDataFilnam->Clear();
  fDataFilnam->SetString(filnam);
  fEvtNum = ((TTree*) GetInputData(0))->GetTree()->GetReadEntry();

  fHistNeventsJPsi->Fill(1);
  fHistNeventsPsi2s->Fill(1);

  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  if(!isMC && !trigger.Contains("CCUP") ) return;
  
  fHistNeventsJPsi->Fill(2);
  fHistNeventsPsi2s->Fill(2);
  
  
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];

  //primary vertex
  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
  fVtxContrib = fESDVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsJPsi->Fill(3);
  fHistNeventsPsi2s->Fill(3);

  //VZERO, AD
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDAD *fADdata = esd->GetADData();
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fADdata){
  	fADAdecision = fADdata->GetADADecision();
  	fADCdecision = fADdata->GetADCDecision();
	}
 
  if(fV0Adecision != AliESDVZERO::kV0Empty || fV0Cdecision != AliESDVZERO::kV0Empty) return;
  if(fADAdecision != AliESDAD::kADEmpty || fADCdecision != AliESDAD::kADEmpty) return;
   
  fHistNeventsJPsi->Fill(4);
  fHistNeventsPsi2s->Fill(4);

   Int_t nGoodTracks=0;
  //Two tracks loop
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  Double_t TrackPt[5]={0,0,0,0,0};
  
  TLorentzVector vMuon[5],vElectron[5],vProton[5],vPion[5], vJPsiCandidate;

  Float_t nSigmaMuon[5], nSigmaElectron[5],nSigmaProton[5],MeanPt;
  Short_t qPion[5];
  TLorentzVector vLepton[5], vDilepton, vPsi2sCandidate;
  Short_t qLepton[5];
  UInt_t nPion = 0,nLepton = 0;


 //Two Track loop
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 70)continue;
      if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;
      Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
      if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
      trk->GetImpactParameters(dca[0],dca[1]);
      if(TMath::Abs(dca[0]) > 2.4) continue;
      if(TMath::Abs(dca[1]) > 3.2) continue;
      
      TrackIndex[nGoodTracks] = itr;
      TrackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
      if(nGoodTracks > 4) break;   
  }//Track loop

  if(nGoodTracks == 2){
  	  fHistNeventsJPsi->Fill(5);
  	  for(Int_t iTrack=0; iTrack<2; iTrack++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[iTrack]);
				   
      		Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    		Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		Float_t fPIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
		
		vElectron[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
    		vMuon[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
    		nSigmaMuon[iTrack] = fPIDTPCMuon;
    		nSigmaElectron[iTrack] = fPIDTPCElectron;
    
    		vPion[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
	
    		vProton[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), protonMass);
		nSigmaProton[iTrack] = fPIDTOFProton;
		
		qLepton[iTrack] = trk->Charge();
				
  		}		
  	  if(qLepton[0]*qLepton[1] > 0) fHistNeventsJPsi->Fill(6);
	  if(qLepton[0]*qLepton[1] < 0){
		fHistNeventsJPsi->Fill(7);
		Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
  		Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
  		Float_t nSigmaDistProton = TMath::Sqrt(TMath::Power(nSigmaProton[0],2) + TMath::Power(nSigmaProton[1],2));
		
		fDiLeptonM = -999; 
		fDiLeptonPt = -999;

 		if(nSigmaDistProton < 5){ 
  	  		fPIDsigma = nSigmaDistProton;
  	  		vJPsiCandidate = vProton[0]+vProton[1];
  	  		fChannel = 2;
			fPt = vJPsiCandidate.Pt();
  			fY = vJPsiCandidate.Rapidity();
  			fM = vJPsiCandidate.M();
  	  		fEveTree->Fill();
  	  		}
  		if(nSigmaDistMuon < nSigmaDistElectron){
  	  		fPIDsigma = nSigmaDistMuon; 
  	  		vJPsiCandidate = vMuon[0]+vMuon[1];
  	  		fChannel = 1;
  	 		fPt = vJPsiCandidate.Pt();
  			fY = vJPsiCandidate.Rapidity();
  			fM = vJPsiCandidate.M();
  	  		fEveTree->Fill();
			fHistDiMuonMass->Fill(fM);
			fHistDiLeptonMass->Fill(fM);
  	 		}
  		if(nSigmaDistMuon > nSigmaDistElectron){ 
  	  		fPIDsigma = nSigmaDistElectron;
  	  		vJPsiCandidate = vElectron[0]+vElectron[1];
  	  		fChannel = -1;
  	  		fPt = vJPsiCandidate.Pt();
  			fY = vJPsiCandidate.Rapidity();
  			fM = vJPsiCandidate.M();
  	  		fEveTree->Fill();
			fHistDiElectronMass->Fill(fM);
			fHistDiLeptonMass->Fill(fM);
  	  		}
		}

  }
  if(nGoodTracks == 4){
    	  MeanPt = GetMedian(TrackPt);
  	  fHistNeventsPsi2s->Fill(6);
  	  for(Int_t iTrack=0; iTrack<4; iTrack++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[iTrack]);

		if(trk->Pt() > MeanPt){
      			qLepton[nLepton] = trk->Charge();
			Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    			Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);

			nSigmaMuon[nLepton] = fPIDTPCMuon;
    			nSigmaElectron[nLepton] = fPIDTPCElectron;
			
			if(TMath::Abs(fPIDTPCMuon) < TMath::Abs(fPIDTPCElectron)){	
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				nLepton++;
				}
      			if(TMath::Abs(fPIDTPCMuon) > TMath::Abs(fPIDTPCElectron)){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				nLepton++;
				}
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}
	}
	if((nLepton == 2) && (nPion == 2)){
		fHistNeventsPsi2s->Fill(7);
		if(qLepton[0]*qLepton[1] > 0) fHistNeventsPsi2s->Fill(8);
		if(qPion[0]*qPion[1] > 0) fHistNeventsPsi2s->Fill(9);
		if((qLepton[0]*qLepton[1] > 0) && (qPion[0]*qPion[1] > 0)) fHistNeventsPsi2s->Fill(10);
		if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0)){
			fHistNeventsPsi2s->Fill(11);
			fHistNeventsPsi2s->Fill(12); 
  			vPsi2sCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  			vDilepton = vLepton[0]+vLepton[1];
			fHistPsi2sMassVsPt->Fill(vPsi2sCandidate.M(),vPsi2sCandidate.Pt());
			
			Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
			Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
			fChannel = 0;
			fPt = vPsi2sCandidate.Pt();
  			fY = vPsi2sCandidate.Rapidity();
  			fM = vPsi2sCandidate.M();
			fDiLeptonM = vDilepton.M();
			fDiLeptonPt = vDilepton.Pt();
			if((nSigmaDistMuon < nSigmaDistElectron)){
				fPIDsigma = nSigmaDistMuon; 
				fChannel = 1;
				fEveTree->Fill();
				}
			if((nSigmaDistMuon > nSigmaDistElectron)){
				fPIDsigma = nSigmaDistElectron;
				fChannel = -1;
				fEveTree->Fill();
				}
			}
		}
  }

  
  PostData(4, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunESDtree()
{

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
  for(Int_t i=0; i<ntrgMB; i++) {
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
  fVtxErr[1] = CovMatx[2];
  fVtxErr[2] = CovMatx[5];
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
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  //Two Track loop
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;
    
    if(fTracking == 0){
      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 70)continue;
      if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;
      Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
      if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
      trk->GetImpactParameters(dca[0],dca[1]);
      if(!isMC){
      	if(TMath::Abs(dca[1]) > 2) continue;
      	Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
        if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
	}
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!(trk->GetStatus() & AliESDtrack::kITSpureSA) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->GetTPCchi2()/trk->GetITSNcls() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
      
      if(nGoodTracks > 2) break;   
  }//Track loop

  fJPsiESDTracks->Clear("C");
  if(nGoodTracks == 2){
  
  	  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
	  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  	  Double_t muonMass = partMuon->Mass();  
          TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
          Double_t electronMass = partElectron->Mass();  
  	  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  	  Double_t pionMass = partPion->Mass();
  
  	  Double_t KFcov[21];
  	  Double_t KFpar[6];
	  Double_t KFmass = pionMass;
	  Double_t fRecTPCsignal;
  	  AliKFParticle *KFpart[2];
  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(esd->GetMagneticField()); 
	  
  	  for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
		if(fESDVertex->UsesTrack(TrackIndex[i]))fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);
				
		new((*fJPsiESDTracks)[i]) AliESDtrack(*trk); 
		
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
  
  if(!isMC) fJPsiTree ->Fill();
  }
  
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
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
    if(fTracking == 1){
      if(!(trk->GetStatus() & AliESDtrack::kITSpureSA) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->GetTPCchi2()/trk->GetITSNcls() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      }
      
      if(nGoodTracks > 4) break;   
  }//Track loop
  
  fPsi2sESDTracks->Clear("C");
  if(nGoodTracks == 4){

  	  AliKFVertex *KFvtx = new AliKFVertex();
  	  KFvtx->SetField(esd->GetMagneticField()); 

  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);

		if(fESDVertex->UsesTrack(TrackIndex[i]))fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fPsi2sESDTracks)[i]) AliESDtrack(*trk);
		
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
  		}
		
  if(!isMC) fPsi2sTree ->Fill();
  }
  
  if(isMC){
  	fJPsiTree ->Fill();
	fPsi2sTree ->Fill();
  }
 
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);

}//RunESD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunESDMC(AliESDEvent* esd)
{
  for(Int_t i=0; i<ntrgMB; i++) fTriggerInputsMC[i] = kFALSE;
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
void AliAnalysisTaskUpcPsi2s::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcPsi2s::GetMedian(Double_t *daArray) {
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

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunAODsystematics(AliAODEvent* aod)
{

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
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vCandidate, vDilepton;
  Short_t qLepton[4];
  UInt_t nLepton=0, nHighPt=0;
  Double_t fRecTPCsignal[5], fRecTPCsignalDist;
  Int_t fChannel = 0;

  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  Double_t muonMass = partMuon->Mass();
  
  TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
  Double_t electronMass = partElectron->Mass();
  

  
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
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[k]));
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
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[k]));
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

}
