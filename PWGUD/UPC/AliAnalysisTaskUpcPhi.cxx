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
#include "AliTriggerAnalysis.h"

// my headers
#include "AliAnalysisTaskUpcPhi.h"

ClassImp(AliAnalysisTaskUpcPhi);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcPhi::AliAnalysisTaskUpcPhi() 
  : AliAnalysisTaskSE(),fType(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fITSTree(0),fTPCTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseITSTracks(0),fNLooseTPCTracks(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fPhiAODTracks(0),fPhiESDTracks(0),fGenPart(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fListHist(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcPhi


//_____________________________________________________________________________
AliAnalysisTaskUpcPhi::AliAnalysisTaskUpcPhi(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fITSTree(0),fTPCTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseITSTracks(0),fNLooseTPCTracks(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fPhiAODTracks(0),fPhiESDTracks(0),fGenPart(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistC0tvxAndCint1TriggersPerRun(0),
    fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0), fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fListHist(0)

{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

}//AliAnalysisTaskUpcPhi

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) fTrigger[i] = kFALSE;
  for(Int_t i=0; i<4; i++) {
	fPIDITSMuon[i] = -666;
	fPIDITSElectron[i] = -666;
	fPIDITSPion[i] = -666;
	fPIDITSKaon[i] = -666;
	fPIDITSProton[i] = -666;
	
	fPIDTPCMuon[i] = -666;
	fPIDTPCElectron[i] = -666;
	fPIDTPCPion[i] = -666;
	fPIDTPCKaon[i] = -666;
	fPIDTPCProton[i] = -666;
	
	fTriggerInputsMC[i] = kFALSE;
	}
  for(Int_t i=0; i<3; i++){
  	fVtxPos[i] = -666; 
	fVtxErr[i] = -666;
	fKfVtxPos[i] = -666;
	fSpdVtxPos[i] = -666;
	}

}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcPhi::~AliAnalysisTaskUpcPhi() 
{
  // Destructor
  if(fITSTree){
     delete fITSTree;
     fITSTree = 0x0;
  }
  if(fTPCTree){
     delete fTPCTree;
     fTPCTree = 0x0;
  }
  if(fListTrig){
     delete fListTrig;
     fListTrig = 0x0;
  }
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
  }

}//~AliAnalysisTaskUpcPhi


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::UserCreateOutputObjects()
{
  //PID response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  //input file
  fDataFilnam = new TObjString();
  fDataFilnam->SetString("");

    //tracks
  fPhiAODTracks = new TClonesArray("AliAODTrack", 1000);
  fPhiESDTracks = new TClonesArray("AliESDtrack", 1000);
  fGenPart = new TClonesArray("TParticle", 1000);

 //output tree with Phi ITSsa candidate events
  fITSTree = new TTree("fITSTree", "fITSTree");
  fITSTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fITSTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fITSTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fITSTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fITSTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fITSTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fITSTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fITSTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fITSTree ->Branch("fNLooseITSTracks", &fNLooseITSTracks, "fNLooseITSTracks/s");
  fITSTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  
  fITSTree ->Branch("fPIDITSMuon", &fPIDITSMuon[0], "fPIDITSMuon[2]/D");
  fITSTree ->Branch("fPIDITSElectron", &fPIDITSElectron[0], "fPIDITSElectron[2]/D");
  fITSTree ->Branch("fPIDITSPion", &fPIDITSPion[0], "fPIDITSPion[2]/D");
  fITSTree ->Branch("fPIDITSKaon", &fPIDITSKaon[0], "fPIDITSKaon[2]/D");
  fITSTree ->Branch("fPIDITSProton", &fPIDITSProton[0], "fPIDITSProton[2]/D");
  
  fITSTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fITSTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fITSTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fITSTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fITSTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fITSTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fITSTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fITSTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fITSTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fITSTree ->Branch("fDataFilnam", &fDataFilnam);
  fITSTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fITSTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fITSTree ->Branch("fPhiESDTracks", &fPhiESDTracks);
  }
  if( fType == 1 ) {
    fITSTree ->Branch("fPhiAODTracks", &fPhiAODTracks);
  }
  if(isMC) {
    fITSTree ->Branch("fGenPart", &fGenPart);
    fITSTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[2]/O");
  }

 //output tree with Phi ITS-TPC candidate events
  fTPCTree = new TTree("fTPCTree", "fTPCTree");
  fTPCTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fTPCTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fTPCTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fTPCTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fTPCTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fTPCTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fTPCTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fTPCTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fTPCTree ->Branch("fNLooseTPCTracks", &fNLooseTPCTracks, "fNLooseTPCTracks/s");
  fTPCTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  
  fTPCTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[2]/D");
  fTPCTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[2]/D");
  fTPCTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[2]/D");
  fTPCTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[2]/D");
  fTPCTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[2]/D");
  
  fTPCTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fTPCTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fTPCTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fTPCTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  fTPCTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fTPCTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fTPCTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fTPCTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fTPCTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fTPCTree ->Branch("fDataFilnam", &fDataFilnam);
  fTPCTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fTPCTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fTPCTree ->Branch("fPhiESDTracks", &fPhiESDTracks);
  }
  if( fType == 1 ) {
    fTPCTree ->Branch("fPhiAODTracks", &fPhiAODTracks);
  }
  if(isMC) {
    fTPCTree ->Branch("fGenPart", &fGenPart);
    fTPCTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[2]/O");
  }

  
  fListTrig = new TList();
  fListTrig ->SetOwner();
  
  fHistCcup4TriggersPerRun = new TH1D("fHistCcup4TriggersPerRun", "fHistCcup4TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCcup4TriggersPerRun);
  
  fHistCcup7TriggersPerRun = new TH1D("fHistCcup7TriggersPerRun", "fHistCcup7TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCcup7TriggersPerRun);
  
  fHistCcup2TriggersPerRun = new TH1D("fHistCcup2TriggersPerRun", "fHistCcup2TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCcup2TriggersPerRun);
  
  fHistCint1TriggersPerRun = new TH1D("fHistCint1TriggersPerRun", "fHistCint1TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCint1TriggersPerRun);
  
  fHistC0tvxAndCint1TriggersPerRun = new TH1D("fHistC0tvxAndCint1TriggersPerRun", "fHistC0tvxAndCint1TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistC0tvxAndCint1TriggersPerRun);
  
  fHistZedTriggersPerRun = new TH1D("fHistZedTriggersPerRun", "fHistZedTriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistZedTriggersPerRun);

  fHistCvlnTriggersPerRun = new TH1D("fHistCvlnTriggersPerRun", "fHistCvlnTriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCvlnTriggersPerRun);
  
  fHistMBTriggersPerRun = new TH1D("fHistMBTriggersPerRun", "fHistMBTriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistMBTriggersPerRun);
  
  fHistCentralTriggersPerRun = new TH1D("fHistCentralTriggersPerRun", "fHistCentralTriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCentralTriggersPerRun);
  
  fHistSemiCentralTriggersPerRun = new TH1D("fHistSemiCentralTriggersPerRun", "fHistSemiCentralTriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistSemiCentralTriggersPerRun);
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
  PostData(1, fITSTree);
  PostData(2, fTPCTree);
  PostData(3, fListTrig);
  PostData(4, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::UserExec(Option_t *) 
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
void AliAnalysisTaskUpcPhi::RunAODtrig()
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
void AliAnalysisTaskUpcPhi::RunAODhist()
{


}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunAODtree()
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
  
  fTrigger[0]  = trigger.Contains("CCUP4-B"); // Central UPC Pb-Pb 2011
  fTrigger[1]  = trigger.Contains("CCUP2-B"); // Double gap
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // Central UPC p-Pb 2013
  fTrigger[3]  = trigger.Contains("CINT1");
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if(!isMC && !isTriggered ) return;

  //trigger inputs
  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  fL1inputs = aod->GetHeader()->GetL1TriggerInputs();  

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
  if(fSPDVertex){
  	fSpdVtxPos[0] = fSPDVertex->GetX();
	fSpdVtxPos[1] = fSPDVertex->GetY();
	fSpdVtxPos[2] = fSPDVertex->GetZ();
	}
  else{
  	fSpdVtxPos[0] = -666;
	fSpdVtxPos[1] = -666;
	fSpdVtxPos[2] = -666;
  	}

  //Tracklets
  fNtracklets = aod->GetTracklets()->GetNumberOfTracklets();

  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  AliAODZDC *fZDCdata = aod->GetZDCData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  fZDCAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  
  fNLooseITSTracks = 0;
  fNLooseTPCTracks = 0;
  
  //Track loop - loose cuts
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if((trk->TestFilterBit(1<<0))) fNLooseTPCTracks++;
    if((trk->TestFilterBit(1<<1))) fNLooseITSTracks++;
      
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
   //ITSsa track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<1))) continue;

      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->Chi2perNDF() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
 
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
      
  fPhiAODTracks->Clear("C");  
  if(nGoodTracks == 2){

  	  for(Int_t i=0; i<2; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));
                if(!trk) AliFatal("Not a standard AOD");
		
		fPIDITSMuon[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kMuon);
		fPIDITSElectron[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kElectron);
		fPIDITSPion[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kPion);
		fPIDITSKaon[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kKaon);
		fPIDITSProton[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kProton);
		
		new((*fPhiAODTracks)[i]) AliAODTrack(*trk);
		}
			
  fITSTree ->Fill();
  }
  
  nGoodTracks=0;
  
   //TPC track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;
    
      if(!(trk->GetStatus() & AliAODTrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
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
      
  fPhiAODTracks->Clear("C");  
  if(nGoodTracks == 2){

  	  for(Int_t i=0; i<2; i++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(TrackIndex[i]));
                if(!trk) AliFatal("Not a standard AOD");
		
		fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		new((*fPhiAODTracks)[i]) AliAODTrack(*trk);
		}
			
  fTPCTree ->Fill();
  }
  
  
  
  PostData(1, fITSTree);
  PostData(2, fTPCTree);

}//RunAOD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunAODMC(AliAODEvent *aod)
{

}//RunAODMC


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDtrig()
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
void AliAnalysisTaskUpcPhi::RunESDhist()
{


}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDtree()
{

//input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;
  
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
  fTrigger[4]  = trigger.Contains("CINT6-B"); // MB trigger
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if(!isMC && !isTriggered ) return;

  //trigger inputs
  fL0inputs = esd->GetHeader()->GetL0TriggerInputs();
  fL1inputs = esd->GetHeader()->GetL1TriggerInputs();
  
  //Event identification
  fPerNum = esd->GetPeriodNumber();
  fOrbNum = esd->GetOrbitNumber();
  fBCrossNum = esd->GetBunchCrossNumber();

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
  fSpdVtxPos[0] = fSPDVertex->GetX();
  fSpdVtxPos[1] = fSPDVertex->GetY();
  fSpdVtxPos[2] = fSPDVertex->GetZ();

  //Tracklets
  fNtracklets = esd->GetMultiplicity()->GetNumberOfTracklets();

  //VZERO, ZDC
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  fZDCAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZNCTowerEnergy()[0];

  
  fNLooseITSTracks = 0;
  fNLooseTPCTracks = 0;
  
  /*/Track loop - loose cuts
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if((trk->TestFilterBit(1<<0))) fNLooseTPCTracks++;
    if((trk->TestFilterBit(1<<1))) fNLooseITSTracks++;
      
  }//Track loop -loose cuts
  /*/
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
   //ITSsa track loop
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 4)continue;
      if(trk->GetTPCchi2()/trk->GetITSNcls() > 2.5)continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))continue;
 
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
      
  fPhiESDTracks->Clear("C");  
  if(nGoodTracks == 2){

  	  for(Int_t i=0; i<2; i++){
                AliESDtrack *trk = dynamic_cast<AliESDtrack*>(esd->GetTrack(TrackIndex[i]));
                if( !trk ) continue;
		
		fPIDITSMuon[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kMuon);
		fPIDITSElectron[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kElectron);
		fPIDITSPion[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kPion);
		fPIDITSKaon[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kKaon);
		fPIDITSProton[i] = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kProton);
		
		new((*fPhiESDTracks)[i]) AliESDtrack(*trk);
		}
			
  fITSTree ->Fill();
  }
  
  nGoodTracks=0;
  
   //TPC track loop
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
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;      
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
      
  fPhiESDTracks->Clear("C");  
  if(nGoodTracks == 2){

  	  for(Int_t i=0; i<2; i++){
                AliESDtrack *trk = dynamic_cast<AliESDtrack*>(esd->GetTrack(TrackIndex[i]));
                if( !trk ) continue;
		
		fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		new((*fPhiESDTracks)[i]) AliESDtrack(*trk);
		}
			
  fTPCTree ->Fill();
  }
  
  PostData(1, fITSTree);
  PostData(2, fTPCTree);


}//RunESD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate
