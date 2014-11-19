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
  : AliAnalysisTaskSE(),fType(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fPhiTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFtrig1(0), fTOFtrig2(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
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
  : AliAnalysisTaskSE(name),fType(0),isMC(kFALSE),fRunTree(kTRUE),fRunHist(kTRUE),fRunSystematics(kFALSE),fPIDResponse(0),fPhiTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFtrig1(0), fTOFtrig2(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
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
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());

}//AliAnalysisTaskUpcPhi

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) fTrigger[i] = kFALSE;
  for(Int_t i=0; i<4; i++) {
  	fTOFphi[i] = -666;
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
  if(fPhiTree){
     delete fPhiTree;
     fPhiTree = 0x0;
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

 //output tree with Phi candidate events
  fPhiTree = new TTree("fPhiTree", "fPhiTree");
  fPhiTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fPhiTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fPhiTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fPhiTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fPhiTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fPhiTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fPhiTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fPhiTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fPhiTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fPhiTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  
  fPhiTree ->Branch("fTOFtrig1", &fTOFtrig1, "fTOFtrig1/O");
  fPhiTree ->Branch("fTOFtrig2", &fTOFtrig2, "fTOFtrig2/O");
  fPhiTree ->Branch("fTOFphi", &fTOFphi[0], "fTOFphi[2]/D");
  
  fPhiTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[2]/D");
  fPhiTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[2]/D");
  fPhiTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[2]/D");
  fPhiTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[2]/D");
  fPhiTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[2]/D");
  
  fPhiTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[2]/D");
  fPhiTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[2]/D");
  fPhiTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[2]/D");
  fPhiTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[2]/D");
  fPhiTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[2]/D");
  
  fPhiTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fPhiTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fPhiTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fPhiTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fPhiTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  fPhiTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fPhiTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fPhiTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fPhiTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fPhiTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fPhiTree ->Branch("fDataFilnam", &fDataFilnam);
  fPhiTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fPhiTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fPhiTree ->Branch("fPhiESDTracks", &fPhiESDTracks);
  }
  if( fType == 1 ) {
    fPhiTree ->Branch("fPhiAODTracks", &fPhiAODTracks);
  }
  if(isMC) {
    fPhiTree ->Branch("fGenPart", &fGenPart);
    fPhiTree ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[2]/O");
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
  
  PostData(1, fPhiTree);
  PostData(2, fListTrig);
  PostData(3, fListHist);

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
  
  if(trigger.Contains("CINT1")) fHistCint1TriggersPerRun->Fill(fRunNum); //CINT1 triggers
  
  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  if(trigger.Contains("CINT1") && (fL0inputs & (1 << 3))) fHistC0tvxAndCint1TriggersPerRun->Fill(fRunNum); //0TVX triggers in CINT1 events
  
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
  
  fNLooseTracks = 0;
  
  //Track loop - loose cuts
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<1))) continue;

      fNLooseTracks++; 
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
   //Two track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<1))) continue;

      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetITSNcls() < 3)continue;
      if(trk->Chi2perNDF() > 4)continue;
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
		
		new((*fPhiAODTracks)[i]) AliAODTrack(*trk);
		}
			
  if(!isMC) fPhiTree ->Fill();
  }
  if(isMC) fPhiTree ->Fill();
  
  
  PostData(1, fPhiTree);

}//RunAOD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunAODMC(AliAODEvent *aod)
{

  fL0inputs = aod->GetHeader()->GetL0TriggerInputs();
  fTriggerInputsMC[0] = kFALSE;//0SM2
  fTriggerInputsMC[1] = fL0inputs & (1 << 0);//0VBA
  fTriggerInputsMC[2] = fL0inputs & (1 << 1);//0VBC
  fTriggerInputsMC[3] = fL0inputs & (1 << 9);//0OMU

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

}//RunAODMC


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDtrig()
{


}
//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDhist()
{


}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDtree()
{


}//RunESD


//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::RunESDMC(AliESDEvent* esd)
{


}//RunESDMC



//_____________________________________________________________________________
void AliAnalysisTaskUpcPhi::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate
