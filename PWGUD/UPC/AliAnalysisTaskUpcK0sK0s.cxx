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
#include <bitset>

// root headers
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TString.h"
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

// my headers
#include "AliAnalysisTaskUpcK0sK0s.h"

ClassImp(AliAnalysisTaskUpcK0sK0s);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcK0sK0s::AliAnalysisTaskUpcK0sK0s() 
  : AliAnalysisTaskSE(),fType(0),fPIDResponse(0),fK0sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fZDCAtime(0),fZDCCtime(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fK0sAODv0s(0),fK0sAODTracks(0),fK0sESDTracks(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcK0sK0s


//_____________________________________________________________________________
AliAnalysisTaskUpcK0sK0s::AliAnalysisTaskUpcK0sK0s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fPIDResponse(0),fK0sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),fSpdVtxContrib(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZNAenergy(0),fZNCenergy(0), fZPAenergy(0),fZPCenergy(0),fZDCAtime(0),fZDCCtime(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fK0sAODv0s(0),fK0sAODTracks(0),fK0sESDTracks(0),
    fListTrig(0),fHistCcup4TriggersPerRun(0), fHistCcup7TriggersPerRun(0), fHistCcup2TriggersPerRun(0),fHistCint1TriggersPerRun(0),fHistCint6TriggersPerRun(0), fHistC0tvxAndCint1TriggersPerRun(0)

{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());

}//AliAnalysisTaskUpcK0sK0s

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) fTrigger[i] = kFALSE;
  
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
  for(Int_t i=0; i<3; i++){
  	fVtxPos[i] = -666; 
	fVtxErr[i] = -666;
	fSpdVtxPos[i] = -666;
	}


}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcK0sK0s::~AliAnalysisTaskUpcK0sK0s() 
{
  // Destructor
  if(fK0sTree){
     delete fK0sTree;
     fK0sTree = 0x0;
  }
  if(fListTrig){
     delete fListTrig;
     fListTrig = 0x0;
  }

}//~AliAnalysisTaskUpcK0sK0s


//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::UserCreateOutputObjects()
{

  //PID response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

    //vertices
  fK0sAODv0s = new TClonesArray("AliAODv0", 1000);
  
    //tracks
  fK0sAODTracks = new TClonesArray("AliAODTrack", 1000);
  fK0sESDTracks = new TClonesArray("AliESDtrack", 1000);

  //output tree with K0s candidate events
  fK0sTree = new TTree("fK0sTree", "fK0sTree");
  fK0sTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fK0sTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fK0sTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fK0sTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fK0sTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fK0sTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fK0sTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fK0sTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fK0sTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fK0sTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fK0sTree ->Branch("fSpdVtxContrib", &fSpdVtxContrib, "fSpdVtxContrib/I");
  
  fK0sTree ->Branch("fPIDTPCMuon", &fPIDTPCMuon[0], "fPIDTPCMuon[4]/D");
  fK0sTree ->Branch("fPIDTPCElectron", &fPIDTPCElectron[0], "fPIDTPCElectron[4]/D");
  fK0sTree ->Branch("fPIDTPCPion", &fPIDTPCPion[0], "fPIDTPCPion[4]/D");
  fK0sTree ->Branch("fPIDTPCKaon", &fPIDTPCKaon[0], "fPIDTPCKaon[4]/D");
  fK0sTree ->Branch("fPIDTPCProton", &fPIDTPCProton[0], "fPIDTPCProton[4]/D");
  
  fK0sTree ->Branch("fPIDTOFMuon", &fPIDTOFMuon[0], "fPIDTOFMuon[4]/D");
  fK0sTree ->Branch("fPIDTOFElectron", &fPIDTOFElectron[0], "fPIDTOFElectron[4]/D");
  fK0sTree ->Branch("fPIDTOFPion", &fPIDTOFPion[0], "fPIDTOFPion[4]/D");
  fK0sTree ->Branch("fPIDTOFKaon", &fPIDTOFKaon[0], "fPIDTOFKaon[4]/D");
  fK0sTree ->Branch("fPIDTOFProton", &fPIDTOFProton[0], "fPIDTOFProton[4]/D");
  
  fK0sTree ->Branch("fIsVtxContributor", &fIsVtxContributor[0], "fIsVtxContributor[4]/O");
  
  fK0sTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fK0sTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fK0sTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fK0sTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fK0sTree ->Branch("fSpdVtxPos", &fSpdVtxPos[0], "fSpdVtxPos[3]/D");
  
  fK0sTree ->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/D");
  fK0sTree ->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/D");
  fK0sTree ->Branch("fZPAenergy", &fZPAenergy, "fZPAenergy/D");
  fK0sTree ->Branch("fZPCenergy", &fZPCenergy, "fZPCenergy/D");
  fK0sTree ->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/D");
  fK0sTree ->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/D");
  fK0sTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fK0sTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fK0sTree ->Branch("fDataFilnam", &fDataFilnam);
  fK0sTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fK0sTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fK0sTree ->Branch("fK0sESDTracks", &fK0sESDTracks);
  }
  if( fType == 1 ) {
    fK0sTree ->Branch("fK0sAODv0s", &fK0sAODv0s);
    fK0sTree ->Branch("fK0sAODTracks", &fK0sAODTracks);
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
  
  fHistCint6TriggersPerRun = new TH1D("fHistCint6TriggersPerRun", "fHistCint6TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistCint6TriggersPerRun);
  
  fHistC0tvxAndCint1TriggersPerRun = new TH1D("fHistC0tvxAndCint1TriggersPerRun", "fHistC0tvxAndCint1TriggersPerRun", 33000, 167000.5, 200000.5);
  fListTrig->Add(fHistC0tvxAndCint1TriggersPerRun);
  
  PostData(1, fK0sTree);
  PostData(2, fListTrig);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ){ 
  	RunESDtrig();
  	RunESDtree();
	}
  if( fType == 1 ) RunAOD();
	

}//UserExec

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunAOD()
{
  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  //input data
  AliAODHeader * header = dynamic_cast<AliAODHeader*>(aod->GetHeader());
if(!header) AliFatal("Not a standard AOD");

  //fDataFilnam = header->GetESDFileName();
  fEvtNum = header->GetEventNumberESDFile();
  fRunNum = aod->GetRunNumber();
  

  //Trigger
  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  fTrigger[0]  = trigger.Contains("CCUP4-B"); // Central UPC Pb-Pb 2011
  fTrigger[1]  = trigger.Contains("CCUP2-B"); // Double gap
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // Central UPC p-Pb 2013
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if(!isTriggered ) return;

  //trigger inputs
  fL0inputs = ((AliVAODHeader*)aod->GetHeader())->GetL0TriggerInputs();
  fL1inputs = ((AliVAODHeader*)aod->GetHeader())->GetL1TriggerInputs();

  //Event identification
  fPerNum = aod ->GetPeriodNumber();
  fOrbNum = aod ->GetOrbitNumber();
  fBCrossNum = aod ->GetBunchCrossNumber();

  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();

  //Tracklets
  fNtracklets = aod->GetTracklets()->GetNumberOfTracklets();


  Int_t nGoodV0s=0;
  Int_t V0Index[3] = {-1,-1,-1};
  Int_t V0TrackID[6] = {-1,-1,-1,-1,-1,-1};

  Int_t nGoodTracks=0;
  Int_t TrackID[5] = {-1,-1,-1,-1,-1};

  //V0s loop
  for(Int_t iV0=0; iV0<aod ->GetNumberOfV0s(); iV0++) {
    AliAODv0 *v0 = aod->GetV0(iV0);
    if( !v0 ) continue;
    Bool_t lOnFlyStatus = v0->GetOnFlyStatus();
    if (lOnFlyStatus) continue;
    
    AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
    AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
    if (!pTrack || !nTrack) continue;

    if ( pTrack->Charge() == nTrack->Charge())continue;

      if(!(pTrack->TestFilterBit(1<<0))) continue;
      if(!(nTrack->TestFilterBit(1<<0))) continue;
      if(!(pTrack->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(nTrack->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(pTrack->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(!(nTrack->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(pTrack->GetTPCNcls() < 50)continue;
      if(nTrack->GetTPCNcls() < 50)continue;
      if(pTrack->Chi2perNDF() > 4)continue;
      if(nTrack->Chi2perNDF() > 4)continue;
      
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      if(!pTrack->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;
      if(!nTrack->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;
      
      V0Index[nGoodV0s] = iV0;
      V0TrackID[2*nGoodV0s] = pTrack->GetID();
      V0TrackID[2*nGoodV0s+1] = nTrack->GetID();
      nGoodV0s++; 
      if(nGoodV0s > 2) break;
  }//V0s loop
  
  //Track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if(!trk) AliFatal("Not a standard AOD");
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      if(!trk->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;
     
      TrackID[nGoodTracks] = trk->GetID();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
    
  if(nGoodV0s == 2 && nGoodTracks == 4){
  	//SortArray(TrackID);
  	//SortArray(V0TrackID);
	//for{Int_t i=0; i<4; i++} if (TrackID[i] != V0TrackID[i]) return;
  	for(Int_t i=0; i<2; i++){
	  	AliAODv0 *v0 = aod->GetV0(V0Index[i]);
		AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
    		AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
				
		new((*fK0sAODv0s)[i]) AliAODv0(*v0);
		new((*fK0sAODTracks)[2*i]) AliAODTrack(*pTrack);
		new((*fK0sAODTracks)[2*i+1]) AliAODTrack(*nTrack);
		}
  fK0sTree ->Fill();
  
  }
  
  PostData(1, fK0sTree);   

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunESDtrig()
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
  if(trigger.Contains("CINT6-B")) fHistCint6TriggersPerRun->Fill(fRunNum); //CINT6 triggers
  
  fL0inputs = esd->GetHeader()->GetL0TriggerInputs();
  if(trigger.Contains("CINT1-B") && (fL0inputs & (1 << 3))) fHistC0tvxAndCint1TriggersPerRun->Fill(fRunNum); //0TVX triggers in CINT1 events
 
PostData(2, fListTrig);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunESDtree()
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
  if(!isTriggered ) return;
  
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
  fSpdVtxContrib = fSPDVertex->GetNContributors();
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

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 20)continue;
      fNLooseTracks++; 
  }//Track loop -loose cuts
  
  Int_t nGoodNegTracks=0;
  Int_t nGoodPosTracks=0;
  Int_t NegTrackIndex[3] = {-1,-1,-1};
  Int_t PosTrackIndex[3] = {-1,-1,-1};
  
  //Four track loop
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
     
      if(trk->Charge()<0){ 
      		NegTrackIndex[nGoodNegTracks] = itr;
      		nGoodNegTracks++;
		}
      if(trk->Charge()>0){ 
      		PosTrackIndex[nGoodPosTracks] = itr;
      		nGoodPosTracks++;
		}
      if(nGoodNegTracks > 2 || nGoodPosTracks > 2) break;   
  }//Track loop
  
  Int_t nGoodV0s=0;
  if(nGoodNegTracks == 2 && nGoodPosTracks == 2){
  
  	Double_t fChi2max=33.; //max chi2
  	// Double_t fDmin=0.01;  //min imp parameter for the daughters
  	Double_t fDCAmax=2.5;  //max DCA between the daughter tracks
  	Double_t fCPAmin=0.80;  //min cosine of V0's pointing angle
  	Double_t fRmin=0.1;    //min radius of the fiducial volume
  	Double_t fRmax=200.;   //max radius of the fiducial volume
	
	for(Int_t i=0; i<2; i++){
  		for(Int_t j=0; j<2; j++){
			AliESDtrack *trkN = (AliESDtrack*) esd->GetTrack(NegTrackIndex[i]);
			AliESDtrack *trkP = (AliESDtrack*) esd->GetTrack(PosTrackIndex[j]);
			Int_t nidx =trkN->GetID();
			Int_t pidx =trkP->GetID();
			Double_t b = esd->GetMagneticField();
			//AliV0vertexer
			// Double_t dxn=TMath::Abs(trkN->GetD(fVtxPos[0],fVtxPos[1],b)), dxp=TMath::Abs(trkP->GetD(fVtxPos[0],fVtxPos[1],b));
			Double_t xn, xp, dca=trkN->GetDCA(trkP,b,xn,xp);
			
        		AliExternalTrackParam nt(*trkN), pt(*trkP);
        		nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);

        		AliESDv0 vertex(nt,nidx,pt,pidx);
        		if (vertex.GetChi2V0() > fChi2max) continue;
	 
        		Double_t x=vertex.Xv(), y=vertex.Yv();
        		Double_t r2=x*x + y*y;

			Float_t cpa=vertex.GetV0CosineOfPointingAngle(fVtxPos[0],fVtxPos[1],fVtxPos[2]);
			//std::cout<<"DCA between = "<<dca<<" Radius = "<<(xn+xp)<<std::endl;
			//std::cout<<"CPA = "<<cpa<<" Radius 2D = "<<r2<<std::endl;
			
			if (cpa < fCPAmin) continue;
			if (dca > fDCAmax) continue;	
			//if(dxn<fDmin || dxp<fDmin) continue;
        		if ((xn+xp) > 2*fRmax) continue;
        		if ((xn+xp) < 2*fRmin) continue;
			if (r2 < fRmin*fRmin) continue;
        		if (r2 > fRmax*fRmax) continue;

			nGoodV0s++;
			}
		}
  }
   
  fK0sESDTracks->Clear("C");
  if(nGoodNegTracks == 2 && nGoodPosTracks == 2 && nGoodV0s>1){
  
  	  for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(NegTrackIndex[i]);
		
		if(fESDVertex->UsesTrack(NegTrackIndex[i]))fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fK0sESDTracks)[i]) AliESDtrack(*trk);
		
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
	for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(PosTrackIndex[i]);
		
		if(fESDVertex->UsesTrack(PosTrackIndex[i]))fIsVtxContributor[i] = kTRUE;
		else fIsVtxContributor[i] = kFALSE;
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fK0sESDTracks)[i+2]) AliESDtrack(*trk);	
		
		fPIDTPCMuon[i+2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		fPIDTPCElectron[i+2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		fPIDTPCPion[i+2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		fPIDTPCKaon[i+2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		fPIDTPCProton[i+2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
		
		fPIDTOFMuon[i+2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
		fPIDTOFElectron[i+2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
		fPIDTOFPion[i+2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
		fPIDTOFKaon[i+2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
		fPIDTOFProton[i+2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
		
  		}
  fK0sTree->Fill();
  }
  
  PostData(1, fK0sTree);

}//RunESDtree

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate

