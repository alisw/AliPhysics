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
  : AliAnalysisTaskSE(),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),fK0sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fK0sAODv0s(0),fK0sESDv0s(0),fK0sAODTracks(0),fK0sESDTracks(0),
    fListHist(0),fHistTriggersPerRun(0),fHistK0sMass(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcK0sK0s


//_____________________________________________________________________________
AliAnalysisTaskUpcK0sK0s::AliAnalysisTaskUpcK0sK0s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),fK0sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fK0sAODv0s(0),fK0sESDv0s(0),fK0sAODTracks(0),fK0sESDTracks(0),
    fListHist(0),fHistTriggersPerRun(0),fHistK0sMass(0)

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

}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcK0sK0s::~AliAnalysisTaskUpcK0sK0s() 
{
  // Destructor
  if(fK0sTree){
     delete fK0sTree;
     fK0sTree = 0x0;
  }
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
  }

}//~AliAnalysisTaskUpcK0sK0s


//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::UserCreateOutputObjects()
{

    //vertices
  fK0sAODv0s = new TClonesArray("AliAODv0", 1000);
  fK0sESDv0s = new TClonesArray("AliESDv0", 1000);
  
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
  fK0sTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fK0sTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fK0sTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fK0sTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fK0sTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fK0sTree ->Branch("fDataFilnam", &fDataFilnam);
  fK0sTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fK0sTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L"); 		       
  if( fType == 0 ) {
    fK0sTree ->Branch("fK0sESDv0s", &fK0sESDv0s);
    fK0sTree ->Branch("fK0sESDTracks", &fK0sESDTracks);
  }
  if( fType == 1 ) {
    fK0sTree ->Branch("fK0sAODv0s", &fK0sAODv0s);
    fK0sTree ->Branch("fK0sAODTracks", &fK0sAODTracks);
  }
 
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
  fHistTriggersPerRun = new TH1D("fHistTriggersPerRun", "fHistTriggersPerRun", 3000, 167000.5, 170000.5);
  fListHist->Add(fHistTriggersPerRun);
    
  fHistK0sMass = new TH2D("fHistK0sMass","fHistK0sMass",20,0.45,0.55,20,0.45,0.55);
  fListHist->Add(fHistK0sMass);
  
  PostData(1, fK0sTree);
  PostData(2, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ) RunESD();
  if( fType == 1 ){ 
  	if(fRunHist) RunAODhist();
	if(fRunTree) RunAODtree();
	}

}//UserExec
//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunAODhist()
{
  
  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;
  //cout<<"File: "<<(aod->GetHeader()->GetESDFileName()).Data()<<" Event: "<<aod->GetHeader()->GetEventNumberESDFile()<<endl;

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if( !trigger.Contains("CCUP4-B") ) return;
  
  fRunNum = aod ->GetRunNumber();
  fHistTriggersPerRun->Fill(fRunNum);


  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();
  if(fVtxContrib < 2) return;


  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;

  AliAODZDC *fZDCdata = aod->GetZDCData();
  
  fZDCAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  if( fZDCAenergy > 8200 || fZDCCenergy > 8200) return;
  
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

      if(!(pTrack->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(!(nTrack->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(!(pTrack->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(nTrack->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
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
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

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
	//for(Int_t i=0; i<4; i++) if (TrackID[i] != V0TrackID[i]) return;
	AliAODv0 *v00 = aod->GetV0(V0Index[0]);	
	AliAODv0 *v01 = aod->GetV0(V0Index[1]);				
	fHistK0sMass->Fill(v00->MassK0Short(),v01->MassK0Short());
			
  }
  
  PostData(2, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunAODtree()
{
  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  //input data
  fDataFilnam = aod->GetHeader()->GetESDFileName();
  fEvtNum = aod->GetHeader()->GetEventNumberESDFile();
  fRunNum = aod->GetRunNumber();

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if( !trigger.Contains("CCUP4-B") ) return;

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

  //Tracklets
  fNtracklets = aod->GetTracklets()->GetNumberOfTracklets();

  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  AliAODZDC *fZDCdata = aod->GetZDCData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  fZDCAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  
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
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

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
void AliAnalysisTaskUpcK0sK0s::SortArray(Double_t *dArray) {
    for (Int_t i = 3; i > 0; --i) {
        for (Int_t j = 0; j < i; ++j) {
            if (dArray[j] > dArray[j+1]) {
                Double_t dTemp = dArray[j];
                dArray[j] = dArray[j+1];
                dArray[j+1] = dTemp;
            }
        }
    }
}

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::RunESD()
{

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcK0sK0s::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate






























