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
#include "AliAnalysisTaskUpcPsi2s.h"

ClassImp(AliAnalysisTaskUpcPsi2s);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s() 
  : AliAnalysisTaskSE(),fType(0),hCounter(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0)
{

//Dummy constructor

}//AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),hCounter(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0)
{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TH1I::Class());

}//AliAnalysisTaskUpcPsi2s

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) fTrigger[i] = kFALSE;

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
  if(hCounter){
     delete hCounter;
     hCounter = 0x0;
  }


}//~AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserCreateOutputObjects()
{
   hCounter = new TH1I("hCounter", "hCounter", 34000, 1., 34001.);

  //input file
  fDataFilnam = new TObjString();
  fDataFilnam->SetString("");

    //tracks
  fJPsiAODTracks = new TClonesArray("AliAODTrack", 1000);
  fJPsiESDTracks = new TClonesArray("AliESDtrack", 1000);
  fPsi2sAODTracks = new TClonesArray("AliAODTrack", 1000);
  fPsi2sESDTracks = new TClonesArray("AliESDtrack", 1000);

  //output tree with JPsi candidate events
  fJPsiTree = new TTree("fJPsiTree", "fJPsiTree");
  fJPsiTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fJPsiTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fJPsiTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fJPsiTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fJPsiTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fJPsiTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fJPsiTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fJPsiTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fJPsiTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fJPsiTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fJPsiTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fJPsiTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fJPsiTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fJPsiTree ->Branch("fDataFilnam", &fDataFilnam);
  fJPsiTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fJPsiTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L"); 		       
  if( fType == 0 ) {
    fJPsiTree ->Branch("fJPsiESDTracks", &fJPsiESDTracks);
  }
  if( fType == 1 ) {
    fJPsiTree ->Branch("fJPsiAODTracks", &fJPsiAODTracks);
  }
 
 //output tree with Psi2s candidate events
  fPsi2sTree = new TTree("fPsi2sTree", "fPsi2sTree");
  fPsi2sTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fPsi2sTree ->Branch("fPerNum", &fPerNum, "fPerNum/i");
  fPsi2sTree ->Branch("fOrbNum", &fOrbNum, "fOrbNum/i");
  
  fPsi2sTree ->Branch("fBCrossNum", &fBCrossNum, "fBCrossNum/s");
  fPsi2sTree ->Branch("fTrigger", &fTrigger[0], Form("fTrigger[%i]/O", ntrg));
  fPsi2sTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fPsi2sTree ->Branch("fL1inputs", &fL1inputs, "fL1inputs/i");
  fPsi2sTree ->Branch("fNtracklets", &fNtracklets, "fNtracklets/s");
  fPsi2sTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  fPsi2sTree ->Branch("fZDCAenergy", &fZDCAenergy, "fZDCAenergy/D");
  fPsi2sTree ->Branch("fZDCCenergy", &fZDCCenergy, "fZDCCenergy/D");
  fPsi2sTree ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fPsi2sTree ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");  
  fPsi2sTree ->Branch("fDataFilnam", &fDataFilnam);
  fPsi2sTree ->Branch("fRecoPass", &fRecoPass, "fRecoPass/S");
  fPsi2sTree ->Branch("fEvtNum", &fEvtNum, "fEvtNum/L");  		       
  if( fType == 0 ) {
    fPsi2sTree ->Branch("fPsi2sESDTracks", &fPsi2sESDTracks);
  }
  if( fType == 1 ) {
    fPsi2sTree ->Branch("fPsi2sAODTracks", &fPsi2sAODTracks);
  }
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);
  PostData(3, hCounter);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ) RunESD();
  if( fType == 1 ) RunAOD();

}//UserExec

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunAOD()
{
  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  //input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  fDataFilnam->Clear();
  fDataFilnam->SetString(filnam);
  fEvtNum = ((TTree*) GetInputData(0))->GetTree()->GetReadEntry();
  fRunNum = aod ->GetRunNumber();

  hCounter->Fill( 1 );

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  fTrigger[0]   = trigger.Contains("CINT7-B");
  fTrigger[1]   = trigger.Contains("CCUP4-B"); // CE 
  fTrigger[2]   = trigger.Contains("CCUP4-E"); // CE 
  /*/
  fTrigger[0]  = trigger.Contains("CCUP6-B"); // CE
  fTrigger[1]  = trigger.Contains("CCUP6-ACE");
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // CE
  fTrigger[3]  = trigger.Contains("CCUP7-ACE");/*/

  Bool_t isTRG = kFALSE;
  for(Int_t i=1; i<ntrg; i++) {
    if( fTrigger[i] ) {isTRG = kTRUE; hCounter->Fill( fRunNum - 167806 + 1 + i*2000 );}
  }
  if( !isTRG ) {PostData(3, hCounter); return;}

  hCounter->Fill( 2 );

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

   //Int_t nLepton=0;
  // Int_t nPion=0;
  // Int_t nHighPt=0;
  
   Int_t nGoodTracks=0;
  //Two tracks loop
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};

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
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
  
  if(nGoodTracks == 2){
  	  for(Int_t i=0; i<2; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		trk->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov);
				
		trk->SetDCA(dca[0],dca[1]); //to get DCAxy trk->DCA(); to get DCAz trk->ZatDCA();
		new((*fJPsiAODTracks)[i]) AliAODTrack(*trk); 
		
  		}
  fJPsiTree ->Fill();
  PostData(1, fJPsiTree);
  }
  
  if(nGoodTracks == 4){
  	  for(Int_t i=0; i<4; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		trk->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov);
		
		trk->SetDCA(dca[0],dca[1]); //to get DCAxy trk->DCA(); to get DCAz trk->ZatDCA();
		new((*fPsi2sAODTracks)[i]) AliAODTrack(*trk); 
		
  		}
  fPsi2sTree ->Fill();
  PostData(2, fPsi2sTree);
  }
    
  PostData(3, hCounter);

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunESD()
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

  hCounter->Fill( 1 );

  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  fTrigger[0]   = trigger.Contains("CINT7-B");
  fTrigger[1]   = trigger.Contains("CCUP4-B"); // CE 
  fTrigger[2]   = trigger.Contains("CCUP4-E"); // CE 
 /*/
  fTrigger[0]  = trigger.Contains("CCUP6-B"); // CE
  fTrigger[1]  = trigger.Contains("CCUP6-ACE");
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // CE
  fTrigger[3]  = trigger.Contains("CCUP7-ACE");/*/


  Bool_t isTRG = kFALSE;
  for(Int_t i=1; i<ntrg; i++) {
    if( fTrigger[i] ) {isTRG = kTRUE; hCounter->Fill( fRunNum - 167806 + 1 + i*2000 );}
  }
  if( !isTRG ) {PostData(3, hCounter); return;}

  hCounter->Fill( 2 );

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

  //Tracklets
  fNtracklets = esd->GetMultiplicity()->GetNumberOfTracklets();

  //VZERO, ZDC
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  fZDCAenergy = fZDCdata->GetZN2TowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZN1TowerEnergy()[0];


   //Int_t nLepton=0;
  // Int_t nPion=0;
  // Int_t nHighPt=0;
  
   Int_t nGoodTracks=0;
  //Track loop
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
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
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 4) break;   
  }//Track loop

  if(nGoodTracks == 2){
  	  for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);
				
		new((*fJPsiESDTracks)[i]) AliESDtrack(*trk); 
		
  		}
  fJPsiTree ->Fill();
  PostData(1, fJPsiTree);
  }
  
  if(nGoodTracks == 4){
  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fPsi2sESDTracks)[i]) AliESDtrack(*trk); 
		
  		}
  fPsi2sTree ->Fill();
  PostData(2, fPsi2sTree);
  }
    
  PostData(3, hCounter);

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate






























