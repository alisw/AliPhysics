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

// my headers
#include "AliAnalysisTaskUpcPsi2s.h"

ClassImp(AliAnalysisTaskUpcPsi2s);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s() 
  : AliAnalysisTaskSE(),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),hCounter(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),
    fListHist(0),fHistTriggersPerRun(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0),
    fHistK0sMass(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),hCounter(0),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),fVtxContrib(0),fBCrossNum(0),fNtracklets(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),
    fListHist(0),fHistTriggersPerRun(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0),
    fHistK0sMass(0)

{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;
  
  Init();

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TH1I::Class());
  DefineOutput(4, TList::Class());

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
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
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
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
  fHistTriggersPerRun = new TH1D("fHistTriggersPerRun", "fHistTriggersPerRun", 3000, 167000.5, 170000.5);
  fListHist->Add(fHistTriggersPerRun);
  
  TString CutNameJPsi[12] = {"Analyzed","Triggered","Vertex cut","V0 decision","Two good tracks",
  				"Like sign","Oposite sign","One p_{T}>1", "Both p_{T}>1","PID","Dimuom","Dielectron"};
  fHistNeventsJPsi = new TH1D("fHistNeventsJPsi","fHistNeventsPsi2s",12,0.5,12.5);
  for (Int_t i = 0; i<12; i++) fHistNeventsJPsi->GetXaxis()->SetBinLabel(i+1,CutNameJPsi[i].Data());
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

  TString CutNamePsi2s[13] = {"Analyzed","Triggered","Vertex cut","V0 decision","Four good tracks",
  				"DiLepton - DiPion","Like sign leptons","Like sign pions","Like sign both","Oposite sign","PID","Dimuom","Dielectron"};

  fHistNeventsPsi2s = new TH1D("fHistNeventsPsi2s","fHistNeventsPsi2s",13,0.5,13.5);
  for (Int_t i = 0; i<13; i++) fHistNeventsPsi2s->GetXaxis()->SetBinLabel(i+1,CutNamePsi2s[i].Data());
  fListHist->Add(fHistNeventsPsi2s);

  fHistPsi2sMassVsPt = new TH2D("fHistPsi2sMassVsPt","Mass vs p_{T} of #psi(2s) candidates",100,3,6,50,0,5);
  fHistPsi2sMassVsPt->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fHistPsi2sMassVsPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fListHist->Add(fHistPsi2sMassVsPt);
  
  fHistPsi2sMassCoherent = new TH1D("fHistPsi2sMassAllCoherent","Invariant mass of coherent #psi(2s) candidates",100,3,6);
  fHistPsi2sMassCoherent->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}#pi^{+}#pi^{-}) (GeV/c)");
  fListHist->Add(fHistPsi2sMassCoherent);
  
  fHistK0sMass = new TH1D("fHistK0sMass","fHistK0sMass",200,0.4,0.6);
  fListHist->Add(fHistK0sMass);
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);
  PostData(3, hCounter);
  PostData(4, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserExec(Option_t *) 
{

  //cout<<"#################### Next event ##################"<<endl;

  if( fType == 0 ) RunESD();
  if( fType == 1 ){ 
  	if(fRunHist) RunAODhist();
	if(fRunTree) RunAODtree();
	}

}//UserExec
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

  fHistNeventsJPsi->Fill(1);
  fHistNeventsPsi2s->Fill(1);

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if( !trigger.Contains("CCUP4-B") ) return;
  
  fRunNum = aod ->GetRunNumber();
  fHistTriggersPerRun->Fill(fRunNum);
  
  fHistNeventsJPsi->Fill(2);
  fHistNeventsPsi2s->Fill(2);

  //primary vertex
  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  fVtxContrib = fAODVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsJPsi->Fill(3);
  fHistNeventsPsi2s->Fill(3);


  //VZERO, ZDC
  AliAODVZERO *fV0data = aod ->GetVZEROData();
  //AliAODZDC *fZDCdata = aod->GetZDCData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;
  
  fHistNeventsJPsi->Fill(4);
  fHistNeventsPsi2s->Fill(4);

   Int_t nGoodTracks=0;
  //Two tracks loop
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4], qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0;
  Double_t jRecTPCsignal[5];
  Int_t mass[3]={-1,-1,-1};

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
  	  fHistNeventsJPsi->Fill(5);
  	  for(Int_t i=0; i<2; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);		
      		if(trk->Pt() > 1) nHighPt++;     
      		jRecTPCsignal[nLepton] = trk->GetTPCsignal();     
      		qLepton[nLepton] = trk->Charge();
      		if(jRecTPCsignal[nLepton] > 40 && jRecTPCsignal[nLepton] < 70){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				mass[nLepton] = 0;
				}
      		if(jRecTPCsignal[nLepton] > 70 && jRecTPCsignal[nLepton] < 100){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				mass[nLepton] = 1;
				}
       		nLepton++;		
  		}		
  	if(nLepton == 2){
	 	if(qLepton[0]*qLepton[1] > 0) fHistNeventsJPsi->Fill(6);
		if(qLepton[0]*qLepton[1] < 0){
			fHistNeventsJPsi->Fill(7);
			if(nHighPt > 0){
				fHistNeventsJPsi->Fill(8);
				fHistTPCsignalJPsi->Fill(jRecTPCsignal[0],jRecTPCsignal[1]);
				if(nHighPt == 2) fHistNeventsJPsi->Fill(9);
				if(mass[0] == mass[1] && mass[0] != -1) {
					fHistNeventsJPsi->Fill(10);
					vCandidate = vLepton[0]+vLepton[1];
					if( vCandidate.M() > 2.8 && vCandidate.M() < 3.2) fHistDiLeptonPtJPsi->Fill(vLepton[0].Pt(),vLepton[1].Pt());
					if(mass[0] == 0) {
						fHistDiMuonMass->Fill(vCandidate.M());
						fHistNeventsJPsi->Fill(11);
						}
  					if(mass[0] == 1) {
						fHistDiElectronMass->Fill(vCandidate.M());
						fHistNeventsJPsi->Fill(12);
   						}
					}
				}
			}
		}
  }
  nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;
  
  if(nGoodTracks == 4){
  	  fHistNeventsPsi2s->Fill(5);
  	  for(Int_t i=0; i<4; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
      		if(trk->Pt() > 1){   
      			jRecTPCsignal[nLepton] = trk->GetTPCsignal();      
      			qLepton[nLepton] = trk->Charge();
      			if(jRecTPCsignal[nLepton] > 40 && jRecTPCsignal[nLepton] < 70){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
					mass[nLepton] = 0;
					}
      			if(jRecTPCsignal[nLepton] > 70 && jRecTPCsignal[nLepton] < 100){
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
		fHistNeventsPsi2s->Fill(6);
		if(qLepton[0]*qLepton[1] > 0) fHistNeventsPsi2s->Fill(7);
		if(qPion[0]*qPion[1] > 0) fHistNeventsPsi2s->Fill(8);
		if((qLepton[0]*qLepton[1] > 0) && (qPion[0]*qPion[1] > 0)) fHistNeventsPsi2s->Fill(9);
		if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0)){
			fHistNeventsPsi2s->Fill(10);
	 		if(mass[0] == mass[1]) {
				fHistNeventsPsi2s->Fill(11); 
  				vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  				vDilepton = vLepton[0]+vLepton[1];
				fHistPsi2sMassVsPt->Fill(vCandidate.M(),vCandidate.Pt());
				if(vCandidate.Pt() < 0.15) fHistPsi2sMassCoherent->Fill(vCandidate.M());
				if(mass[0] == 0) fHistNeventsPsi2s->Fill(12);	
  				if(mass[0] == 1) fHistNeventsPsi2s->Fill(13);
				}
			}
		}
  }
  
  
   //---------------------------------K0s + K0s loop - very experimental-------------------- 
   nGoodTracks = 0;
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
      if(pTrack->GetTPCNcls() < 50)continue;
      if(nTrack->GetTPCNcls() < 50)continue;
      if(pTrack->Chi2perNDF() > 4)continue;
      if(nTrack->Chi2perNDF() > 4)continue;
      
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      if(!pTrack->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;
      if(!nTrack->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;
      
      TrackIndex[nGoodTracks] = iV0;
      nGoodTracks++; 
      if(nGoodTracks > 2) break;
  }//V0s loop
  if(nGoodTracks == 2){
  	for(Int_t i=0; i<2; i++){
	  	AliAODv0 *v0 = aod->GetV0(TrackIndex[i]);
  		fHistK0sMass->Fill(v0->MassK0Short());
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
  
  Int_t nGoodTracks=0;
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
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  //Track loop
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






























