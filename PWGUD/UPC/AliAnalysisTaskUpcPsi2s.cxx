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

// my headers
#include "AliAnalysisTaskUpcPsi2s.h"

ClassImp(AliAnalysisTaskUpcPsi2s);

using std::cout;
using std::endl;

//trees for UPC analysis,
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s() 
  : AliAnalysisTaskSE(),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFtrig1(0), fTOFtrig2(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),
    fListTrig(0),fHistUpcTriggersPerRun(0),fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0),
    fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fListHist(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcPsi2s


//_____________________________________________________________________________
AliAnalysisTaskUpcPsi2s::AliAnalysisTaskUpcPsi2s(const char *name) 
  : AliAnalysisTaskSE(name),fType(0),fRunTree(kTRUE),fRunHist(kTRUE),fJPsiTree(0),fPsi2sTree(0),
    fRunNum(0),fPerNum(0),fOrbNum(0),fL0inputs(0),fL1inputs(0),
    fTOFtrig1(0), fTOFtrig2(0),
    fVtxContrib(0),fVtxChi2(0),fVtxNDF(0),
    fBCrossNum(0),fNtracklets(0),fNLooseTracks(0),
    fZDCAenergy(0),fZDCCenergy(0),fV0Adecision(0),fV0Cdecision(0),
    fDataFilnam(0),fRecoPass(0),fEvtNum(0),
    fJPsiAODTracks(0),fJPsiESDTracks(0),fPsi2sAODTracks(0),fPsi2sESDTracks(0),
    fListTrig(0),fHistUpcTriggersPerRun(0),fHistZedTriggersPerRun(0),fHistCvlnTriggersPerRun(0),
    fHistMBTriggersPerRun(0),fHistCentralTriggersPerRun(0),fHistSemiCentralTriggersPerRun(0),
    fListHist(0),fHistNeventsJPsi(0),fHistTPCsignalJPsi(0),fHistDiLeptonPtJPsi(0),fHistDiElectronMass(0),fHistDiMuonMass(0),
    fHistNeventsPsi2s(0),fHistPsi2sMassVsPt(0),fHistPsi2sMassCoherent(0)

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
  
  for(Int_t i=0; i<ntrg; i++) fTrigger[i] = kFALSE;
  for(Int_t i=0; i<4; i++) fTOFphi[i] = -666;
  for(Int_t i=0; i<3; i++){
  	fVtxPos[i] = -666; 
	fVtxErr[i] = -666;
	fKfVtxPos[i] = -666;
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
  fJPsiTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fJPsiTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  
  fJPsiTree ->Branch("fTOFtrig1", &fTOFtrig1, "fTOFtrig1/O");
  fJPsiTree ->Branch("fTOFtrig2", &fTOFtrig2, "fTOFtrig2/O");
  fJPsiTree ->Branch("fTOFphi", &fTOFphi[0], "fTOFphi[4]/D");
  
  fJPsiTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fJPsiTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fJPsiTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fJPsiTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fJPsiTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  
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
  fPsi2sTree ->Branch("fNLooseTracks", &fNLooseTracks, "fNLooseTracks/s");
  fPsi2sTree ->Branch("fVtxContrib", &fVtxContrib, "fVtxContrib/I");
  
  fPsi2sTree ->Branch("fTOFtrig1", &fTOFtrig1, "fTOFtrig1/O");
  fPsi2sTree ->Branch("fTOFtrig2", &fTOFtrig2, "fTOFtrig2/O");
  fPsi2sTree ->Branch("fTOFphi", &fTOFphi[0], "fTOFphi[4]/D");
  
  fPsi2sTree ->Branch("fVtxPos", &fVtxPos[0], "fVtxPos[3]/D");
  fPsi2sTree ->Branch("fVtxErr", &fVtxErr[0], "fVtxErr[3]/D");
  fPsi2sTree ->Branch("fVtxChi2", &fVtxChi2, "fVtxChi2/D");
  fPsi2sTree ->Branch("fVtxNDF", &fVtxNDF, "fVtxNDF/D");
  
  fPsi2sTree ->Branch("fKfVtxPos", &fKfVtxPos[0], "fKfVtxPos[3]/D");
  
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
  
  fListTrig = new TList();
  fListTrig ->SetOwner();
  
  fHistUpcTriggersPerRun = new TH1D("fHistUpcTriggersPerRun", "fHistUpcTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistUpcTriggersPerRun);
  
  fHistZedTriggersPerRun = new TH1D("fHistZedTriggersPerRun", "fHistZedTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistZedTriggersPerRun);

  fHistCvlnTriggersPerRun = new TH1D("fHistCvlnTriggersPerRun", "fHistCvlnTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistCvlnTriggersPerRun);
  
  fHistMBTriggersPerRun = new TH1D("fHistMBTriggersPerRun", "fHistMBTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistMBTriggersPerRun);
  
  fHistCentralTriggersPerRun = new TH1D("fHistCentralTriggersPerRun", "fHistCentralTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistCentralTriggersPerRun);
  
  fHistSemiCentralTriggersPerRun = new TH1D("fHistSemiCentralTriggersPerRun", "fHistSemiCentralTriggersPerRun", 3000, 167000.5, 170000.5);
  fListTrig->Add(fHistSemiCentralTriggersPerRun);
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
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
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);
  PostData(3, fListTrig);
  PostData(4, fListHist);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::UserExec(Option_t *) 
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
void AliAnalysisTaskUpcPsi2s::RunAODtrig()
{

  //input event
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  fRunNum = aod ->GetRunNumber();
  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if(trigger.Contains("CCUP4-B")) fHistUpcTriggersPerRun->Fill(fRunNum); //Upc triggers
  
  if(trigger.Contains("CVLN_B2-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - synchronously downscaled
  if(trigger.Contains("CVLN_R1-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - randomly downscaled
  
  fL1inputs = aod->GetHeader()->GetL1TriggerInputs();
  if(fL1inputs & (1 << 18)) fHistZedTriggersPerRun->Fill(fRunNum); //1ZED trigger inputs
  
  //MB, Central and SemiCentral triggers
  AliCentrality *centrality = aod->GetCentrality();
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  Double_t percentile = centrality->GetCentralityPercentileUnchecked("V0M");
  //Double_t percentile = centrality->GetCentralityPercentile("V0M");
  
  if(((selectionMask & AliVEvent::kMB) == AliVEvent::kMB) && percentile<80 && percentile>0) fHistMBTriggersPerRun->Fill(fRunNum);
  
  if(((selectionMask & AliVEvent::kCentral) == AliVEvent::kCentral) && percentile<6 && percentile>0) fHistCentralTriggersPerRun->Fill(fRunNum);

  if(((selectionMask & AliVEvent::kSemiCentral) == AliVEvent::kSemiCentral) && percentile<50 && percentile>15) fHistSemiCentralTriggersPerRun->Fill(fRunNum);

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

  fHistNeventsJPsi->Fill(1);
  fHistNeventsPsi2s->Fill(1);

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if( !trigger.Contains("CCUP") ) return;
  
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
  AliAODZDC *fZDCdata = aod->GetZDCData();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;
  
  fZDCAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZNCTowerEnergy()[0];

  if( fZDCAenergy > 8200 || fZDCCenergy > 8200) return;
  
  fHistNeventsJPsi->Fill(4);
  fHistNeventsPsi2s->Fill(4);

  //Two tracks loop
  Int_t nGoodTracks = 0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4], qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0;
  Double_t fRecTPCsignal[5];
  Int_t mass[3]={-1,-1,-1};
  
   
  //Four track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;

      if(TMath::Abs(dca[1]) > 2) continue;
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
  
  nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;
  
  if(nGoodTracks == 4){
  	  fHistNeventsPsi2s->Fill(5);
  	  for(Int_t i=0; i<4; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
      		if(trk->Pt() > 1){   
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
  
  nGoodTracks = 0;
  //Two track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

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
      if(TMath::Abs(dca[0]) > 0.2) continue;
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
  
   nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;

  if(nGoodTracks == 2){
  	  fHistNeventsJPsi->Fill(5);
  	  for(Int_t i=0; i<2; i++){
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);		
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
	 	if(qLepton[0]*qLepton[1] > 0) fHistNeventsJPsi->Fill(6);
		if(qLepton[0]*qLepton[1] < 0){
			fHistNeventsJPsi->Fill(7);
			if(nHighPt > 0){
				fHistNeventsJPsi->Fill(8);
				fHistTPCsignalJPsi->Fill(fRecTPCsignal[0],fRecTPCsignal[1]);
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

  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  fTrigger[0]  = trigger.Contains("CCUP4-B"); // Central UPC Pb-Pb 2011
  fTrigger[1]  = trigger.Contains("CCUP2-B"); // Double gap
  fTrigger[2]  = trigger.Contains("CCUP7-B"); // Central UPC p-Pb 2013
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if( !isTriggered ) return;

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
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliAODTrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 20)continue;
      fNLooseTracks++; 
  }//Track loop -loose cuts
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  //Two track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

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
      if(TMath::Abs(dca[0]) > 0.2) continue;
     
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop

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
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      		if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      		delete trk_clone;
				
		new((*fJPsiAODTracks)[i]) AliAODTrack(*trk); 
		((AliAODTrack*)((*fJPsiAODTracks)[i]))->SetDCA(dca[0],dca[1]);//to get DCAxy trk->DCA(); to get DCAz trk->ZAtDCA();
		
		trk->GetPosition(KFpar);
    		trk->PxPyPz(KFpar+3);
    		trk->GetCovMatrix(KFcov);
		
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
		
		
		Double_t pos[3]={0,0,0};
    		AliExternalTrackParam *parTrk = new AliExternalTrackParam();
		parTrk->CopyFromVTrack((AliVTrack*) trk);
      		if(!parTrk->GetXYZAt(378,aod->GetMagneticField(),pos)) fTOFphi[i] = -666;
    		else {
		     fTOFphi[i] =  TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
    		     if(fTOFphi[i] < 0) fTOFphi[i]+=(2*TMath::Pi()*TMath::RadToDeg());
		     }
		delete parTrk;		
  		}
  fKfVtxPos[0]= KFvtx->GetX();
  fKfVtxPos[1]= KFvtx->GetY();
  fKfVtxPos[2]= KFvtx->GetZ();
  for(UInt_t i=0; i<2; i++)delete KFpart[i];
  delete KFvtx; 

  fJPsiTree ->Fill();
  }
  
   nGoodTracks = 0;
   //Four track loop
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = aod->GetTrack(itr);
    if( !trk ) continue;

      if(!(trk->GetStatus() & AliAODTrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliAODTrack::kITSrefit) ) continue;
      if(trk->GetTPCNcls() < 50)continue;
      if(trk->Chi2perNDF() > 4)continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      if(!trk->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      if(TMath::Abs(dca[1]) > 2) continue;

      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
      
    
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
	  	AliAODTrack *trk = aod->GetTrack(TrackIndex[i]);
		
		Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
		AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      		if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      		delete trk_clone;
		
		new((*fPsi2sAODTracks)[i]) AliAODTrack(*trk);
		((AliAODTrack*)((*fPsi2sAODTracks)[i]))->SetDCA(dca[0],dca[1]);//to get DCAxy trk->DCA(); to get DCAz trk->ZAtDCA();
				
		trk->GetPosition(KFpar);
    		trk->PxPyPz(KFpar+3);
    		trk->GetCovMatrix(KFcov);
		
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
				
		Double_t pos[3]={0,0,0};
    		AliExternalTrackParam *parTrk = new AliExternalTrackParam();
		parTrk->CopyFromVTrack((AliVTrack*) trk);
      		if(!parTrk->GetXYZAt(378,aod->GetMagneticField(),pos)) fTOFphi[i] = -666;
    		else {
		     fTOFphi[i] =  TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
    		     if(fTOFphi[i] < 0) fTOFphi[i]+=(2*TMath::Pi()*TMath::RadToDeg());
		     }
		delete parTrk; 		
  		}
  fKfVtxPos[0]= KFvtx->GetX();
  fKfVtxPos[1]= KFvtx->GetY();
  fKfVtxPos[2]= KFvtx->GetZ();
  for(UInt_t i=0; i<4; i++)delete KFpart[i];
  delete KFvtx; 
  fPsi2sTree ->Fill();
  }
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunESDtrig()
{

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  fRunNum = esd ->GetRunNumber();
  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  if(trigger.Contains("CCUP4-B")) fHistUpcTriggersPerRun->Fill(fRunNum); //Upc triggers
  
  if(trigger.Contains("CVLN_B2-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - synchronously downscaled
  if(trigger.Contains("CVLN_R1-B")) fHistCvlnTriggersPerRun->Fill(fRunNum); //CVLN triggers - randomly downscaled
  
  if(esd->GetHeader()->IsTriggerInputFired("1ZED")) fHistZedTriggersPerRun->Fill(fRunNum); //1ZED trigger inputs
  
   //MB, Central and SemiCentral triggers
  AliCentrality *centrality = esd->GetCentrality();
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  //Double_t percentile = centrality->GetCentralityPercentile("V0M");
  Double_t percentile = centrality->GetCentralityPercentileUnchecked("V0M");
  
  if(((selectionMask & AliVEvent::kMB) == AliVEvent::kMB) && percentile<80 && percentile>0) fHistMBTriggersPerRun->Fill(fRunNum);
  
  if(((selectionMask & AliVEvent::kCentral) == AliVEvent::kCentral) && percentile<6 && percentile>0) fHistCentralTriggersPerRun->Fill(fRunNum);

  if(((selectionMask & AliVEvent::kSemiCentral) == AliVEvent::kSemiCentral) && percentile<50 && percentile>15) fHistSemiCentralTriggersPerRun->Fill(fRunNum);

  
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

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  fHistNeventsJPsi->Fill(1);
  fHistNeventsPsi2s->Fill(1);

  //Trigger
  TString trigger = esd->GetFiredTriggerClasses();
  
  if( !trigger.Contains("CCUP") ) return;
  
  fHistNeventsJPsi->Fill(2);
  fHistNeventsPsi2s->Fill(2);

  //primary vertex
  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
  fVtxContrib = fESDVertex->GetNContributors();
  if(fVtxContrib < 2) return;
  
  fHistNeventsJPsi->Fill(3);
  fHistNeventsPsi2s->Fill(3);

  //VZERO, ZDC
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  if(fV0Adecision != AliESDVZERO::kV0Empty || fV0Cdecision != AliESDVZERO::kV0Empty) return;
  
  fZDCAenergy = fZDCdata->GetZN2TowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZN1TowerEnergy()[0];
  if( fZDCAenergy > 8200 || fZDCCenergy > 8200) return;
  
  fHistNeventsJPsi->Fill(4);
  fHistNeventsPsi2s->Fill(4);

   Int_t nGoodTracks=0;
  //Two tracks loop
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4], qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0;
  Double_t fRecTPCsignal[5];
  Int_t mass[3]={-1,-1,-1};

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
      if(TMath::Abs(dca[1]) > 2) continue;
      if(TMath::Abs(dca[1]) > 0.2) continue;
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 2) break;   
  }//Track loop

  
  if(nGoodTracks == 2){
  	  fHistNeventsJPsi->Fill(5);
  	  for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);		
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
	 	if(qLepton[0]*qLepton[1] > 0) fHistNeventsJPsi->Fill(6);
		if(qLepton[0]*qLepton[1] < 0){
			fHistNeventsJPsi->Fill(7);
			if(nHighPt > 0){
				fHistNeventsJPsi->Fill(8);
				fHistTPCsignalJPsi->Fill(fRecTPCsignal[0],fRecTPCsignal[1]);
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
  nGoodTracks = 0; nLepton=0; nPion=0; nHighPt=0;
  mass[0]= -1; mass[1]= -1, mass[2]= -1;
  
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
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 4) break;   
  }//Track loop
  
  if(nGoodTracks == 4){
  	  fHistNeventsPsi2s->Fill(5);
  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
      		if(trk->Pt() > 1){   
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
  
  PostData(4, fListHist);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::RunESDtree()
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
  
  Bool_t isTriggered = kFALSE;
  for(Int_t i=0; i<ntrg; i++) {
    if( fTrigger[i] ) isTriggered = kTRUE;
  }
  if( !isTriggered ) return;
  
  //trigger inputs
  fL0inputs = esd->GetHeader()->GetL0TriggerInputs();
  fL1inputs = esd->GetHeader()->GetL1TriggerInputs();
  
  //TOF trigger info (0OMU)
  fTOFtrig1 = esd->GetHeader()->IsTriggerInputFired("0OMU");
  fTOFtrig2 = esd->GetHeader()->GetActiveTriggerInputs().Contains("0OMU") ? ((esd->GetHeader()) ? esd->GetHeader()->IsTriggerInputFired("0OMU") : kFALSE) : TESTBIT(esd->GetHeader()->GetL0TriggerInputs(), (kFALSE) ? 21 : 9);

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

  //Tracklets
  fNtracklets = esd->GetMultiplicity()->GetNumberOfTracklets();

  //VZERO, ZDC
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  fZDCAenergy = fZDCdata->GetZN2TowerEnergy()[0];
  fZDCCenergy = fZDCdata->GetZN1TowerEnergy()[0];

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
  
  Int_t nGoodTracks=0;
  Int_t TrackIndex[5] = {-1,-1,-1,-1,-1};
  
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
      if(TMath::Abs(dca[1]) > 2) continue;
      if(TMath::Abs(dca[1]) > 0.2) continue;
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 2) break;   
  }//Track loop

  if(nGoodTracks == 2){
  	  for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);
				
		new((*fJPsiESDTracks)[i]) AliESDtrack(*trk); 
		
		Double_t pos[3]={0,0,0};
    		if(!trk->GetXYZAt(378,esd->GetMagneticField(),pos)) fTOFphi[i] = -666;
    		else {
		     fTOFphi[i] =  TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
    		     if(fTOFphi[i] < 0) fTOFphi[i]+=(2*TMath::Pi()*TMath::RadToDeg());
		     } 		
  		}
  fJPsiTree ->Fill();
  }
  
  nGoodTracks = 0;
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
      
      TrackIndex[nGoodTracks] = itr;
      nGoodTracks++;
      if(nGoodTracks > 4) break;   
  }//Track loop
  
  if(nGoodTracks == 4){
  	  for(Int_t i=0; i<4; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
		
		AliExternalTrackParam cParam;
      		trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam);// to get trk->GetImpactParameters(DCAxy,DCAz);

		new((*fPsi2sESDTracks)[i]) AliESDtrack(*trk);
		 
		Double_t pos[3]={0,0,0};
    		if(!trk->GetXYZAt(378,esd->GetMagneticField(),pos)) fTOFphi[i] = -666;
    		else {
		     fTOFphi[i] =  TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
    		     if(fTOFphi[i] < 0) fTOFphi[i]+=(2*TMath::Pi()*TMath::RadToDeg());
		     } 		
  		}
  fPsi2sTree ->Fill();
  }
  
  PostData(1, fJPsiTree);
  PostData(2, fPsi2sTree);

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcPsi2s::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate

