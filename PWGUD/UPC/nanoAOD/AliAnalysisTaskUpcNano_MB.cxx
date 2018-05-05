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
#include "TROOT.h"
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TColor.h"
#include "TRandom.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliVVZERO.h"
#include "AliAODZDC.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliTOFTriggerMask.h"
#include "TObjArray.h"
#include "AliDataFile.h"

#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"


// my headers
#include "AliAnalysisTaskUpcNano_MB.h"

ClassImp(AliAnalysisTaskUpcNano_MB);

using std::cout;
using std::endl;

//analysis skeleton of UPC nano AODs,

//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::AliAnalysisTaskUpcNano_MB() 
  : AliAnalysisTaskSE(),fPIDResponse(0), fTrackCutsBit0(0),fTrackCutsBit1(0),fTrackCutsBit5(0),isMC(kFALSE), isESD(kFALSE),cutEta(0.9),fOutputList(0),
    	fHistEvents(0),
	fHistMCTriggers(0),
    	fTreePhi(0),
	fTreeJPsi(0),
	fTreePsi2s(0),
	fTreeRho(0),
	fTreeGen(0),
	hTPCPIDMuonCorr(0),
    	hTPCPIDMuon(0),
        hTPCPIDElectronCorr(0),
        hTPCPIDElectron(0),
	hTPCPIDPion(0),
	hTPCPIDPionCorr(0),
	hTOFPIDProton(0),
	hTOFPIDProtonCorr(0),
	hITSPIDKaon(0),
	hITSPIDKaonCorr(0),
	hTPCdEdxCorr(0),
	fSPDfile(0),
  	hBCmod4(0),
  	hSPDeff(0) 
	

{

//Dummy constructor

}//AliAnalysisTaskUpcNano_MB


//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::AliAnalysisTaskUpcNano_MB(const char *name) 
  : AliAnalysisTaskSE(name),fPIDResponse(0),fTrackCutsBit0(0),fTrackCutsBit1(0),fTrackCutsBit5(0),isMC(kFALSE), isESD(kFALSE),cutEta(0.9),fOutputList(0),
    	fHistEvents(0),
	fHistMCTriggers(0),
    	fTreePhi(0),
	fTreeJPsi(0),
	fTreePsi2s(0),
	fTreeRho(0),
	fTreeGen(0),
	hTPCPIDMuonCorr(0),
    	hTPCPIDMuon(0),
        hTPCPIDElectronCorr(0),
        hTPCPIDElectron(0),
	hTPCPIDPion(0),
	hTPCPIDPionCorr(0),
	hTOFPIDProton(0),
	hTOFPIDProtonCorr(0),
	hITSPIDKaon(0),
	hITSPIDKaonCorr(0),
	hTPCdEdxCorr(0),
	fSPDfile(0),
  	hBCmod4(0),
  	hSPDeff(0) 

{
  for(Int_t i = 0; i<10; i++) fTriggerInputsMC[i] = kFALSE;
  DefineOutput(1, TList::Class());

}//AliAnalysisTaskUpcNano_MB

//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::~AliAnalysisTaskUpcNano_MB() 
{
  // Destructor
  
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;
     fOutputList = 0x0;
  }

}//~AliAnalysisTaskUpcNano_MB


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::UserCreateOutputObjects()
{
  
AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTrackCutsBit0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackCutsBit5 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  
  fTrackCutsBit1 = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
  fTrackCutsBit1->SetRequireITSStandAlone(kTRUE);
  
  fOutputList = new TList();
  fOutputList ->SetOwner();

  TString gCutName[9] = {"Analyzed","2 TPC"," 2TPC+2ITS","4TPC","3TPC+1ITS","4TPC+ITS","2TPC+1ITS","2ITS","4ITS"};		  
  fHistEvents = new TH1D("fHistEvents","NEvents",9,0.5,9.5);
  for (Int_t i = 0; i<9; i++) fHistEvents->GetXaxis()->SetBinLabel(i+1,gCutName[i].Data());
  fOutputList->Add(fHistEvents);
  
  TString gTriggerName[12] = {"0VBA","0VBC","0UBA","0UBC","0OMU","0OM2","0SMB","0SM2","0STP","0SH1","Offline OM2","Offline OMU"};		  
  fHistMCTriggers = new TH1D("fHistMCTriggers","Fired MC trigger inputs",12,0.5,12.5);
  for (Int_t i = 0; i<12; i++) fHistMCTriggers->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
  fOutputList->Add(fHistMCTriggers);
  
  fTreeJPsi = new TTree("fTreeJPsi", "fTreeJPsi");
  fTreeJPsi ->Branch("fPt", &fPt, "fPt/F");
  fTreeJPsi ->Branch("fY", &fY, "fY/F");
  fTreeJPsi ->Branch("fM", &fM, "fM/F");
  fTreeJPsi ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreeJPsi ->Branch("fSign", &fSign, "fSign/I");
  fTreeJPsi ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreeJPsi ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreeJPsi ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreeJPsi ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreeJPsi ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreeJPsi ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreeJPsi ->Branch("fNFiredMaxiPads", &fNFiredMaxiPads, "fNFiredMaxiPads/I");
  fTreeJPsi ->Branch("fNTOFtrgPads", &fNTOFtrgPads, "fNTOFtrgPads/I");
  fTreeJPsi ->Branch("fInEtaRec", &fInEtaRec, "fInEtaRec/O");
  fTreeJPsi ->Branch("fTrackIndex", &fTrackIndex[0],"fTrackIndex[2]/I");
  fTreeJPsi ->Branch("fTOFhits", &fTOFhits);
  fTreeJPsi ->Branch("fTrackIndices", &fTrackIndices);
  if(isMC){
	fTreeJPsi ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[10]/O");
	}
  fOutputList->Add(fTreeJPsi);
  
  fTreePhi = new TTree("fTreePhi", "fTreePhi");
  fTreePhi ->Branch("fPt", &fPt, "fPt/F");
  fTreePhi ->Branch("fY", &fY, "fY/F");
  fTreePhi ->Branch("fM", &fM, "fM/F");
  fTreePhi ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreePhi ->Branch("fSign", &fSign, "fSign/I");
  fTreePhi ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreePhi ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreePhi ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreePhi ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreePhi ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreePhi ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  if(isMC) fTreePhi ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[10]/O");
  fOutputList->Add(fTreePhi);
  
  fTreeRho = new TTree("fTreeRho", "fTreeRho");
  fTreeRho ->Branch("fPt", &fPt, "fPt/F");
  fTreeRho ->Branch("fY", &fY, "fY/F");
  fTreeRho ->Branch("fM", &fM, "fM/F");
  fTreeRho ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreeRho ->Branch("fSign", &fSign, "fSign/I");
  fTreeRho ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreeRho ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreeRho ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreeRho ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreeRho ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreeRho ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreeRho ->Branch("fNFiredMaxiPads", &fNFiredMaxiPads, "fNFiredMaxiPads/I");
  fTreeRho ->Branch("fNTOFtrgPads", &fNTOFtrgPads, "fNTOFtrgPads/I");
  fTreeRho ->Branch("fTrackIndex", &fTrackIndex[0],"fTrackIndex[2]/I");
  fTreeRho ->Branch("fTOFhits", &fTOFhits);
  fTreeRho ->Branch("fTrackIndices", &fTrackIndices);
  fTreeRho ->Branch("fInEtaRec", &fInEtaRec, "fInEtaRec/O");
  
  if(isMC) fTreeRho ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[10]/O");
  fOutputList->Add(fTreeRho);

  fTreePsi2s = new TTree("fTreePsi2s", "fTreePsi2s");
  fTreePsi2s ->Branch("fPt", &fPt, "fPt/F");
  fTreePsi2s ->Branch("fY", &fY, "fY/F");
  fTreePsi2s ->Branch("fM", &fM, "fM/F");
  fTreePsi2s ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreePsi2s ->Branch("fSign", &fSign, "fSign/I");
  fTreePsi2s ->Branch("fDiLeptonM", &fDiLeptonM, "fDiLeptonM/F");
  fTreePsi2s ->Branch("fDiLeptonPt", &fDiLeptonPt, "fDiLeptonPt/F");
  fTreePsi2s ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreePsi2s ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreePsi2s ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreePsi2s ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreePsi2s ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreePsi2s ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreePsi2s ->Branch("fNFiredMaxiPads", &fNFiredMaxiPads, "fNFiredMaxiPads/I");
  if(isMC){ 
  	fTreePsi2s ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[10]/O");
	}
  fOutputList->Add(fTreePsi2s);
    
  fTreeGen = new TTree("fTreeGen", "fTreeGen");
  fTreeGen ->Branch("fPt", &fPt, "fPt/F");
  fTreeGen ->Branch("fY", &fY, "fY/F");
  fTreeGen ->Branch("fM", &fM, "fM/F");
  fTreeGen ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreeGen ->Branch("fInEtaGen", &fInEtaGen, "fInEtaGen/O");
  if(isMC) fOutputList->Add(fTreeGen);
   
  hTPCPIDMuonCorr = new TH2D("hTPCPIDMuonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDMuonCorr->GetXaxis()->SetTitle("TPC PID N #sigma (#mu)");
  hTPCPIDMuonCorr->GetYaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hTPCPIDMuonCorr);
  
  hTPCPIDElectronCorr = new TH2D("hTPCPIDElectronCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDElectronCorr->GetXaxis()->SetTitle("TPC PID N #sigma (e)");
  hTPCPIDElectronCorr->GetYaxis()->SetTitle("TPC PID N #sigma (e)");
  fOutputList->Add(hTPCPIDElectronCorr);

  hTPCPIDPionCorr = new TH2D("hTPCPIDPionCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDPionCorr->GetXaxis()->SetTitle("TPC PID N #sigma (#pi)");
  hTPCPIDPionCorr->GetYaxis()->SetTitle("TPC PID N #sigma (#pi)");
  fOutputList->Add(hTPCPIDPionCorr);
  
  hTOFPIDProtonCorr = new TH2D("hTOFPIDProtonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTOFPIDProtonCorr->GetXaxis()->SetTitle("TOF PID N #sigma (p)");
  hTOFPIDProtonCorr->GetYaxis()->SetTitle("TOF PID N #sigma (p)");
  fOutputList->Add(hTOFPIDProtonCorr);
  
  hITSPIDKaonCorr = new TH2D("hITSPIDKaonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hITSPIDKaonCorr->GetXaxis()->SetTitle("ITS PID N #sigma (K)");
  hITSPIDKaonCorr->GetYaxis()->SetTitle("ITS PID N #sigma (K)");
  fOutputList->Add(hITSPIDKaonCorr);
  
  hTPCPIDMuon = new TH1D("hTPCPIDMuon"," ",200,-10,10);
  hTPCPIDMuon->GetXaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hTPCPIDMuon);
  
  hTPCPIDElectron = new TH1D("hTPCPIDElectron"," ",200,-10,10);
  hTPCPIDElectron->GetXaxis()->SetTitle("TPC PID N #sigma (#e)");
  fOutputList->Add(hTPCPIDElectron);
  
  hTPCPIDPion = new TH1D("hTPCPIDPion"," ",200,-10,10);
  hTPCPIDPion->GetXaxis()->SetTitle("TPC PID N #sigma (#pi)");
  fOutputList->Add(hTPCPIDPion);
  
  hTOFPIDProton = new TH1D("hTOFPIDProton"," ",200,-10,10);
  hTOFPIDProton->GetXaxis()->SetTitle("TOF PID N #sigma (p)");
  fOutputList->Add(hTOFPIDProton);
  
  hITSPIDKaon = new TH1D("hITSPIDKaon"," ",200,-10,10);
  hITSPIDKaon->GetXaxis()->SetTitle("ITS PID N #sigma (K)");
  fOutputList->Add(hITSPIDKaon);
  
  hTPCdEdxCorr = new TH2D("hTPCdEdxCorr"," ",200,0,200,200,0,200);
  hTPCdEdxCorr->GetXaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  hTPCdEdxCorr->GetYaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  fOutputList->Add(hTPCdEdxCorr);
  
  fSPDfile = AliDataFile::OpenOADB("PWGUD/UPC/SPDFOEfficiency_run245067.root");
  fSPDfile->Print();
  fSPDfile->Map();
  hSPDeff = (TH2D*) fSPDfile->Get("hEff");
  hSPDeff->SetDirectory(0);
  TH2D *hBCmod4_2D = (TH2D*) fSPDfile->Get("hCounts");
  hBCmod4 = hBCmod4_2D->ProjectionY();
  fSPDfile->Close();
 
    
  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::UserExec(Option_t *) 
{

  AliVEvent *fEvent = InputEvent();
  if(!fEvent) return;
  
  TString trigger = fEvent->GetFiredTriggerClasses();
  if(!isMC && !trigger.Contains("CCUP8-B") && !trigger.Contains("CCUP9-B"))return;
   
  fRunNumber = fEvent->GetRunNumber();
  
  AliVVZERO *fV0data = fEvent->GetVZEROData();
  AliVAD *fADdata = fEvent->GetADData();
  
  if(isMC) RunMC(fEvent);
  
  for(Int_t i = 0; i<10; i++)if(fTriggerInputsMC[i])fHistMCTriggers->Fill(i+1);
  
  if(isESD){
  	AliESDZDC *fZDCdata = (AliESDZDC*)fEvent->GetZDCData();
  	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];

	Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
        Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
	if (fEvent->GetRunNumber()>=245726 && fEvent->GetRunNumber()<=245793) detChZNA = 10;
  	for (Int_t i=0;i<4;i++){ 
  		fZNAtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
  		fZNCtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
		}
  }
  else{
  	AliAODZDC *fZDCdata = (AliAODZDC*)fEvent->GetZDCData();
  	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];

  	for (Int_t i=0;i<4;i++){ 
  		fZNAtime[i] = fZDCdata->GetZNATDCm(i);
  		fZNCtime[i] = fZDCdata->GetZNCTDCm(i);
		}
  } 
 
  Int_t fV0Adecision = fV0data->GetV0ADecision();
  Int_t fV0Cdecision = fV0data->GetV0CDecision();
  if( fV0Adecision != 0 || fV0Cdecision != 0) return;
  
  Int_t fADAdecision = fADdata->GetADADecision();
  Int_t fADCdecision = fADdata->GetADCDecision();
  if( fADAdecision != 0 || fADCdecision != 0) return;
  
  const AliTOFHeader *tofH = fEvent->GetTOFHeader();
  AliTOFTriggerMask *fTOFmask = tofH->GetTriggerMask();
  fNFiredMaxiPads = fTOFmask->GetNumberMaxiPadOn();
  fNTOFtrgPads = tofH->GetNumberOfTOFtrgPads();
  
  if(isESD){
  
  TClonesArray* tofClusters = ((AliESDEvent*)fInputEvent)->GetESDTOFClusters();

  Int_t nTOFhits = 0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
     AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
     nTOFhits+=cl->GetNTOFhits();
   }

   fTOFhits.Reset();
   fTrackIndices.Reset();
   fTOFhits.Set(nTOFhits);
   fTrackIndices.Set(nTOFhits);

   Int_t hitCounts=0;
   for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
     AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
     for (Int_t ihit=0;ihit<cl->GetNTOFhits();ihit++){
       AliESDTOFHit* hit = (AliESDTOFHit*) cl->GetTOFHit(ihit);
       Int_t channel = hit->GetTOFchannel();
       Int_t trackIndex = (cl->GetNMatchableTracks()==1) ? cl->GetTrackIndex(0) : -1;
       fTOFhits.AddAt(channel,hitCounts);
       fTrackIndices.AddAt(trackIndex,hitCounts);
       hitCounts++;
     }
   }
  
  }
    
  fHistEvents->Fill(1);
  
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Float_t kaonMass = partKaon->Mass();
  
  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  Float_t muonMass = partMuon->Mass();
  
  TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
  Float_t electronMass = partElectron->Mass();
  
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Float_t pionMass = partPion->Mass();
  
  TParticlePDG *partProton = pdgdat->GetParticle( 2212 );
  Float_t protonMass = partProton->Mass();
  
  Short_t qTrack[5];
  TLorentzVector vMuon[5],vElectron[5],vProton[5], vJPsiCandidate;
  TLorentzVector vKaon[5], vPhiCandidate;
  TLorentzVector vPion[5], vRhoCandidate;

  Float_t nSigmaMuon[5], nSigmaElectron[5], nSigmaPion[5], nSigmaProton[5],  nSigmaKaon[5], dEdx[5];
  Short_t qPion[5];
  TLorentzVector vLepton[5], vDilepton, vPsi2sCandidate;
  Short_t qLepton[5];
  UInt_t nPion = 0, nLepton = 0, nHighPt = 0;
  UInt_t nGoodTracksTPC=0;
  UInt_t nGoodTracksITS=0;
  UInt_t nGoodTracksSPD=0;
  UInt_t nGoodTracksLoose=0;
  Int_t TrackIndexTPC[5] = {-1,-1,-1,-1,-1};
  Int_t TrackIndexITS[5] = {-1,-1,-1,-1,-1};
  Int_t TrackIndexALL[7] = {-1,-1,-1,-1,-1,-1,-1};
  Double_t TrackPtTPC[5]={0,0,0,0,0};
  Double_t TrackPtALL[7]={0,0,0,0,0,0,0};
  Double_t MeanPt = -1;
   
  //Track loop
  for(Int_t iTrack=0; iTrack<fEvent ->GetNumberOfTracks(); iTrack++) {
  Bool_t goodTPCTrack = kTRUE;
  Bool_t goodITSTrack = kTRUE;
    if(isESD){ 
    	AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrack));
	if( !trk ) continue;
	
	if(fTrackCutsBit0->AcceptTrack(trk))nGoodTracksLoose++;
    	if(!fTrackCutsBit5->AcceptTrack(trk))goodTPCTrack = kFALSE;
    	else{
    		if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1) && trk->GetTOFsignal()<99998)nGoodTracksSPD++;
    		}
    
    	if(!fTrackCutsBit1->AcceptTrack(trk)) goodITSTrack = kFALSE;
	}
    else{ 
    	AliAODTrack *trk = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));
    	if( !trk ) continue;
    
    	if(trk->TestFilterBit(1<<0) && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1)))nGoodTracksLoose++;
    	if(!(trk->TestFilterBit(1<<5)))goodTPCTrack = kFALSE;
    	else{
    		if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1) && trk->GetTOFsignal()<99998)nGoodTracksSPD++;
    		}
    
    	if(!(trk->TestFilterBit(1<<1))) goodITSTrack = kFALSE;
	}
    AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(iTrack));

    if(goodTPCTrack){
    	TrackIndexTPC[nGoodTracksTPC] = iTrack;
    	TrackPtTPC[nGoodTracksTPC] = trk->Pt();
	TrackIndexALL[nGoodTracksTPC+nGoodTracksITS] = iTrack;
    	TrackPtALL[nGoodTracksTPC+nGoodTracksITS] = trk->Pt();
    	nGoodTracksTPC++;
	}
     	
    if(goodITSTrack){
    	TrackIndexITS[nGoodTracksITS] = iTrack;
	TrackIndexALL[nGoodTracksTPC+nGoodTracksITS] = iTrack;
    	TrackPtALL[nGoodTracksTPC+nGoodTracksITS] = trk->Pt();
    	nGoodTracksITS++;
    	}
     
	
    if(nGoodTracksTPC + nGoodTracksITS > 6) break;
    //if(nGoodTracksTPC > 4) break;
    }
    
  //{"Analyzed","2 TPC"," 2TPC+2ITS","4TPC","3TPC+1ITS","4TPC+ITS","2TPC+1ITS","2ITS","4ITS"};
  if(nGoodTracksTPC == 2 && nGoodTracksITS == 0)fHistEvents->Fill(2);
  if(nGoodTracksTPC == 2 && nGoodTracksITS == 2)fHistEvents->Fill(3);
  if(nGoodTracksTPC == 4 && nGoodTracksITS == 0)fHistEvents->Fill(4);
  if(nGoodTracksTPC == 3 && nGoodTracksITS == 1)fHistEvents->Fill(5);
  if(nGoodTracksTPC == 4 && nGoodTracksITS != 0)fHistEvents->Fill(6);
  if(nGoodTracksTPC == 2 && nGoodTracksITS != 0 && nGoodTracksITS != 2)fHistEvents->Fill(7);
  if(nGoodTracksTPC == 0 && nGoodTracksITS == 2)fHistEvents->Fill(8);
  if(nGoodTracksTPC == 0 && nGoodTracksITS == 4)fHistEvents->Fill(9);
     
  if(nGoodTracksTPC+nGoodTracksITS == 4 && nGoodTracksTPC > 1 && nGoodTracksSPD > 1 && (isMC || trigger.Contains("CCUP8-B"))){
    	MeanPt = GetMedian(TrackPtALL);
  	for(Int_t iTrack=0; iTrack<4; iTrack++) {
	AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(TrackIndexALL[iTrack]));
	if(trk->Pt() > 1.0) nHighPt++;
      		if(trk->Pt() > MeanPt){
			if(!trk->HasPointOnITSLayer(0) || !trk->HasPointOnITSLayer(1))continue;    
      			qLepton[nLepton] = trk->Charge();
			Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    			Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
			Float_t fPIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);

			nSigmaMuon[nLepton] = fPIDTPCMuon;
    			nSigmaElectron[nLepton] = fPIDTPCElectron;
			
			nSigmaProton[nLepton] = fPIDTOFProton;
			vProton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), protonMass);
			
			dEdx[nLepton] = trk->GetTPCsignal();
			
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
		vPsi2sCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fDiLeptonM = vDilepton.M();
  		fDiLeptonPt = vDilepton.Pt();
		if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0)) fSign = -1;
		else fSign = 1;
		Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
		Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
		Float_t nSigmaDistProton = TMath::Sqrt(TMath::Power(nSigmaProton[0],2) + TMath::Power(nSigmaProton[1],2));
		fChannel = 0;
		if((nSigmaDistMuon < nSigmaDistElectron)){
			fPIDsigma = nSigmaDistMuon; 
			if(nGoodTracksITS == 0)fChannel = 1;
			else fChannel = 10;
			}
		if((nSigmaDistMuon > nSigmaDistElectron)){
			fPIDsigma = nSigmaDistElectron;
			if(nGoodTracksITS == 0)fChannel = -1;
			else fChannel = -10;
			}
		if(nSigmaDistProton < 3){ 
			vJPsiCandidate = vProton[0]+vProton[1]+vPion[0]+vPion[1];
			fPIDsigma = nSigmaDistProton;
			fChannel = 3;
			FillTree(fTreeJPsi,vJPsiCandidate);
			}
		
		if(nHighPt > 0 )FillTree(fTreePsi2s,vPsi2sCandidate);
		}
  	}
	
  //Two track loop
  nHighPt = 0;
  fInEtaRec = kTRUE;
  if(nGoodTracksTPC == 2 && nGoodTracksLoose == 2 && nGoodTracksSPD == 2 && nGoodTracksITS == 0){
  	for(Int_t iTrack=0; iTrack<2; iTrack++) {
    	AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(TrackIndexTPC[iTrack]));
	fTrackIndex[iTrack] = TrackIndexTPC[iTrack]; 
	
	if(TMath::Abs(trk->Eta())>cutEta) fInEtaRec = kFALSE;
	
	if(trk->Pt() > 1.0) nHighPt++;
	
	Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    	Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
	Float_t fPIDTPCPion = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
	
	Float_t fPIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
	
    	qTrack[iTrack] = trk->Charge();
	
	vElectron[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
    	vMuon[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
    	nSigmaMuon[iTrack] = fPIDTPCMuon;
    	nSigmaElectron[iTrack] = fPIDTPCElectron;
    
    	vPion[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
    	nSigmaPion[iTrack] = fPIDTPCPion;
	
    	vProton[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), protonMass);
	nSigmaProton[iTrack] = fPIDTOFProton;
 	
	dEdx[iTrack] = trk->GetTPCsignal();
    	}
 
  fChannel = 0;
  if(qTrack[0]*qTrack[1]<0)fSign = -1;
  if(qTrack[0]*qTrack[1]>0)fSign = 1;
  vRhoCandidate = vPion[0]+vPion[1];

  Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
  Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
  Float_t nSigmaDistProton = TMath::Sqrt(TMath::Power(nSigmaProton[0],2) + TMath::Power(nSigmaProton[1],2));

  if(nSigmaDistProton < 5 && (isMC || trigger.Contains("CCUP8-B"))){ 
  	  fPIDsigma = nSigmaDistProton;
  	  vJPsiCandidate = vProton[0]+vProton[1];
  	  fChannel = 2;
  	  FillTree(fTreeJPsi,vJPsiCandidate);
  	  }
  if(nSigmaDistMuon < nSigmaDistElectron && (isMC || trigger.Contains("CCUP8-B"))){
  	  fPIDsigma = nSigmaDistMuon; 
  	  vJPsiCandidate = vMuon[0]+vMuon[1];
  	  fChannel = 1;
  	  FillTree(fTreeJPsi,vJPsiCandidate);
  	  }
  if(nSigmaDistMuon > nSigmaDistElectron && (isMC || trigger.Contains("CCUP8-B"))){ 
  	  fPIDsigma = nSigmaDistElectron;
  	  vJPsiCandidate = vElectron[0]+vElectron[1];
  	  fChannel = -1;
  	  FillTree(fTreeJPsi,vJPsiCandidate);
  	  }

  if(nSigmaDistMuon > nSigmaDistElectron &&(isMC || trigger.Contains("CCUP9-B"))){
  	  fPIDsigma = nSigmaDistElectron;
  	  vRhoCandidate = vElectron[0]+vElectron[1];
  	  fChannel = -1;
  	  FillTree(fTreeRho,vRhoCandidate);
	  }

  hTPCPIDPionCorr->Fill(nSigmaPion[0],nSigmaPion[1]);
  hTPCPIDPion->Fill(nSigmaPion[0]);hTPCPIDPion->Fill(nSigmaPion[1]);

  if(vJPsiCandidate.M()>2.8 && vJPsiCandidate.M()<3.2){ 
  	  hTPCdEdxCorr->Fill(dEdx[0],dEdx[1]);
  	  hTPCdEdxCorr->Fill(dEdx[1],dEdx[0]);
  	  hTPCPIDMuonCorr->Fill(nSigmaMuon[0],nSigmaMuon[1]);
  	  hTPCPIDMuonCorr->Fill(nSigmaMuon[1],nSigmaMuon[0]);
  	  hTPCPIDMuon->Fill(nSigmaMuon[0]);hTPCPIDMuon->Fill(nSigmaMuon[1]);
  	  hTPCPIDElectronCorr->Fill(nSigmaElectron[0],nSigmaElectron[1]);
  	  hTPCPIDElectronCorr->Fill(nSigmaElectron[1],nSigmaElectron[0]);
  	  hTPCPIDElectron->Fill(nSigmaElectron[0]);hTPCPIDElectron->Fill(nSigmaElectron[1]); 
  	  hTOFPIDProtonCorr->Fill(nSigmaProton[0],nSigmaProton[1]);
  	  hTOFPIDProtonCorr->Fill(nSigmaProton[1],nSigmaProton[0]);
  	  hTOFPIDProton->Fill(nSigmaProton[0]);hTOFPIDProton->Fill(nSigmaProton[1]);
  	  }
  } 
  
  
  //Two track loop
  if(nGoodTracksITS == 2 && nGoodTracksTPC== 0){
  	for(Int_t iTrack=0; iTrack<2; iTrack++) {
    	AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(TrackIndexITS[iTrack]));   
    
    	Float_t fPIDITSElectron = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kElectron);
    	Float_t fPIDITSPion = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kPion);
    	Float_t fPIDITSKaon = fPIDResponse->NumberOfSigmasITS(trk,AliPID::kKaon);
	        
    	qTrack[iTrack] = trk->Charge();
   
    	vKaon[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);    
    	vPion[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
	vElectron[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
	
	nSigmaKaon[iTrack] = fPIDITSKaon;
	nSigmaPion[iTrack] = fPIDITSPion;
	nSigmaElectron[iTrack] = fPIDITSElectron;
	 
    	}
	Float_t nSigmaDistKaon = TMath::Sqrt(TMath::Power(nSigmaKaon[0],2) + TMath::Power(nSigmaKaon[1],2));
	Float_t nSigmaDistPion = TMath::Sqrt(TMath::Power(nSigmaPion[0],2) + TMath::Power(nSigmaPion[1],2));
	Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
	
  	if(qTrack[0]*qTrack[1]<0)fSign = -1;
	if(qTrack[0]*qTrack[1]>0)fSign = 1;
	fChannel = 0;
	vRhoCandidate = vPion[0]+vPion[1];
	//FillTree(fTreeRho,vRhoCandidate);
	
	if(trigger.Contains("CCUP9-B")){ 
		if(nSigmaDistPion > 4 && nSigmaDistElectron > 4 && nSigmaDistKaon < 4){
			fChannel = 1;
			vPhiCandidate = vKaon[0]+vKaon[1];
			FillTree(fTreePhi,vPhiCandidate);
			}
		if(nSigmaDistPion > 4 && nSigmaDistElectron < 4 && nSigmaDistKaon > 4){
			fChannel = -1;
			vPhiCandidate = vElectron[0]+vElectron[1];
			FillTree(fTreePhi,vPhiCandidate);
			}
		}

	hITSPIDKaonCorr->Fill(nSigmaKaon[0],nSigmaKaon[1]);
	hITSPIDKaon->Fill(nSigmaKaon[0]);hITSPIDKaon->Fill(nSigmaKaon[1]);
  } 
  
  PostData(1, fOutputList);

}//UserExec


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::RunMC(AliVEvent *fEvent)
{
  
  for(Int_t i=0; i<10; i++) fTriggerInputsMC[i] = kFALSE;
  UShort_t fTriggerAD = fEvent->GetADData()->GetTriggerBits();
  UShort_t fTriggerVZERO = fEvent->GetVZEROData()->GetTriggerBits();
  UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();
  
  //fTriggerInputsMC[0] = fL0inputs & (1 << 9);   //0VBA VZERO A
  //fTriggerInputsMC[1] = fL0inputs & (1 << 10);   //0VBC VZERO C
  fTriggerInputsMC[0] = fTriggerVZERO & (1 << 12); //0VBA VZERO A
  fTriggerInputsMC[1] = fTriggerVZERO & (1 << 13); //0VBC VZERO C
  fTriggerInputsMC[2] = fTriggerAD & (1 << 12);   //0UBA ADA
  fTriggerInputsMC[3] = fTriggerAD & (1 << 13);   //0UBC ADC
  fTriggerInputsMC[4] = fL0inputs & (1 << 22);  //0OMU TOF two hits with topology
  fTriggerInputsMC[5] = fL0inputs & (1 << 19);	//0OM2 TOF two hits
  		
		
  const Int_t bcMod4 = TMath::Nint(hBCmod4->GetRandom());			
  //SPD inputs
  const AliVMultiplicity *mult = fEvent->GetMultiplicity();
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
  
  TLorentzVector vGenerated, vDecayProduct;
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  vGenerated.SetXYZM(0.,0.,0.,0.);
  Bool_t motherFound = kFALSE;
  fInEtaGen = kTRUE;
  
  AliMCEvent *mc = MCEvent();
  if(!mc) return;

  for(Int_t imc=0; imc<mc->GetNumberOfTracks(); imc++) {
    AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
    if(!mcPart) continue;
    
    if((TMath::Abs(mcPart->PdgCode()) == 443 || TMath::Abs(mcPart->PdgCode()) == 100443) && mcPart->GetMother() == -1){//Mother found
    	fPt = mcPart->Pt();
  	fY  = mcPart->Y();
  	fM  = mcPart->Particle()->GetCalcMass();
	fTreeGen->Fill();
	motherFound = kTRUE;
	break;
        }
    
    if(TMath::Abs(mcPart->PdgCode()) == 11 || 
       TMath::Abs(mcPart->PdgCode()) == 13 ||
       TMath::Abs(mcPart->PdgCode()) == 2212 ||
       TMath::Abs(mcPart->PdgCode()) == 211||
       TMath::Abs(mcPart->PdgCode()) == 22){
       
       if(mcPart->GetMother() != -1)continue;
       if(TMath::Abs(mcPart->PdgCode()) != 22 && TMath::Abs(mcPart->Eta())>cutEta) fInEtaGen = kFALSE;
    
       TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
       vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
    
       vGenerated += vDecayProduct;
       }
       
    
  }//loop over mc particles 
  if(!motherFound)FillTree(fTreeGen,vGenerated);

}//RunMC



//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::Terminate(Option_t *) 
{
  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::FillTree(TTree *t, TLorentzVector v) {

  fPt      = v.Pt();
  if(v.E() != v.Pz())fY = v.Rapidity();
  else fY = -999;
  fM       = v.M();

  t->Fill();

}

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcNano_MB::GetMedian(Double_t *daArray) {
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
