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

// my headers
#include "AliAnalysisTaskUpcRho0.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
// #include "AliAODEvent.h"
// #include "AliMCEvent.h"
// #include "AliAODVZERO.h"
// #include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
// #include "AliAODTrack.h"
// #include "AliAODPid.h"
// #include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
// #include "AliESDMuonTrack.h"
// #include "AliAODMCParticle.h"
// #include "AliMCParticle.h"
// #include "AliCentrality.h"
// #include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
// #include "AliTriggerAnalysis.h"
// #include "AliAODMCHeader.h"

ClassImp(AliAnalysisTaskUpcRho0);

AliAnalysisTaskUpcRho0::AliAnalysisTaskUpcRho0()
  : AliAnalysisTaskSE(),
    fPIDResponse(0),
  	fRhoTree(0), 
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), 
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0), DeltaPhi_T(0), 
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), ITSModule_T(0),
	fListHist(0),
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0)
{
//Dummy constructor
}

AliAnalysisTaskUpcRho0::AliAnalysisTaskUpcRho0(const char *name)
  : AliAnalysisTaskSE(name),
    fPIDResponse(0),
  	fRhoTree(0), 
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), 
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0), DeltaPhi_T(0), 
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), ITSModule_T(0),
	fListHist(0),
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0)
{
  Init();
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
}

AliAnalysisTaskUpcRho0::~AliAnalysisTaskUpcRho0() 
{
  // Destructor
  if(fRhoTree){
	delete fRhoTree;
	fRhoTree = 0x0;
  }
  if(fListHist){
	delete fListHist;
	fListHist = 0x0;
  }
}

void AliAnalysisTaskUpcRho0::Init()
{
	for (Int_t i=0;i<2;i++){
		PIDTPCPion_T[i] = -666;
		PIDTPCElectron_T[i] = -666;
		TPCsignal_T[i] = -666;
		TrackP_T[i] = -666;
		TrackEta_T[i] = -666;
		TrackPhi_T[i] = -666;
	}
	for (Int_t i=0;i<3;i++){
		Vertex_T[i] = -666;
		SpdVertex_T[i] = -666;
	}
	for (Int_t i=0;i<4;i++){
		ZDCAtime_T[i] = -666;
		ZDCCtime_T[i] = -666;
	}
}

void AliAnalysisTaskUpcRho0::UserCreateOutputObjects() 
{
  	//PID response
 	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  	AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  	fPIDResponse = inputHandler->GetPIDResponse();

	fRhoTree = new TTree("Selected","Selected Rho0 events");
	//define branches
	fRhoTree->Branch("RunNum_T",&RunNum_T,"RunNum_T/I");
	fRhoTree->Branch("LikeSign_T",&LikeSign_T,"LikeSign_T/O");
	fRhoTree->Branch("Mass_T",&Mass_T,"Mass_T/F");
	fRhoTree->Branch("Pt_T",&Pt_T,"Pt_T/F");
	fRhoTree->Branch("Rapidity_T",&Rapidity_T,"Rapidity_T/F");
	fRhoTree->Branch("Phi_T",&Phi_T,"Phi_T/F");
	fRhoTree->Branch("ZNAenergy_T",&ZNAenergy_T,"ZNAenergy_T/F");
	fRhoTree->Branch("ZNCenergy_T",&ZNCenergy_T,"ZNCenergy_T/F");
	fRhoTree->Branch("ZPAenergy_T",&ZPAenergy_T,"ZPAenergy_T/F");
	fRhoTree->Branch("ZPCenergy_T",&ZPCenergy_T,"ZPCenergy_T/F");
	fRhoTree->Branch("ZDCAtime_T",&ZDCAtime_T,"ZDCAtime_T[4]/F");
	fRhoTree->Branch("ZDCCtime_T",&ZDCCtime_T,"ZDCCtime_T[4]/F");
	fRhoTree->Branch("PIDTPCPion_T",&PIDTPCPion_T,"PIDTPCPion_T[2]/F");
	fRhoTree->Branch("PIDTPCElectron_T",&PIDTPCElectron_T,"PIDTPCElectron_T[2]/F");
	// fRhoTree->Branch("TPCsignal_T",&TPCsignal_T,"TPCsignal_T[2]/I");
	// fRhoTree->Branch("TrackP_T",&TrackP_T,"TrackP_T[2]/F");
	// fRhoTree->Branch("TrackEta_T",&TrackEta_T,"TrackEta_T[2]/F");
	// fRhoTree->Branch("TrackPhi_T",&TrackPhi_T,"TrackPhi_T[2]/F");
	fRhoTree->Branch("VtxX_T",&Vertex_T[0],"VtxX_T/F");
	fRhoTree->Branch("VtxY_T",&Vertex_T[1],"VtxY_T/F");
	fRhoTree->Branch("VtxZ_T",&Vertex_T[2],"VtxZ_T/F");
	// fRhoTree->Branch("SpdVtxX_T",&SpdVertex_T[0],"SpdVtxX_T/F");
	// fRhoTree->Branch("SpdVtxY_T",&SpdVertex_T[1],"SpdVtxY_T/F");
	// fRhoTree->Branch("SpdVtxZ_T",&SpdVertex_T[2],"SpdVtxZ_T/F");
	fRhoTree->Branch("V0Adecision_T",&V0Adecision_T,"V0Adecision_T/I");
	fRhoTree->Branch("V0Cdecision_T",&V0Cdecision_T,"V0Cdecision_T/I");
	fRhoTree->Branch("ADAdecision_T",&ADAdecision_T,"ADAdecision_T/I");
	fRhoTree->Branch("ADCdecision_T",&ADCdecision_T,"ADCdecision_T/I");
	// fRhoTree->Branch("DeltaPhi_T",&DeltaPhi_T,"DeltaPhi_T/F");
	fRhoTree->Branch("Ntracklets_T",&Ntracklets_T,"Ntracklets_T/I");
	//fRhoTree->Branch("SpdVtxContrib_T",&fSpdVtxContrib,"SpdVtxContrib_T/I");
	fRhoTree->Branch("ITSModule_T",&ITSModule_T,"ITSModule_T/I");
	fRhoTree->Branch("ChipCut_T",&ChipCut_T,"ChipCut_T/B");

	fListHist = new TList();
  	fListHist ->SetOwner();
  
  	fHistTriggersPerRun = new TH1I("fHistTriggersPerRun", "fHistTriggersPerRun", 50000, 240000.5, 290000.5);
  	fListHist->Add(fHistTriggersPerRun);
  	fITSmodule = new TH1I("fITSmodule","fITSmodule",240,0,240);
  	fListHist->Add(fITSmodule);
  	fFOchip = new TH1I("fFOchip","fFOchip",1200,0,1200);
  	fListHist->Add(fFOchip);
  	fFOcount = new TH1I("fFOcount","fFOcount",30,0,30);
  	fListHist->Add(fFOcount);

	// TPC clusters
	TPCclustersP = new TH1F("TPCclustersP","TPCclustersP",181,0,180); fListHist->Add(TPCclustersP);
	TPCclustersN = new TH1F("TPCclustersN","TPCclustersN",181,0,180); fListHist->Add(TPCclustersN);
	// dE/dx
	dEdx = new TH2F("dEdx","dEdx",500,0.1,1.5,100,0,200); fListHist->Add(dEdx);
	// eta-phi
	EtaPhiP = new TH2F("EtaPhiP","EtaPhiP",100,-1,1,100,0,2*3.14159); fListHist->Add(EtaPhiP);
	EtaPhiN = new TH2F("EtaPhiN","EtaPhiN",100,-1,1,100,0,2*3.14159); fListHist->Add(EtaPhiN);

	PostData(1, fRhoTree);
	PostData(2, fListHist);
}

void AliAnalysisTaskUpcRho0::UserExec(Option_t *) 
{

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  // data
  RunNum_T = esd->GetRunNumber();

  // trigger
  TString trigger = esd->GetFiredTriggerClasses();

  if(trigger.Contains("CCUP9-B")) fHistTriggersPerRun->Fill(RunNum_T); //CCUP9 triggers

  // CCUP9-B - *0VBA *0VBC *0UBA *0UBC 0STP
  if (!trigger.Contains("CCUP9-B")) return; 

  // VZERO, ZDC, AD
  AliESDVZERO *fV0data = esd->GetVZEROData();
  AliESDZDC *fZDCdata = esd->GetESDZDC();
  AliESDAD *fADdata = esd->GetADData();
  
  V0Adecision_T = fV0data->GetV0ADecision();
  V0Cdecision_T = fV0data->GetV0CDecision();
  if(fADdata){
  	ADAdecision_T = fADdata->GetADADecision();
  	ADCdecision_T = fADdata->GetADCDecision();
	}

  // ZN energy
  ZNAenergy_T = fZDCdata->GetZNATowerEnergy()[0];
  ZNCenergy_T = fZDCdata->GetZNCTowerEnergy()[0];
  ZPAenergy_T = fZDCdata->GetZPATowerEnergy()[0];
  ZPCenergy_T = fZDCdata->GetZPCTowerEnergy()[0];

  // neutron ZDC time
  Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
  Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
  if (esd->GetRunNumber()>=245726 && esd->GetRunNumber()<=245793) detChZNA = 10;
  for (Int_t i=0;i<4;i++){ 
	ZDCAtime_T[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
	ZDCCtime_T[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
  }

  // primary vertex
  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
  // VtxContrib = fESDVertex->GetNContributors();
  Vertex_T[0] = fESDVertex->GetX();
  Vertex_T[1] = fESDVertex->GetY();
  Vertex_T[2] = fESDVertex->GetZ();

  // Tracklets
  Ntracklets_T = esd->GetMultiplicity()->GetNumberOfTracklets();

  Int_t nGoodTracks=0;
  Int_t TrackIndex[2] = {-1,-1};

  //Track loop - cuts
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
    Bool_t isMC = kFALSE;
    if(!isMC){
      if(TMath::Abs(dca[1]) > 2) continue;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
	}

	// store good track index
	TrackIndex[nGoodTracks] = itr;

    nGoodTracks++;
	if(nGoodTracks > 2) break;

  }//Track loop end

  if(nGoodTracks == 2){

  	TDatabasePDG *pdgdat = TDatabasePDG::Instance(); 
  	TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  	Double_t pionMass = partPion->Mass();

  	Double_t charge[2]; // charge
	TLorentzVector lv[2];
	TLorentzVector lvSum; // pair-4vector

	Int_t fFOmodules[240];
	for (Int_t i=0;i<240;i++) fFOmodules[i] = 0;

  	for(Int_t chipkey=0;chipkey<1200;chipkey++){
  		if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)){
  			fFOmodules[chipkey/5]++;
 		}
 	}

	// loop over two good tracks
  	for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);

	  	// chips cut
	  	//if (fFOmodules[trk->GetITSModuleIndex(0)/1000000] == 0) return;
	  	//if (fFOmodules[trk->GetITSModuleIndex(1)/1000000] == 0) return;
		
		// contributor to Vertex
		// if(fESDVertex->UsesTrack(TrackIndex[i]))fIsVtxContributor[i] = kTRUE;
		// else fIsVtxContributor[i] = kFALSE;
		
		// TPC PID n-sigma
		PIDTPCElectron_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		PIDTPCPion_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);

		charge[i] = trk->Charge();
		TPCsignal_T[i] = trk->GetTPCsignal();
		TrackP_T[i] = trk->P();
		TrackPhi_T[i] = trk->Phi();
		TrackEta_T[i] = trk->Eta();
		ITSModule_T = trk->GetITSModuleIndex(0);

		fITSmodule->Fill(trk->GetITSModuleIndex(0)/1000000);
		fITSmodule->Fill(trk->GetITSModuleIndex(1)/1000000);
		dEdx->Fill(trk->Pt(),trk->GetTPCsignal());
		if (trk->Charge()>0) {
			TPCclustersP->Fill(trk->GetTPCNcls());
			EtaPhiP->Fill(trk->Eta(),trk->Phi());
		}
		else {
			TPCclustersN->Fill(trk->GetTPCNcls());
			EtaPhiN->Fill(trk->Eta(),trk->Phi());
		}

		lv[i].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
  	}

  	lvSum = lv[0]+lv[1];

  	//RunNum_T = fRunNum;
	if (charge[0]*charge[1]>0) LikeSign_T = 1;
	else LikeSign_T = 0;
	Mass_T = lvSum.M();
	Pt_T = lvSum.Pt();
	Rapidity_T = lvSum.Rapidity();
	Phi_T = lvSum.Phi();

	// virtual cut on FO chip matching
	ChipCut_T = 0;
	if ((fFOmodules[esd->GetTrack(TrackIndex[0])->GetITSModuleIndex(0)/1000000] == 0)
		||(fFOmodules[esd->GetTrack(TrackIndex[0])->GetITSModuleIndex(1)/1000000] == 0)
		||(fFOmodules[esd->GetTrack(TrackIndex[1])->GetITSModuleIndex(0)/1000000] == 0)
		||(fFOmodules[esd->GetTrack(TrackIndex[1])->GetITSModuleIndex(1)/1000000] == 0)
		) ChipCut_T = 1;
  
  	//fill
  	fRhoTree ->Fill();

    Int_t fFOcounter = 0;
  	for(Int_t chipkey=0;chipkey<1200;chipkey++){
  		if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)){
  			fFOchip->Fill(chipkey);
  			fFOcounter++;
 		}
  	}
  	fFOcount->Fill(fFOcounter);

  }

  PostData(1, fRhoTree);
  PostData(2, fListHist);

}//UserExec