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
#include <iostream>
// my headers
#include "AliAnalysisTaskUpcRho0.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TRandom.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
// #include "AliAODEvent.h"
#include "AliMCEvent.h"
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
#include "AliMCParticle.h"
// #include "AliCentrality.h"
// #include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
// #include "AliTriggerAnalysis.h"
// #include "AliAODMCHeader.h"
#include "AliDataFile.h"

ClassImp(AliAnalysisTaskUpcRho0);

AliAnalysisTaskUpcRho0::AliAnalysisTaskUpcRho0()
  : AliAnalysisTaskSE(),
    fPIDResponse(0), isMC(0), isUsingEffi(0), fTriggerName(0),
  	fRhoTree(0), fMCTree(0),
	BunchCrossNumber_T(0), OrbitNumber_T(0), PeriodNumber_T(0),
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), 
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), UBAfired_T(0), UBCfired_T(0), 
  	VBAfired_T(0), VBCfired_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0), VtxContrib_T(0), SpdVtxContrib_T(0),
  	VtxChi2_T(0),VtxNDF_T(0),
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), GenPart_T(0),
  	RunNum_MC_T(0), Mass_MC_T(0), Pt_MC_T(0), Rapidity_MC_T(0), Phi_MC_T(0), 
	fListHist(0),fSPDfile(0), hBCmod4(0), hSPDeff(0), fEfficiencyFileName(0), 
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0), fFOcorr(0), fGoodTracks(0), fTrackChi2(0) 
{
//Dummy constructor
}

AliAnalysisTaskUpcRho0::AliAnalysisTaskUpcRho0(const char *name, Bool_t _isMC)
  : AliAnalysisTaskSE(name),
    fPIDResponse(0), isMC(0),isUsingEffi(0), fTriggerName(0),
  	fRhoTree(0), fMCTree(0),
  	BunchCrossNumber_T(0), OrbitNumber_T(0), PeriodNumber_T(0),
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), 
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), UBAfired_T(0), UBCfired_T(0), 
  	VBAfired_T(0), VBCfired_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0),VtxContrib_T(0), SpdVtxContrib_T(0),
  	VtxChi2_T(0),VtxNDF_T(0),
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), GenPart_T(0),
  	RunNum_MC_T(0), Mass_MC_T(0), Pt_MC_T(0), Rapidity_MC_T(0), Phi_MC_T(0), 
	fListHist(0),fSPDfile(0), hBCmod4(0), hSPDeff(0), fEfficiencyFileName(0), 
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0), fFOcorr(0), fGoodTracks(0), fTrackChi2(0)
{
  Init();
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
  if (_isMC) DefineOutput(3, TTree::Class());

}

AliAnalysisTaskUpcRho0::~AliAnalysisTaskUpcRho0() 
{
  // Destructor
  if(fRhoTree){
	delete fRhoTree;
	fRhoTree = 0x0;
  }
   if(fMCTree){
	delete fMCTree;
	fMCTree = 0x0;
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
	for(Int_t i=0; i<9; i++){
    	fHistdEdxVsP[i]=0x0;
  	}
}

void AliAnalysisTaskUpcRho0::UserCreateOutputObjects() 
{
  	//PID response
 	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  	AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  	fPIDResponse = inputHandler->GetPIDResponse();

  	GenPart_T = new TClonesArray("TParticle", 1000);

	fRhoTree = new TTree("Selected","Selected Rho0 events");
	//define branches
	fRhoTree->Branch("RunNum_T",&RunNum_T,"RunNum_T/I");
	fRhoTree->Branch("PeriodNumber_T",&PeriodNumber_T,"PeriodNumber_T/i");
	fRhoTree->Branch("OrbitNumber_T",&OrbitNumber_T,"OrbitNumber_T/i");
	fRhoTree->Branch("BunchCrossNumber_T",&BunchCrossNumber_T,"BunchCrossNumber_T/s");
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
	fRhoTree->Branch("TPCsignal_T",&TPCsignal_T,"TPCsignal_T[2]/I");
	fRhoTree->Branch("TrackP_T",&TrackP_T,"TrackP_T[2]/F");
	fRhoTree->Branch("TrackEta_T",&TrackEta_T,"TrackEta_T[2]/F");
	fRhoTree->Branch("TrackPhi_T",&TrackPhi_T,"TrackPhi_T[2]/F");
	fRhoTree->Branch("TrackPx_T",&TrackPx_T,"TrackPx_T[2]/F");
	fRhoTree->Branch("TrackPy_T",&TrackPy_T,"TrackPy_T[2]/F");
	fRhoTree->Branch("TrackPz_T",&TrackPz_T,"TrackPz_T[2]/F");
	fRhoTree->Branch("VtxX_T",&Vertex_T[0],"VtxX_T/F");
	fRhoTree->Branch("VtxY_T",&Vertex_T[1],"VtxY_T/F");
	fRhoTree->Branch("VtxZ_T",&Vertex_T[2],"VtxZ_T/F");
	fRhoTree->Branch("VtxContrib_T",&VtxContrib_T,"VtxContrib_T/I");
	fRhoTree->Branch("VtxChi2_T",&VtxChi2_T,"VtxChi2_T/F");
	fRhoTree->Branch("VtxNDF_T",&VtxNDF_T,"VtxNDF_T/F");
	fRhoTree->Branch("SpdVtxX_T",&SpdVertex_T[0],"SpdVtxX_T/F");
	fRhoTree->Branch("SpdVtxY_T",&SpdVertex_T[1],"SpdVtxY_T/F");
	fRhoTree->Branch("SpdVtxZ_T",&SpdVertex_T[2],"SpdVtxZ_T/F");
	fRhoTree->Branch("SpdVtxContrib_T",&SpdVtxContrib_T,"SpdVtxContrib_T/I");
	fRhoTree->Branch("V0Adecision_T",&V0Adecision_T,"V0Adecision_T/I");
	fRhoTree->Branch("V0Cdecision_T",&V0Cdecision_T,"V0Cdecision_T/I");
	fRhoTree->Branch("ADAdecision_T",&ADAdecision_T,"ADAdecision_T/I");
	fRhoTree->Branch("ADCdecision_T",&ADCdecision_T,"ADCdecision_T/I");
	fRhoTree->Branch("UBAfired_T",&UBAfired_T,"UBAfired_T/O");
	fRhoTree->Branch("UBCfired_T",&UBCfired_T,"UBCfired_T/O");
	fRhoTree->Branch("VBAfired_T",&VBAfired_T,"VBAfired_T/O");
	fRhoTree->Branch("VBCfired_T",&VBCfired_T,"VBCfired_T/O");
	fRhoTree->Branch("Ntracklets_T",&Ntracklets_T,"Ntracklets_T/I");
	// fRhoTree->Branch("ITSModule_T",&ITSModule_T,"ITSModule_T/I");
	fRhoTree->Branch("ChipCut_T",&ChipCut_T,"ChipCut_T/O");

	// MC tree
	if (isMC){
	fMCTree = new TTree("Generated","Generated Rho0 events");
	//define branches
	fMCTree->Branch("RunNum_MC_T",&RunNum_MC_T,"RunNum_MC_T/I");
	fMCTree->Branch("Mass_MC_T",&Mass_MC_T,"Mass_MC_T/F");
	fMCTree->Branch("Pt_MC_T",&Pt_MC_T,"Pt_MC_T/F");
	fMCTree->Branch("Rapidity_MC_T",&Rapidity_MC_T,"Rapidity_MC_T/F");
	fMCTree->Branch("Phi_MC_T",&Phi_MC_T,"Phi_MC_T/F");
	}

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
  	fFOcorr = new TH2F("fFOcorr","fFOcorr",240,0,240,240,0,240);
  	fListHist->Add(fFOcorr);
  	fGoodTracks = new TH1F("fGoodTracks","fGoodTracks",5,0,5);
  	fListHist->Add(fGoodTracks);
  	fTrackChi2 = new TH1F("fTrackChi2","fTrackChi2",100,0,10);
  	fListHist->Add(fTrackChi2);

	// TPC clusters
	TPCclustersP = new TH1F("TPCclustersP","TPCclustersP",181,0,180); fListHist->Add(TPCclustersP);
	TPCclustersN = new TH1F("TPCclustersN","TPCclustersN",181,0,180); fListHist->Add(TPCclustersN);
	// dE/dx
	dEdx = new TH2F("dEdx","dEdx",500,0.1,1.5,100,0,200); fListHist->Add(dEdx);
	// eta-phi
	EtaPhiP = new TH2F("EtaPhiP","EtaPhiP",100,-1,1,100,0,2*3.14159); fListHist->Add(EtaPhiP);
	EtaPhiN = new TH2F("EtaPhiN","EtaPhiN",100,-1,1,100,0,2*3.14159); fListHist->Add(EtaPhiN);

	TString pNames[9]={"Elec","Muon","Pion","Kaon","Proton","Deuteron","Triton","He3","Alpha"};
	for(Int_t jsp=0; jsp<9; jsp++){
    	fHistdEdxVsP[jsp] = new TH2F(Form("hdEdxVsP%s",pNames[jsp].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    	fListHist->Add(fHistdEdxVsP[jsp]);
 	}

	// load SPD effi
	if (isUsingEffi) {
		std::cout<<"Using efficiency file: "<<fEfficiencyFileName<<std::endl;
		fSPDfile = AliDataFile::OpenOADB(fEfficiencyFileName.Data());
		fSPDfile->Print();
		fSPDfile->Map();
		hSPDeff = (TH2D*) fSPDfile->Get("hEff");
		hSPDeff->SetDirectory(0);
		TH2D *hBCmod4_2D = (TH2D*) fSPDfile->Get("hCounts");
		hBCmod4_2D->SetDirectory(0);
		hBCmod4 = hBCmod4_2D->ProjectionY();
		fSPDfile->Close();
	}

	PostData(1, fRhoTree);
	PostData(2, fListHist);
	if (isMC) PostData(3, fMCTree);
}

void AliAnalysisTaskUpcRho0::UserExec(Option_t *) 
{

  //input event
  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  // MC generated particles
  if(isMC){
	GenPart_T->Clear("C");

	AliMCEvent *mc = MCEvent();
	if(!mc) return;

	Int_t nmc = 0;
	// loop over mc particles
	for(Int_t imc=0; imc<mc->GetNumberOfTracks(); imc++) {
		AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
		if(!mcPart) continue;

		if(mcPart->GetMother() >= 0) continue;

		TParticle *part = (TParticle*) GenPart_T->ConstructedAt(nmc++);
		part->SetMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->E());
		part->SetPdgCode(mcPart->PdgCode());
		part->SetUniqueID(imc);
	  } // loop over mc particles
	  if(nmc == 2){
	  	TParticle *mcp0,*mcp1;
		TLorentzVector lv0,lv1,lvSum;
		//load particle
		mcp0 = (TParticle*) GenPart_T->At(0);
		mcp1 = (TParticle*) GenPart_T->At(1);
		//create fourvector
		lv0.SetPxPyPzE(mcp0->Px(), mcp0->Py(), mcp0->Pz(), mcp0->Energy());
		lv1.SetPxPyPzE(mcp1->Px(), mcp1->Py(), mcp1->Pz(), mcp1->Energy());
		lvSum = lv0+lv1;
		//connect variables
		RunNum_MC_T = esd->GetRunNumber();
		Mass_MC_T = lvSum.M();
		Pt_MC_T = lvSum.Pt();
		Rapidity_MC_T = lvSum.Rapidity();
		Phi_MC_T = lvSum.Phi();
	  }
	fMCTree->Fill();
  } // end of MC generated particles

  // trigger
  TString trigger = esd->GetFiredTriggerClasses();

  // triggered in data for lumi scalling
  if(!isMC && trigger.Contains(fTriggerName.Data())) {
  	fHistTriggersPerRun->Fill(esd->GetRunNumber());
  	// cout<<trigger<<endl;
  }

  // CCUP9-B - *0VBA *0VBC *0UBA *0UBC 0STP
  if (!isMC) { // data
  	if (!trigger.Contains(fTriggerName.Data())) return;
  }
  else { // MC
  	if (!IsTriggered(esd)) return;
  } // end of MC trigger

  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();

  Int_t nGoodTracks=0;
  Int_t TrackIndex[2] = {-1,-1};
// cout<<"starting track loop"<<endl;
  //Track loop - cuts
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    if( !trk ) continue;
 	if( trk->IsOn(AliESDtrack::kITSpureSA) ) continue;
    if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
    if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
    if(trk->GetTPCNcls() < 50)continue;
    // if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
    if(!((trk->HasPointOnITSLayer(0))&&(trk->HasPointOnITSLayer(1)))) continue;
    Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
    if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
    trk->GetImpactParameters(dca[0],dca[1]);
    if(TMath::Abs(dca[1]) > 2) continue;
    Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
    if(TMath::Abs(dca[0]) > cut_DCAxy) continue;

	// store good track index
	TrackIndex[nGoodTracks] = itr;

    nGoodTracks++;
    // cout<<nGoodTracks<<" good tracks"<<endl;
	if(nGoodTracks > 5) break; // just to know how many nGoodTrack are there

  }//Track loop end

  	fGoodTracks->Fill(nGoodTracks);

  if(nGoodTracks == 2){ // fill tree variables
 // cout<<"two good tracks"<<endl;
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
  			fFOmodules[(chipkey/5)]++;
 		}
 	}

	// event info
	RunNum_T = esd->GetRunNumber();
	OrbitNumber_T = esd->GetOrbitNumber();
	PeriodNumber_T = esd->GetPeriodNumber();
	BunchCrossNumber_T = esd->GetBunchCrossNumber();

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

	UBAfired_T = esd->GetHeader()->IsTriggerInputFired("0UBA");
	UBCfired_T = esd->GetHeader()->IsTriggerInputFired("0UBC");
	VBAfired_T = esd->GetHeader()->IsTriggerInputFired("0VBA");
	VBCfired_T = esd->GetHeader()->IsTriggerInputFired("0VBC");

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
	VtxContrib_T = fESDVertex->GetNContributors();
	Vertex_T[0] = fESDVertex->GetX();
	Vertex_T[1] = fESDVertex->GetY();
	Vertex_T[2] = fESDVertex->GetZ();
	VtxChi2_T = fESDVertex->GetChi2();
	VtxNDF_T = fESDVertex->GetNDF();

	//SPD primary vertex
	AliESDVertex *fSPDVertex = (AliESDVertex*) esd->GetPrimaryVertexSPD();
	SpdVtxContrib_T = fSPDVertex->GetNContributors();
	SpdVertex_T[0] = fSPDVertex->GetX();
	SpdVertex_T[1] = fSPDVertex->GetY();
	SpdVertex_T[2] = fSPDVertex->GetZ();

	// Tracklets
	Ntracklets_T = esd->GetMultiplicity()->GetNumberOfTracklets();

	// loop over two good tracks
  	for(Int_t i=0; i<2; i++){
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);

		ITSModuleInner_T[i] = trk->GetITSModuleIndex(0)/1000000;
		ITSModuleOuter_T[i] = trk->GetITSModuleIndex(1)/1000000;

		// TPC PID n-sigma
		PIDTPCElectron_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		PIDTPCPion_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);

		// separated PID
		Double_t ptrackTPC=-999.;
		const AliExternalTrackParam* ippar=trk->GetInnerParam();
		if(ippar) ptrackTPC=ippar->P();
		Double_t dedx=trk->GetTPCsignal();
		Int_t  pidtr=trk->GetPIDForTracking();
		if(pidtr>=0 && pidtr<9) fHistdEdxVsP[pidtr]->Fill(ptrackTPC,dedx);

		charge[i] = trk->Charge();
		TPCsignal_T[i] = trk->GetTPCsignal();
		TrackP_T[i] = trk->P();
		TrackPhi_T[i] = trk->Phi();
		TrackEta_T[i] = trk->Eta();
		TrackPx_T[i] = trk->Px();
		TrackPy_T[i] = trk->Py();
		TrackPz_T[i] = trk->Pz();

		fTrackChi2->Fill((Float_t)trk->GetTPCchi2()/trk->GetTPCNcls());

		fITSmodule->Fill(ITSModuleInner_T[i]);
		fITSmodule->Fill(ITSModuleOuter_T[i]);

		for(Int_t i=0;i<240;i++){
			if (fFOmodules[i] > 0){
				fFOcorr->Fill(ITSModuleInner_T[i],i);
				fFOcorr->Fill(ITSModuleOuter_T[i],i);
			}
		}

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
  	} // end loop over two good tracks

  	lvSum = lv[0]+lv[1];

	if (charge[0]*charge[1]>0) LikeSign_T = 1;
	else LikeSign_T = 0;
	Mass_T = lvSum.M();
	Pt_T = lvSum.Pt();
	Rapidity_T = lvSum.Rapidity();
	Phi_T = lvSum.Phi();

	// virtual cut on FO chip matching
	Int_t SPDInner[20]; for (Int_t i=0; i<20; ++i) SPDInner[i]=0;
	Int_t SPDOuter[40]; for (Int_t i=0; i<40; ++i) SPDOuter[i]=0;

	SPDInner[ITSModuleInner_T[0]/4]++;
	SPDInner[ITSModuleInner_T[1]/4]++;
	SPDOuter[(ITSModuleOuter_T[0]-80)/4]++;
	SPDOuter[(ITSModuleOuter_T[1]-80)/4]++;
  
	ChipCut_T = 0;
	if ((fTriggerName == "CCUP9-B") &&
		((fFOmodules[ITSModuleInner_T[0]] == 0)||(fFOmodules[ITSModuleOuter_T[0]] == 0)
		||(fFOmodules[ITSModuleInner_T[1]] == 0)||(fFOmodules[ITSModuleOuter_T[1]] == 0)
		|| !Is0STPfired(SPDInner,SPDOuter))) ChipCut_T = 1;

	if ((fTriggerName == "CCUP2-B") &&
		((fFOmodules[ITSModuleInner_T[0]] == 0)||(fFOmodules[ITSModuleOuter_T[0]] == 0)
		||(fFOmodules[ITSModuleInner_T[1]] == 0)||(fFOmodules[ITSModuleOuter_T[1]] == 0)
		)) ChipCut_T = 1;
	if ((fTriggerName == "CCUP4-B") &&
		((fFOmodules[ITSModuleInner_T[0]] == 0)||(fFOmodules[ITSModuleOuter_T[0]] == 0)
		||(fFOmodules[ITSModuleInner_T[1]] == 0)||(fFOmodules[ITSModuleOuter_T[1]] == 0)
		)) ChipCut_T = 1;
	if ((fTriggerName == "C1ZED") &&
		((fFOmodules[ITSModuleInner_T[0]] == 0)||(fFOmodules[ITSModuleOuter_T[0]] == 0)
		||(fFOmodules[ITSModuleInner_T[1]] == 0)||(fFOmodules[ITSModuleOuter_T[1]] == 0)
		)) ChipCut_T = 1;

    Int_t fFOcounter = 0;
  	for(Int_t chipkey=0;chipkey<1200;chipkey++){
  		if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)){
  			fFOchip->Fill(chipkey);
  			fFOcounter++;
 		}
  	}
  	fFOcount->Fill(fFOcounter);
  	// fFOcorr->Fill();

  //fill
  fRhoTree ->Fill();

  } // end 2 good tracks
// cout<<"saving data"<<endl;
  PostData(1, fRhoTree);
  PostData(2, fListHist);
  if (isMC) PostData(3, fMCTree);
  

}//UserExec

// fuction that get two arrays and return if 0STP trigger was fired
Bool_t AliAnalysisTaskUpcRho0::Is0STPfired(Int_t *vPhiInner, Int_t *vPhiOuter) // array 20, 40
{
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
	if (fired != 0) return kTRUE;
	else return kFALSE;
}

Bool_t AliAnalysisTaskUpcRho0::IsTriggered(AliESDEvent *esd)
// return kTRUE if CCUP9 triggered was fired
{
	Bool_t V0A = kFALSE;
	Bool_t V0C = kFALSE;
	Bool_t ADA = kFALSE;
	Bool_t ADC = kFALSE;
	Bool_t STP = kFALSE;
	Bool_t SMB = kFALSE;
	Bool_t SM2 = kFALSE;
	Bool_t SH1 = kFALSE;
	Bool_t OM2 = kFALSE;
	Bool_t OMU = kFALSE;
	//SPD inputs
	Int_t bcMod4 = 0;
	if (isUsingEffi) bcMod4 = TMath::Nint(hBCmod4->GetRandom());
	AliMultiplicity *mult = esd->GetMultiplicity();
	Int_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=0;
	Int_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=0;

	Int_t nInner(0), nOuter(0);
	for (Int_t i(0); i<1200; ++i) {
		Double_t eff = 1;
		if (isUsingEffi) eff = hSPDeff->GetBinContent(1+i, 1+bcMod4);
		Bool_t isFired = (mult->TestFastOrFiredChips(i)) && (gRandom->Uniform(0,1) < eff);
		if (i<400) {
			vPhiInner[i/20] += isFired;
			nInner += isFired;
		} else {
			vPhiOuter[(i-400)/20] += isFired;
			nOuter += isFired;
		}
		}
	// 0STP
	STP = Is0STPfired(vPhiInner,vPhiOuter);
	// 0SMB - At least one hit in SPD
	if (nOuter > 0 || nInner > 0) SMB = kTRUE;
	// 0SM2 - Two hits on outer layer
	if (nOuter > 1) SM2 = kTRUE;
	// 0SH1 - More then 6 hits on outer layer
	// if (nOuter >= 7) SH1 = kTRUE;
	//0SH1 2017 - Two hits on inner and outer layer
	if (nInner >= 2 && nOuter >= 2) SH1 = kTRUE;
	// V0
	V0A = esd->GetHeader()->IsTriggerInputFired("0VBA");
	V0C = esd->GetHeader()->IsTriggerInputFired("0VBC");
	// AD
	ADA = esd->GetHeader()->IsTriggerInputFired("0UBA");
	ADC = esd->GetHeader()->IsTriggerInputFired("0UBC");
	// TOF
	OM2 = esd->GetHeader()->IsTriggerInputFired("0OM2");
	OMU = esd->GetHeader()->IsTriggerInputFired("0OMU");
	  
	if ((fTriggerName == "CCUP9-B") && (!V0A && !V0C && !ADA && !ADC && STP)) return kTRUE; // CCUP9 is fired
	if ((fTriggerName == "CCUP2-B") && (!V0A && !V0C && SM2 && OM2)) return kTRUE; // CCUP2 is fired works only in 2015
	if ((fTriggerName == "CCUP4-B") && (!V0A && !V0C && SM2 && OMU)) return kTRUE; // CCUP4 is fired works only in 2015

	else return kFALSE;
} // end of MC trigger