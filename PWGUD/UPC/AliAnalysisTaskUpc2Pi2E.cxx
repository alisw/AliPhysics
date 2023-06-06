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
#include <string.h>
#include <bitset>
// my headers
#include "AliAnalysisTaskUpc2Pi2E.h"

#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include <TObjString.h>
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TRandom.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
//#include "AliAODEvent.h"
#include "AliMCEvent.h"
//#include "AliAODVZERO.h"
//#include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
//#include "AliAODTrack.h"
//#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDMuonTrack.h"
//#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
// #include "AliCentrality.h"
// #include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
// #include "AliTriggerAnalysis.h"
// #include "AliAODMCHeader.h"
#include "AliDataFile.h"
#include "AliOADBContainer.h"

ClassImp(AliAnalysisTaskUpc2Pi2E);

AliAnalysisTaskUpc2Pi2E::AliAnalysisTaskUpc2Pi2E()
  : AliAnalysisTaskSE(),
    fPIDResponse(0), fTrackCutsBit4(0), isMC(0), debugMode(0), isUsingEffi(0), fTriggerName(0),
  	f2Pi2ETree(0), f2Pi2ETree1(0), fMCTree(0), fTPCNcls(50), fOption(0),
	BunchCrossNumber_T(0), OrbitNumber_T(0), PeriodNumber_T(0),
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), NTracks_T(0), 
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), UBAfired_T(0), UBCfired_T(0), 
  	VBAfired_T(0), VBCfired_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0), VtxContrib_T(0), SpdVtxContrib_T(0),
  	VtxChi2_T(0),VtxNDF_T(0), TriggerTOF_T(0), TriggerSPD_T(0),
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), GenPart_T(0),ntrk(0),
  	RunNum_MC_T(0), Mass_MC_T(0), Pt_MC_T(0), Rapidity_MC_T(0), Phi_MC_T(0), 
	fListHist(0),fSPDfile(0), hBCmod4(0), hSPDeff(0), fEfficiencyFileName(0), 
	fTOFfile(0),fTOFFileName(0), hTOFeff(0), fTOFmask(0), isUsingTOFeff(0),
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0), fFOcorr(0), fGoodTracks(0), fTrackChi2(0),
	fHistdEdxVsP1(0),fHistdEdxVsP2(0),fHistdEdxVsP3(0),fHistdEdxVsP4(0),fHistdEdxVsP5(0),
	fHistdEdxVsP6(0),fHistdEdxVsP7(0),fHistdEdxVsP8(0),fHistdEdxVsP9(0) , fDeltaPhiRho(0), 
	fDeltaPhiEe(0), FORChip(0),Trigger_T(0)
{
//Dummy constructor
}

AliAnalysisTaskUpc2Pi2E::AliAnalysisTaskUpc2Pi2E(const char *name, Bool_t _isMC)
  : AliAnalysisTaskSE(name),
    fPIDResponse(0), fTrackCutsBit4(0), isMC(0), debugMode(0), isUsingEffi(0), fTriggerName(0),
  	f2Pi2ETree(0), f2Pi2ETree1(0), fMCTree(0), fTPCNcls(50), fOption(0),
  	BunchCrossNumber_T(0), OrbitNumber_T(0), PeriodNumber_T(0),
  	RunNum_T(0), LikeSign_T(0), Mass_T(0), Pt_T(0), Rapidity_T(0), V0Adecision_T(0), NTracks_T(0),
  	V0Cdecision_T(0), ADAdecision_T(0), ADCdecision_T(0), UBAfired_T(0), UBCfired_T(0), 
  	VBAfired_T(0), VBCfired_T(0), ZNAenergy_T(0), ZNCenergy_T(0), 
  	ZPAenergy_T(0), ZPCenergy_T(0),VtxContrib_T(0), SpdVtxContrib_T(0),
  	VtxChi2_T(0),VtxNDF_T(0), TriggerTOF_T(0), TriggerSPD_T(0),
  	Ntracklets_T(0), Phi_T(0), ChipCut_T(0), GenPart_T(0),ntrk(0),
  	RunNum_MC_T(0), Mass_MC_T(0), Pt_MC_T(0), Rapidity_MC_T(0), Phi_MC_T(0), 
	fListHist(0),fSPDfile(0), hBCmod4(0), hSPDeff(0), fEfficiencyFileName(0),
	fTOFfile(0), fTOFFileName(0), hTOFeff(0), fTOFmask(0), isUsingTOFeff(0),
	fHistTriggersPerRun(0),fITSmodule(0),fFOchip(0),fFOcount(0),TPCclustersP(0),
	TPCclustersN(0),dEdx(0),EtaPhiP(0),EtaPhiN(0), fFOcorr(0), fGoodTracks(0), fTrackChi2(0),
	fHistdEdxVsP1(0),fHistdEdxVsP2(0),fHistdEdxVsP3(0),fHistdEdxVsP4(0),fHistdEdxVsP5(0),
	fHistdEdxVsP6(0),fHistdEdxVsP7(0),fHistdEdxVsP8(0),fHistdEdxVsP9(0), fDeltaPhiRho(0), 
	fDeltaPhiEe(0), FORChip(0),Trigger_T(0)
{
  if(debugMode) std::cout<<"Initialization..."<<std::endl;
  Init();
  if(debugMode) std::cout<<"Defining output..."<<std::endl;
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
//  DefineOutput(3, TTree::Class());
  if (_isMC) DefineOutput(4, TTree::Class());
  if(debugMode) std::cout<<"Initialization done."<<std::endl;
}

AliAnalysisTaskUpc2Pi2E::~AliAnalysisTaskUpc2Pi2E() 
{
  // Destructor
  if(f2Pi2ETree){
	delete f2Pi2ETree;
	f2Pi2ETree = 0x0;
  }
  if(f2Pi2ETree1){
	delete f2Pi2ETree1;
	f2Pi2ETree1 = 0x0;
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

void AliAnalysisTaskUpc2Pi2E::Init()
{
	Trigger_T=0;
	for (Int_t i=0;i<Maxtrk;i++){
	//for (Int_t i=0;i<4;i++){
		PIDTPCPion_T[i] = -20;
		PIDTPCMuon_T[i] = -20;
		PIDTPCKaon_T[i] = -20;
		PIDTPCProton_T[i] = -20;
		PIDTPCElectron_T[i] = -20;
		TPCsignal_T[i] = -1;
		TPCrefit_T[i] = -1;
		ITSrefit_T[i] = -1;
		ITSI_T[i] = -1;
		ITSO_T[i] = -1;
		ITSSA_T[i] = -1;
		TrackCuts_T[i] = -1;
		ITSNcls_T[i] = -1;
		TPCNcls_T[i] = -1;
		TPCchi2_T[i] = -1;
		TrackP_T[i] = -1;
		TrackC_T[i] = 0;
		TrackEta_T[i] = -9;
		TrackPhi_T[i] = -9;
		fMatchingSPD_T[i]=-1;
	}
	for (Int_t i=0;i<3;i++){
		Vertex_T[i] = -99;
		SpdVertex_T[i] = -99;
	}
	for (Int_t i=0;i<4;i++){
		ZDCAtime_T[i] = -99;
		ZDCCtime_T[i] = -99;
	}
	FORChip.reserve(1084);
}

void AliAnalysisTaskUpc2Pi2E::UserCreateOutputObjects() 
{
//	debugMode=1;
	if(debugMode) std::cout<<"Starting UserCreateOutputObjects..."<<std::endl;
  	//PID response
 	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  	AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  	fPIDResponse = inputHandler->GetPIDResponse();
	fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);

  	GenPart_T = new TClonesArray("TParticle", 1000);

  	if(debugMode) std::cout<<"Defining ttree..."<<std::endl;
	f2Pi2ETree1 = new TTree("SelectedRhoJpsi","Selected RhoJpsi events");
	f2Pi2ETree1->Branch("RunNum1_T",&RunNum_T,"RunNum_T/I");
	f2Pi2ETree1->Branch("PeriodNumber1_T",&PeriodNumber_T,"PeriodNumber_T/i");
	f2Pi2ETree1->Branch("OrbitNumber1_T",&OrbitNumber_T,"OrbitNumber_T/i");
	f2Pi2ETree1->Branch("BunchCrossNumber1_T",&BunchCrossNumber_T,"BunchCrossNumber_T/s");
//	f4PiTree1->Branch("LikeSign1_T",&LikeSign_T,"LikeSign_T/O");
//	f4PiTree1->Branch("Mass1_T",&Mass_T,"Mass_T/F");
//	f4PiTree1->Branch("Pt1_T",&Pt_T,"Pt_T/F");
//	f4PiTree1->Branch("Rapidity1_T",&Rapidity_T,"Rapidity_T/F");
//	f4PiTree1->Branch("Phi1_T",&Phi_T,"Phi_T/F");
	f2Pi2ETree1->Branch("ZNAenergy1_T",&ZNAenergy_T,"ZNAenergy_T/F");
	f2Pi2ETree1->Branch("ZNCenergy1_T",&ZNCenergy_T,"ZNCenergy_T/F");
	f2Pi2ETree1->Branch("ZPAenergy1_T",&ZPAenergy_T,"ZPAenergy_T/F");
	f2Pi2ETree1->Branch("ZPCenergy1_T",&ZPCenergy_T,"ZPCenergy_T/F");
	f2Pi2ETree1->Branch("ZDCAtime1_T",&ZDCAtime_T,"ZDCAtime_T[4]/F");
	f2Pi2ETree1->Branch("ZDCCtime1_T",&ZDCCtime_T,"ZDCCtime_T[4]/F");
	f2Pi2ETree1->Branch("PIDTPCPion1_T",&PIDTPCPion_T,"PIDTPCPion_T[20]/F");
	f2Pi2ETree1->Branch("PIDTPCMuon1_T",&PIDTPCMuon_T,"PIDTPCMuon_T[20]/F");
	f2Pi2ETree1->Branch("PIDTPCElectron1_T",&PIDTPCElectron_T,"PIDTPCElectron_T[20]/F");
	f2Pi2ETree1->Branch("PIDTPCKaon1_T",&PIDTPCKaon_T,"PIDTPCKaon_T[20]/F");
	f2Pi2ETree1->Branch("PIDTPCProton1_T",&PIDTPCProton_T,"PIDTPCProton_T[20]/F");
	f2Pi2ETree1->Branch("TPCsignal1_T",&TPCsignal_T,"TPCsignal_T[20]/I");
	f2Pi2ETree1->Branch("TPCrefit1_T",&TPCrefit_T,"TPCrefit_T[20]/O");
	f2Pi2ETree1->Branch("TPCNcls1_T",&TPCNcls_T,"TPCNcls_T[20]/I");
	f2Pi2ETree1->Branch("TPCchi21_T",&TPCchi2_T,"TPCchi2_T[20]/F");
	f2Pi2ETree1->Branch("ISTrefit1_T",&ITSrefit_T,"ITSrefit_T[20]/O");
	f2Pi2ETree1->Branch("ITSI1_T",&ITSI_T,"ITSI_T[20]/O");
	f2Pi2ETree1->Branch("ITSO1_T",&ITSO_T,"ITSO_T[20]/O");
	f2Pi2ETree1->Branch("ITSSA1_T",&ITSSA_T,"ITSSA_T[20]/O");
	f2Pi2ETree1->Branch("TrackCuts_T",&TrackCuts_T,"TrackCuts_T[20]/O");
	f2Pi2ETree1->Branch("ITSNcls1_T",&ITSNcls_T,"ITSNcls_T[20]/I");
	f2Pi2ETree1->Branch("TrackP1_T",&TrackP_T,"TrackP_T[20]/F");
	f2Pi2ETree1->Branch("TrackDCAxy1_T",&TrackDCAxy_T,"TrackDCAxy_T[20]/F");
	f2Pi2ETree1->Branch("TrackDCAz1_T",&TrackDCAz_T,"TrackDCAz_T[20]/F");
	f2Pi2ETree1->Branch("TrackC1_T",&TrackC_T,"TrackC_T[20]/I");
	f2Pi2ETree1->Branch("TrackEta1_T",&TrackEta_T,"TrackEta_T[20]/F");
	f2Pi2ETree1->Branch("TrackPhi1_T",&TrackPhi_T,"TrackPhi_T[20]/F");
	f2Pi2ETree1->Branch("Trigger1_T",&Trigger_T,"Trigger_T/I");
//	f4PiTree1->Branch("TrackPx1_T",&TrackPx_T,"TrackPx_T[8]/F");
//	f4PiTree1->Branch("TrackPy1_T",&TrackPy_T,"TrackPy_T[8]/F");
//	f4PiTree1->Branch("TrackPz1_T",&TrackPz_T,"TrackPz_T[8]/F");
	f2Pi2ETree1->Branch("VtxX1_T",&Vertex_T[0],"VtxX_T/F");
	f2Pi2ETree1->Branch("VtxY1_T",&Vertex_T[1],"VtxY_T/F");
	f2Pi2ETree1->Branch("VtxZ1_T",&Vertex_T[2],"VtxZ_T/F");
//	f4PiTree1->Branch("VtxContrib1_T",&VtxContrib_T,"VtxContrib_T/I");
//	f4PiTree1->Branch("VtxChi21_T",&VtxChi2_T,"VtxChi2_T/F");
//	f4PiTree1->Branch("VtxNDF1_T",&VtxNDF_T,"VtxNDF_T/F");
//	f2Pi2ETree1->Branch("SpdVtxX1_T",&SpdVertex_T[0],"SpdVtxX_T/F");
//	f2Pi2ETree1->Branch("SpdVtxY1_T",&SpdVertex_T[1],"SpdVtxY_T/F");
//	f2Pi2ETree1->Branch("SpdVtxZ1_T",&SpdVertex_T[2],"SpdVtxZ_T/F");
	f2Pi2ETree1->Branch("SpdVtxContrib1_T",&SpdVtxContrib_T,"SpdVtxContrib_T/I");
	f2Pi2ETree1->Branch("V0Adecision1_T",&V0Adecision_T,"V0Adecision_T/I");
	f2Pi2ETree1->Branch("V0Cdecision1_T",&V0Cdecision_T,"V0Cdecision_T/I");
	f2Pi2ETree1->Branch("ADAdecision1_T",&ADAdecision_T,"ADAdecision_T/I");
	f2Pi2ETree1->Branch("ADCdecision1_T",&ADCdecision_T,"ADCdecision_T/I");
	f2Pi2ETree1->Branch("NTracks1_T",&NTracks_T,"NTracks_T/I");
	f2Pi2ETree1->Branch("fMatchingSPD_T", &fMatchingSPD_T, "fMatchingSPD_T[20]/O");
//	f4PiTree1->Branch("UBAfired1_T",&UBAfired_T,"UBAfired_T/O");
//	f4PiTree1->Branch("UBCfired1_T",&UBCfired_T,"UBCfired_T/O");
//	f4PiTree1->Branch("VBAfired1_T",&VBAfired_T,"VBAfired_T/O");
//	f4PiTree1->Branch("VBCfired1_T",&VBCfired_T,"VBCfired_T/O");
//	f4PiTree1->Branch("Ntracklets1_T",&Ntracklets_T,"Ntracklets_T/I");
	f2Pi2ETree1->Branch("ITSModuleInner_T",&ITSModuleInner_T,"ITSModuleInner_T[20]/I");
	f2Pi2ETree1->Branch("ITSModuleOuter_T",&ITSModuleOuter_T,"ITSModuleOuter_T[20]/I");
	f2Pi2ETree1->Branch("FOmodules_T",&fFOmodules_T,"FOmodules_T[240]/I");
	f2Pi2ETree1->Branch("ChipCut1_T",&ChipCut_T,"ChipCut_T/O");
	//f4PiTree1->Branch("FORChip1_T",&FORChip);
//	f4PiTree1->Branch("TriggerSPD1_T",&TriggerSPD_T,"TriggerSPD_T/O");
//	f4PiTree1->Branch("TriggerTOF1_T",&TriggerTOF_T,"TriggerTOF_T/O");

/*	f2Pi2ETree = new TTree("Selected","Selected Rho0 events");
	//define branches
	f2Pi2ETree->Branch("RunNum_T",&RunNum_T,"RunNum_T/I");
	f2Pi2ETree->Branch("PeriodNumber_T",&PeriodNumber_T,"PeriodNumber_T/i");
	f2Pi2ETree->Branch("OrbitNumber_T",&OrbitNumber_T,"OrbitNumber_T/i");
	f2Pi2ETree->Branch("BunchCrossNumber_T",&BunchCrossNumber_T,"BunchCrossNumber_T/s");
	f2Pi2ETree->Branch("LikeSign_T",&LikeSign_T,"LikeSign_T/O");
	f2Pi2ETree->Branch("Mass_T",&Mass_T,"Mass_T/F");
	f2Pi2ETree->Branch("Pt_T",&Pt_T,"Pt_T/F");
	f2Pi2ETree->Branch("Rapidity_T",&Rapidity_T,"Rapidity_T/F");
	f2Pi2ETree->Branch("Phi_T",&Phi_T,"Phi_T/F");
	f2Pi2ETree->Branch("ZNAenergy_T",&ZNAenergy_T,"ZNAenergy_T/F");
	f2Pi2ETree->Branch("ZNCenergy_T",&ZNCenergy_T,"ZNCenergy_T/F");
	f2Pi2ETree->Branch("ZPAenergy_T",&ZPAenergy_T,"ZPAenergy_T/F");
	f2Pi2ETree->Branch("ZPCenergy_T",&ZPCenergy_T,"ZPCenergy_T/F");
	f2Pi2ETree->Branch("ZDCAtime_T",&ZDCAtime_T,"ZDCAtime_T[4]/F");
	f2Pi2ETree->Branch("ZDCCtime_T",&ZDCCtime_T,"ZDCCtime_T[4]/F");
	f2Pi2ETree->Branch("PIDTPCPion_T",&PIDTPCPion_T,"PIDTPCPion_T[18]/F");
	f2Pi2ETree->Branch("PIDTPCElectron_T",&PIDTPCElectron_T,"PIDTPCElectron_T[18]/F");
	f2Pi2ETree->Branch("TPCsignal_T",&TPCsignal_T,"TPCsignal_T[18]/I");
	f2Pi2ETree->Branch("TPCrefit_T",&TPCrefit_T,"TPCrefit_T[18]/O");
	f2Pi2ETree->Branch("TPCNcls_T",&TPCNcls_T,"TPCNcls_T[18]/I");
	f2Pi2ETree->Branch("TPCchi2_T",&TPCchi2_T,"TPCchi2_T[18]/F");
	f2Pi2ETree->Branch("ISTrefit_T",&ITSrefit_T,"ITSrefit_T[18]/O");
	f2Pi2ETree->Branch("NTracks_T",&NTracks_T,"NTracks_T/I");
	f2Pi2ETree->Branch("TrackP_T",&TrackP_T,"TrackP_T[18]/F");
	f2Pi2ETree->Branch("TrackEta_T",&TrackEta_T,"TrackEta_T[18]/F");
	f2Pi2ETree->Branch("TrackPhi_T",&TrackPhi_T,"TrackPhi_T[18]/F");
	f2Pi2ETree->Branch("TrackPx_T",&TrackPx_T,"TrackPx_T[18]/F");
	f2Pi2ETree->Branch("TrackPy_T",&TrackPy_T,"TrackPy_T[18]/F");
	f2Pi2ETree->Branch("TrackPz_T",&TrackPz_T,"TrackPz_T[18]/F");
	f2Pi2ETree->Branch("VtxX_T",&Vertex_T[0],"VtxX_T/F");
	f2Pi2ETree->Branch("VtxY_T",&Vertex_T[1],"VtxY_T/F");
	f2Pi2ETree->Branch("VtxZ_T",&Vertex_T[2],"VtxZ_T/F");
	f2Pi2ETree->Branch("VtxContrib_T",&VtxContrib_T,"VtxContrib_T/I");
	f2Pi2ETree->Branch("VtxChi2_T",&VtxChi2_T,"VtxChi2_T/F");
	f2Pi2ETree->Branch("VtxNDF_T",&VtxNDF_T,"VtxNDF_T/F");
	f2Pi2ETree->Branch("SpdVtxX_T",&SpdVertex_T[0],"SpdVtxX_T/F");
	f2Pi2ETree->Branch("SpdVtxY_T",&SpdVertex_T[1],"SpdVtxY_T/F");
	f2Pi2ETree->Branch("SpdVtxZ_T",&SpdVertex_T[2],"SpdVtxZ_T/F");
	f2Pi2ETree->Branch("SpdVtxContrib_T",&SpdVtxContrib_T,"SpdVtxContrib_T/I");
	f2Pi2ETree->Branch("V0Adecision_T",&V0Adecision_T,"V0Adecision_T/I");
	f2Pi2ETree->Branch("V0Cdecision_T",&V0Cdecision_T,"V0Cdecision_T/I");
	f2Pi2ETree->Branch("ADAdecision_T",&ADAdecision_T,"ADAdecision_T/I");
	f2Pi2ETree->Branch("ADCdecision_T",&ADCdecision_T,"ADCdecision_T/I");
	f2Pi2ETree->Branch("UBAfired_T",&UBAfired_T,"UBAfired_T/O");
	f2Pi2ETree->Branch("UBCfired_T",&UBCfired_T,"UBCfired_T/O");
	f2Pi2ETree->Branch("VBAfired_T",&VBAfired_T,"VBAfired_T/O");
	f2Pi2ETree->Branch("VBCfired_T",&VBCfired_T,"VBCfired_T/O");
	f2Pi2ETree->Branch("Ntracklets_T",&Ntracklets_T,"Ntracklets_T/I");
	// fRhoTree->Branch("ITSModule_T",&ITSModule_T,"ITSModule_T/I");
	f2Pi2ETree->Branch("ChipCut_T",&ChipCut_T,"ChipCut_T/O");
	f2Pi2ETree->Branch("TriggerSPD_T",&TriggerSPD_T,"TriggerSPD_T/O");
	f2Pi2ETree->Branch("TriggerTOF_T",&TriggerTOF_T,"TriggerTOF_T/O");
	//f4PiTree->Branch("FORChip_T",&FORChip);
*/
	if(debugMode) std::cout<<"Defining MC ttree..."<<std::endl;
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

	if(debugMode) std::cout<<"Defining TList..."<<std::endl;
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
	// phitest
	fDeltaPhiRho = new TH1F("fDeltaPhiRho","fDeltaPhiRho",100,-6,6); fListHist->Add(fDeltaPhiRho);
	fDeltaPhiEe = new TH1F("fDeltaPhiEe","fDeltaPhiEe",100,-6,6); fListHist->Add(fDeltaPhiEe);

	fHistdEdxVsP1 = new TH2F("fHistdEdxVsP_Electron","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP1);
	fHistdEdxVsP2 = new TH2F("fHistdEdxVsP_Muon","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP2);
	fHistdEdxVsP3 = new TH2F("fHistdEdxVsP_Pion","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP3);
	fHistdEdxVsP4 = new TH2F("fHistdEdxVsP_Kaon","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP4);
	fHistdEdxVsP5 = new TH2F("fHistdEdxVsP_Proton","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP5);
	fHistdEdxVsP6 = new TH2F("fHistdEdxVsP_Deuteron","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP6);
	fHistdEdxVsP7 = new TH2F("fHistdEdxVsP_Triton","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP7);
	fHistdEdxVsP8 = new TH2F("fHistdEdxVsP_He3","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP8);
	fHistdEdxVsP9 = new TH2F("fHistdEdxVsP_Alpha","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
	fListHist->Add(fHistdEdxVsP9);


	// load SPD effi
	if (isUsingEffi) {
		if(debugMode) std::cout<<"Using efficiency file: "<<fEfficiencyFileName<<std::endl;
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

	// load TOF effi
	if (isUsingTOFeff) {
		if(debugMode) std::cout<<"Using TOF efficiency file: "<<fTOFFileName<<std::endl;
		if(!fTOFfile)fTOFfile = AliDataFile::OpenOADB(fTOFFileName.Data());
	    AliOADBContainer* fTOFcont = (AliOADBContainer*)fTOFfile->Get("TOFTriggerEfficiency");
	    hTOFeff  = (TH2F*)fTOFcont->GetObject(280235,"Default");
	    hTOFeff->SetDirectory(0);
	    fTOFfile->Close();
	}

	if(debugMode) std::cout<<"Post data..."<<std::endl;
	PostData(1, f2Pi2ETree1);
	PostData(2, fListHist);
//	PostData(3, f4PiTree1);
	if (isMC) PostData(4, fMCTree);
	if(debugMode) std::cout<<"Post data done."<<std::endl;
}

void AliAnalysisTaskUpc2Pi2E::UserExec(Option_t *) 
{
  	TDatabasePDG *pdgdat = TDatabasePDG::Instance(); 
  	TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  	Double_t pionMass = partPion->Mass();
  // if(debugMode) std::cout<<"Starting UserExec..."<<std::endl;
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
  //if(!isMC && trigger.Contains(fTriggerName.Data())) {
  //	if(debugMode)std::cout<<trigger<<std::endl;
 // }

  // CCUP9-B - *0VBA *0VBC *0UBA *0UBC 0STP
	Trigger_T=0;
  if (!isMC) { // data
	if (esd->GetRunNumber() < 295828){
  	if (trigger.Contains("CCUP9-B")) Trigger_T = 2;
//  	else if (trigger.Contains("CCUP8")) Trigger_T = 1;
	}
      //if(esd->GetRunNumber()>=295881 && esd->GetRunNumber()<296594)
      if(esd->GetRunNumber()>=295881)
      //if(trigger.Contains("CCUP29-B-SPD2-CENTNOTRD") || trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") || trigger.Contains("CCUP31-B-SPD2-CENTNOTRD")) Trigger_T = 3;
      if(trigger.Contains("CCUP31-B-SPD2-CENTNOTRD")) Trigger_T = 3;
//    if(esd->GetRunNumber()>=296594)
      //if(trigger.Contains("CCUP29-U-SPD2-CENTNOTRD") || trigger.Contains("CCUP30-B-SPD2-CENTNOTRD")  || trigger.Contains("CCUP31-B-SPD2-CENTNOTRD")) Trigger_T = 4;
    //if(esd->GetRunNumber()<295881 && esd->GetRunNumber() > 295828)
    if(esd->GetRunNumber()<295881)
      //if(trigger.Contains("CCUP29-B-NOPF-CENTNOTRD") || trigger.Contains("CCUP30-B-NOPF-CENTNOTRD") || trigger.Contains("CCUP31-B-NOPF-CENTNOTRD")) Trigger_T = 5;
      if(trigger.Contains("CCUP31-B-NOPF-CENTNOTRD")) Trigger_T = 3;
	if (!Trigger_T) return;
	else fHistTriggersPerRun->Fill(esd->GetRunNumber());
  }
  else { // MC
  	if (!IsTriggered(esd)) return;
  } // end of MC trigger
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

	//SPD primary vertex
//	AliESDVertex *fSPDVertex = (AliESDVertex*) esd->GetPrimaryVertexSPD();
//	SpdVtxContrib_T = fSPDVertex->GetNContributors();
//	SpdVertex_T[0] = fSPDVertex->GetX();
//	SpdVertex_T[1] = fSPDVertex->GetY();
//	SpdVertex_T[2] = fSPDVertex->GetZ();

	// Tracklets
	Ntracklets_T = esd->GetMultiplicity()->GetNumberOfTracklets();

	Double_t PhiPlus = 0;
	// loop over four good tracks

  AliESDVertex *fESDVertex = (AliESDVertex*) esd->GetPrimaryVertex();
	// primary vertex
	VtxContrib_T = fESDVertex->GetNContributors();
	Vertex_T[0] = fESDVertex->GetX();
	Vertex_T[1] = fESDVertex->GetY();
	Vertex_T[2] = fESDVertex->GetZ();
	VtxChi2_T = fESDVertex->GetChi2();
	VtxNDF_T = fESDVertex->GetNDF();

if (fESDVertex->GetNContributors()<2) return;

  Int_t nGoodTracks=0;
  Int_t TrackIndex[Maxtrk] = {0};
//  Int_t TrackIndex[4];
for (Int_t it=0;it < Maxtrk; it++){
	TrackIndex[it]=-1;
}
  if(debugMode)std::cout<<"starting track loop"<<std::endl;
  //Track loop - cuts
  for(Int_t itr=0; itr<esd ->GetNumberOfTracks(); itr++) {
    AliESDtrack *trk = esd->GetTrack(itr);
    // geometrical cut
//    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
//	if (fOption.Contains("GeoCut")) esdTrackCuts->SetCutGeoNcrNcl(3.,130.,1.5,0.85,0.7);

    if( !trk ) continue;
 	if( trk->IsOn(AliESDtrack::kITSpureSA) ) continue;
//    if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
    if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
//    if(trk->GetTPCNcls() < fTPCNcls)continue;
//    if(trk->GetTPCchi2()/trk->GetTPCNcls() > 4)continue;
//    if(!((trk->HasPointOnITSLayer(0))&&(trk->HasPointOnITSLayer(1)))) continue;
    Float_t dca[2] = {0.0,0.0}; AliExternalTrackParam cParam;
    if(!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(),300.,&cParam)) continue;
    trk->GetImpactParameters(dca[0],dca[1]);
    //if(TMath::Abs(dca[1]) > 2) continue;
    if(TMath::Abs(dca[1]) > 3) continue;
    Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
//    if(TMath::Abs(dca[0]) > cut_DCAxy) continue;
    if(TMath::Abs(dca[0]) > 3) continue;

	// store good track index
	TrackIndex[nGoodTracks] = itr;
	TPCrefit_T[nGoodTracks] = trk->GetStatus() & AliESDtrack::kTPCrefit;
	ITSrefit_T[nGoodTracks] = trk->GetStatus() & AliESDtrack::kITSrefit;
	TPCNcls_T[nGoodTracks] = trk->GetTPCNcls();
	TPCchi2_T[nGoodTracks] = trk->GetTPCchi2();
	ITSI_T[nGoodTracks] = trk->HasPointOnITSLayer(0);
	ITSO_T[nGoodTracks] = trk->HasPointOnITSLayer(1);
	ITSSA_T[nGoodTracks] = trk->IsPureITSStandalone();
	ITSNcls_T[nGoodTracks] = trk->GetNumberOfITSClusters();
	TrackCuts_T[nGoodTracks] = fTrackCutsBit4->AcceptTrack(trk);
	TrackDCAxy_T[nGoodTracks] = dca[0];
	TrackDCAz_T[nGoodTracks] = dca[1];
	Int_t crossedFO[4];
	TBits fFOCrossedChips(1200); 
	const AliVMultiplicity *mult = esd->GetMultiplicity();
	TBits fFOFiredChips = mult->GetFastOrFiredChips();

	fFOCrossedChips.ResetAllBits(kFALSE);
        crossedFO[0] = trk->GetITSModuleIndex(0);
        crossedFO[1] = trk->GetITSModuleIndex(1);
        crossedFO[2] = trk->GetITSModuleIndex(6);
        crossedFO[3] = trk->GetITSModuleIndex(7);
        SetCrossed(crossedFO, fFOCrossedChips);
	fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
	fMatchingSPD_T[nGoodTracks] = IsSTGFired(fFOCrossFiredChips,RunNum_T >= 295753 ? 9 : 3);
	
	ntrk++;
    nGoodTracks++;
    // std::cout<<nGoodTracks<<" good tracks"<<endl;
	if(nGoodTracks > Maxtrk) break; // just to know how many nGoodTrack are there

  }//Track loop end

  	fGoodTracks->Fill(nGoodTracks);
	NTracks_T = nGoodTracks;
  if(nGoodTracks > 1 && nGoodTracks < Maxtrk){ // fill tree variables
	if(debugMode)std::cout<<"two or four good tracks"<<std::endl;

  	Double_t charge[Maxtrk]; // charge
	TLorentzVector lv[Maxtrk];
//  	Double_t charge[4]; // charge
//	TLorentzVector lv[4];
	TLorentzVector lvSum; // pair-4vector

//	Int_t fFOmodules[240];
	for (Int_t i=0;i<240;i++) fFOmodules_T[i] = 0;

  	for(Int_t chipkey=0;chipkey<1200;chipkey++){
  		if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)){
  			fFOmodules_T[(chipkey/5)]++;
 		}
 	}

  	for(Int_t i=0; i<nGoodTracks; i++){				//CHANGE i here to look at different track number
	//	if (TrackIndex[i] <= 0) continue; 
  if(debugMode) std::cout<<"0"<<std::endl;
	  	AliESDtrack *trk = esd->GetTrack(TrackIndex[i]);
//		if (!trk) continue;
		ITSModuleInner_T[i] = trk->GetITSModuleIndex(0)/1000000;
		ITSModuleOuter_T[i] = trk->GetITSModuleIndex(1)/1000000;

		// phi test
		if (trk->Charge()>0) PhiPlus = trk->Phi();

  if(debugMode) std::cout<<"1"<<std::endl;
		// TPC PID n-sigma
		PIDTPCElectron_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
		PIDTPCPion_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
		PIDTPCMuon_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
		PIDTPCKaon_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
		PIDTPCProton_T[i] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
  if(debugMode) std::cout<<"2"<<std::endl;

		// separated PID
		Double_t ptrackTPC=-1.;
		const AliExternalTrackParam* ippar=trk->GetInnerParam();
		if(ippar) ptrackTPC=ippar->P();
		Double_t dedx=trk->GetTPCsignal();
		Int_t  pidtr=trk->GetPIDForTracking();
  if(debugMode) std::cout<<"3"<<std::endl;
		if(pidtr==0) fHistdEdxVsP1->Fill(ptrackTPC,dedx);
		if(pidtr==1) fHistdEdxVsP2->Fill(ptrackTPC,dedx);
		if(pidtr==2) fHistdEdxVsP3->Fill(ptrackTPC,dedx);
		if(pidtr==3) fHistdEdxVsP4->Fill(ptrackTPC,dedx);
		if(pidtr==4) fHistdEdxVsP5->Fill(ptrackTPC,dedx);
		if(pidtr==5) fHistdEdxVsP6->Fill(ptrackTPC,dedx);
		if(pidtr==6) fHistdEdxVsP7->Fill(ptrackTPC,dedx);
		if(pidtr==7) fHistdEdxVsP8->Fill(ptrackTPC,dedx);
		if(pidtr==8) fHistdEdxVsP9->Fill(ptrackTPC,dedx);
  if(debugMode) std::cout<<"4"<<std::endl;

		charge[i] = trk->Charge();
		TPCsignal_T[i] = trk->GetTPCsignal();
		TrackP_T[i] = trk->P();
		TrackC_T[i] = charge[i];
		TrackPhi_T[i] = trk->Phi();
		TrackEta_T[i] = trk->Eta();
		TrackPx_T[i] = trk->Px();
		TrackPy_T[i] = trk->Py();
		TrackPz_T[i] = trk->Pz();
  if(debugMode) std::cout<<"5"<<std::endl;

		fTrackChi2->Fill((Float_t)trk->GetTPCchi2()/trk->GetTPCNcls());

//		fITSmodule->Fill(ITSModuleInner_T[i]);
//		fITSmodule->Fill(ITSModuleOuter_T[i]);
  if(debugMode) std::cout<<"6"<<std::endl;

//		for(Int_t j=0;j<240;j++){
//			if (fFOmodules[j] > 0){
//				fFOcorr->Fill(ITSModuleInner_T[i],j);
//				fFOcorr->Fill(ITSModuleOuter_T[i],j);
//			}
//		}

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

  	} // end loop over four good tracks
  if(debugMode) std::cout<<"00"<<std::endl;

  	lvSum = lv[0]+lv[1]+lv[2]+lv[3];

	if (charge[0]+charge[1]+charge[2]+charge[3]!=0) LikeSign_T = 1;
	else LikeSign_T = 0;
	Mass_T = lvSum.M();
	Pt_T = lvSum.Pt();
	Rapidity_T = lvSum.Rapidity();
	Phi_T = lvSum.Phi();
	// phi test
	if (!LikeSign_T && fabs(Rapidity_T)<0.8 && Pt_T<0.1 && (pow(PIDTPCPion_T[0],2)+pow(PIDTPCPion_T[1],2)+pow(PIDTPCPion_T[2],2)+pow(PIDTPCPion_T[3],2)<18) ) fDeltaPhiRho->Fill(Phi_T-PhiPlus);
	if (!LikeSign_T && fabs(Rapidity_T)<0.8 && Pt_T<0.1 && (pow(PIDTPCElectron_T[0],2)+pow(PIDTPCElectron_T[1],2)+pow(PIDTPCElectron_T[2],2)+pow(PIDTPCElectron_T[3],2)<18) ) fDeltaPhiEe->Fill(Phi_T-PhiPlus);

	// virtual cut on FO chip matching
	Int_t SPDInner[20]; for (Int_t i=0; i<20; ++i) SPDInner[i]=0;
	Int_t SPDOuter[40]; for (Int_t i=0; i<40; ++i) SPDOuter[i]=0;

	SPDInner[ITSModuleInner_T[0]/4]++;
	SPDInner[ITSModuleInner_T[1]/4]++;
	SPDInner[ITSModuleInner_T[2]/4]++;
	SPDInner[ITSModuleInner_T[3]/4]++;
	SPDOuter[(ITSModuleOuter_T[0]-80)/4]++;
	SPDOuter[(ITSModuleOuter_T[1]-80)/4]++;
	SPDOuter[(ITSModuleOuter_T[2]-80)/4]++;
	SPDOuter[(ITSModuleOuter_T[3]-80)/4]++;
	ChipCut_T = 0;
	if ((fTriggerName == "CCUP9-B") &&
//		((fFOmodules_T[ITSModuleInner_T[0]] == 0)||(fFOmodules_T[ITSModuleOuter_T[0]] == 0)
//		||(fFOmodules_T[ITSModuleInner_T[1]] == 0)||(fFOmodules_T[ITSModuleOuter_T[1]] == 0)
		/*||*/ !Is0STPfired(SPDInner,SPDOuter)) ChipCut_T = 1;
//  }
    
//	Int_t fFOcounter = 0;
  //	for(Int_t chipkey=0;chipkey<1200;chipkey++){
  //		if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)){
  //			fFOchip->Fill(chipkey);
  //			fFOcounter++;
//			FORChip.push_back(chipkey);
 //		}
  //	}
//	Int_t crossedFO[4];
//	TBits fFOCrossedChips(1200);
//	const AliVMultiplicity *mult = esd->GetMultiplicity();
//	TBits fFOFiredChips = mult->GetFastOrFiredChips();
//
//	fFOCrossedChips.ResetAllBits(kFALSE);
//		
  //	fFOcount->Fill(fFOcounter);
//  	fFOcorr->Fill();

  //fill
 // f2Pi2ETree ->Fill();
// if (/*nGoodTracks == 4 &&*/ LikeSign_T==0)
 f2Pi2ETree1->Fill();

  } // end four good tracks
  if(debugMode) std::cout<<"saving data"<<std::endl;
  PostData(1, f2Pi2ETree1);
  PostData(2, fListHist);
//  PostData(3, f4PiTree1);
  if (isMC) PostData(4, fMCTree);
  

}//UserExec

// fuction that get two arrays and return if 0STP trigger was fired
Bool_t AliAnalysisTaskUpc2Pi2E::Is0STPfired(Int_t *vPhiInner, Int_t *vPhiOuter) // array 20, 40
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
void AliAnalysisTaskUpc2Pi2E::SetCrossed(Int_t spd[4], TBits &crossed){ 
    // from the macro by MB
    Int_t chipId2;
    for(Int_t iLayer = 0; iLayer < 4; iLayer++)
    if(spd[iLayer] > 0){
        crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); 
        crossed.SetBitNumber(chipId2); 
    }
}
Int_t AliAnalysisTaskUpc2Pi2E::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
    // from the macro by MB
    Int_t status   = (index%1000000)/100000;
    Int_t iModule  = index/1000000;           // 0 - 239
    Int_t iPhi     = iModule/4;               // 0-19 - inner, 20-59 outer
    Int_t iModuleZ = iModule%4;               // 0-3
    Int_t iSign    = (index%100000)/10000;    // 1-4
    Int_t signZ    = iPhi<20 ? (iSign%2==1 ? 1 : -1) : (iSign%2==0 ? 1 : -1); // 1 or -1
    Int_t iX       = (index%10000)/100;       // ??
    Int_t iZ       = index%100;               // 0-36 [mm]
    Int_t signZiZ  = (36-signZ*iZ);
    Int_t chipId   = iModule*5+signZiZ*5/72;
    if (chipId<0) return 1200;
    if (chipId>=1200) return 1201;
    if (signZiZ<0) return 1202;
    if (signZiZ>72) return 1203;
    if (signZiZ==72 && chipId%20==0 && chipId>=400) return 1204;
    chipId2=chipId;

    if (signZiZ==0  && chipId%20!=0)  chipId2=chipId-1;
    if (signZiZ==72 && chipId%20!=19) chipId2=chipId+1;
    if (signZiZ==13)  chipId2=chipId+1;
    if (signZiZ==14)  chipId2=chipId+1;
    if (signZiZ==15)  chipId2=chipId-1;
    if (signZiZ==16)  chipId2=chipId-1;
    if (signZiZ==27)  chipId2=chipId+1;
    if (signZiZ==28)  chipId2=chipId+1;
    if (signZiZ==29)  chipId2=chipId-1;
    if (signZiZ==30)  chipId2=chipId-1;
    if (signZiZ==42)  chipId2=chipId+1;
    if (signZiZ==43)  chipId2=chipId+1;
    if (signZiZ==44)  chipId2=chipId-1;
    if (signZiZ==45)  chipId2=chipId-1;
    if (signZiZ==56)  chipId2=chipId+1;
    if (signZiZ==57)  chipId2=chipId+1;
    if (signZiZ==58)  chipId2=chipId-1;
    if (signZiZ==59)  chipId2=chipId-1;
    if (debug) printf("%4i %4i %3i %3i %3i\n",chipId,chipId2,iX,signZiZ,iSign);
    return chipId;
}
Bool_t AliAnalysisTaskUpc2Pi2E::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
    // from the macro by MB
    Int_t n1 = bits.CountBits(400);
    Int_t n0 = bits.CountBits() - n1;
    //cout<<n0<<" "<<n1<<endl;
    if (n0<1 || n1<1) return 0;
    Bool_t stg = 0;
    Bool_t l0[20]={0};
    Bool_t l1[40]={0};
    Bool_t phi[20]={0};
    for (Int_t i=0;   i< 400; ++i) if (bits.TestBitNumber(i)) l0[      i/20] = 1;
    for (Int_t i=400; i<1200; ++i) if (bits.TestBitNumber(i)) l1[(i-400)/20] = 1;
    for (Int_t i=0; i<20; ++i){
        if (tolerance) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
        else           phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40]);
    }
    for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++) for (Int_t i=0; i<20; ++i) stg |= phi[i] & phi[(i+dphi)%20];
    return stg;
}
Bool_t AliAnalysisTaskUpc2Pi2E::IsTriggered(AliESDEvent *esd)
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
	if (nInner >= 2 && nOuter >= 2) {
		SH1 = kTRUE;
	}

	// V0
	V0A = esd->GetHeader()->IsTriggerInputFired("0VBA");
	V0C = esd->GetHeader()->IsTriggerInputFired("0VBC");
	// AD
	ADA = esd->GetHeader()->IsTriggerInputFired("0UBA");
	ADC = esd->GetHeader()->IsTriggerInputFired("0UBC");
	// TOF
	OMU = esd->GetHeader()->IsTriggerInputFired("0OMU");

	// OM2
	if (isUsingTOFeff) {
	const AliTOFHeader *tofH = esd->GetTOFHeader();
	fTOFmask = tofH->GetTriggerMask();
  
	Bool_t firedMaxiPhi[36] = {0};
	Int_t NfiredMaxiPads = 0;
 
	for(Int_t ltm=0;ltm<72;ltm++){
		Int_t ip = ltm%36;
		for(Int_t cttm=0;cttm<23;cttm++){
			if(fTOFmask->IsON(ltm,cttm) && gRandom->Rndm(1.0)<hTOFeff->GetBinContent(ltm+1,cttm+1)){
				firedMaxiPhi[ip] = kTRUE;
				NfiredMaxiPads++;
			}
		}
	}
	if(NfiredMaxiPads >= 2) {
		OM2 = kTRUE; //0OM2 TOF two hits
	}

	}
	else OM2 = esd->GetHeader()->IsTriggerInputFired("0OM2");

	// save spd and tof trigger decisions to tree
	TriggerSPD_T = SH1;
	TriggerTOF_T = OM2;

	if ((fTriggerName == "CCUP9-B") && (!V0A && !V0C && !ADA && !ADC && STP)) return kTRUE; // CCUP9 is fired
//	else
 //if ((fTriggerName == "CCUP8") && (!V0A && !V0C && !ADA && !ADC && STP && OMU)) return kTRUE; // CCUP8 is fired
// if ((fTriggerName == "CCUP31") && (!V0A && !V0C && !ADA && !ADC && STP && OMU)) return kTRUE; // CCUP8 is fired

	else return kFALSE;
} // end of MC trigger


//_____________________________________________________________________________
void AliAnalysisTaskUpc2Pi2E::Terminate(Option_t *)
{

  std::cout<<"Analysis complete."<<std::endl;
}

