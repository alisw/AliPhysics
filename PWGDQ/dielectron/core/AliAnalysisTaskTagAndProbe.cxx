#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "THnSparse.h"
#include "THashList.h"
#include "TMath.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"

#include "AliMultSelection.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisFilter.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronPID.h"
#include "AliDielectronHistos.h"

#include "AliAnalysisTaskTagAndProbe.h"

//Author: Daiki Sekihata (Center for Nuclear Study, the University of Tokyo)
//daiki.sekihata@cern.ch
//This analysis task is for data-driven efficiency for a single track with tag-and-probe method.
//Since this task aims to measure the efficiency for a single track,
//electrons from a primary vertex or electrons from gamma conversions do not matter.
//In dielectron analyses, a hit on 1st SPD layer is required for electrons.
//Thus, converted electrons are from the beam pipe.
//They are considered as electrons from the primary vertex which have the similar PID/tracking efficiency.

ClassImp(AliAnalysisTaskTagAndProbe)

//________________________________________________________________________
AliAnalysisTaskTagAndProbe::AliAnalysisTaskTagAndProbe():
  AliAnalysisTaskSE(),
  fOutputContainer(0x0),
	fTriggerMask(0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentrality(0),
  fCentralityMin(0.),
  fCentralityMax(101.),
  fNMixed(10),
  fZvtxBin(-1),
	fPIDResponse(0x0),
	fPIDCalibinPU(kTRUE),
	fTagTrackArray(0x0),
	fProbeTrackArray(0x0),
	fPassingProbeTrackArray(0x0),
	fEventFilter(0x0),
	fTagFilter(0x0),
	fProbeFilter(0x0),
	fPassingProbeFilter(0x0),
	fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	fMmax(-1),
	fPhiVmin(3.2)
{
  // Constructor

	for(Int_t i=0;i<3;i++)  fVertex[i] = 0;

	for(Int_t i=0;i<2;i++){
		for(Int_t j=0;j<10;j++){
			fEventList[i][j] = 0x0;
		}
	}

	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<15;j++){
			fPostPIDCntrdCorrPU[i][j] = 0x0;
			fPostPIDWdthCorrPU[i][j]  = 0x0;
		}
	}

}
//________________________________________________________________________
AliAnalysisTaskTagAndProbe::AliAnalysisTaskTagAndProbe(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(0x0),
	fTriggerMask(0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentrality(0),
  fCentralityMin(0.),
  fCentralityMax(101.),
  fNMixed(10),
  fZvtxBin(-1),
	fPIDResponse(0x0),
	fPIDCalibinPU(kTRUE),
	fTagTrackArray(0x0),
	fProbeTrackArray(0x0),
	fPassingProbeTrackArray(0x0),
	fEventFilter(0x0),
	fTagFilter(0x0),
	fProbeFilter(0x0),
	fPassingProbeFilter(0x0),
	fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	fMmax(-1),
	fPhiVmin(3.2)
{
  // Constructor

  for(Int_t i=0;i<3;i++)  fVertex[i] = 0;
	for(Int_t i=0;i<2;i++){
		for(Int_t j=0;j<10;j++){
			fEventList[i][j] = 0x0;
		}
	}

	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<15;j++){
			fPostPIDCntrdCorrPU[i][j] = 0x0;
			fPostPIDWdthCorrPU[i][j]  = 0x0;
		}
	}

	fTagTrackArray = new TObjArray(1000);
	fTagTrackArray->SetOwner(kFALSE);

	fProbeTrackArray = new TObjArray(1000);
	fProbeTrackArray->SetOwner(kFALSE);

	fPassingProbeTrackArray = new TObjArray(1000);
	fPassingProbeTrackArray->SetOwner(kFALSE);

	fEventFilter = new AliAnalysisFilter("fEventFilter","fEventFilter");
	fTagFilter   = new AliAnalysisFilter("fTagFilter"  ,"fTagFilter");
	fProbeFilter = new AliAnalysisFilter("fProbeFilter","fProbeFilter");
	fPassingProbeFilter = new AliAnalysisFilter("fPassingProbeFilter"  ,"fPassingProbeFilter");

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());
}
//________________________________________________________________________
AliAnalysisTaskTagAndProbe::~AliAnalysisTaskTagAndProbe()
{
	delete fEventFilter;
	delete fTagFilter;
	delete fProbeFilter;
	delete fPassingProbeFilter;
	delete fTagTrackArray;
	delete fProbeTrackArray;
	delete fPassingProbeTrackArray;
	delete fUsedVars;
}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  TH1F *hEventSummary = new TH1F("hEventSummary","Event Summary",10,0.5,10.5);
  hEventSummary->GetXaxis()->SetBinLabel(1 ,"all");
  hEventSummary->GetXaxis()->SetBinLabel(2 ,"after P.S.");
  hEventSummary->GetXaxis()->SetBinLabel(3 ,"selected");
  fOutputContainer->Add(hEventSummary);

  //event character histogram
  fOutputContainer->Add(new TH1F("hVertexZ","VertexZ;Zvtx (cm)",1000,-50.,50.));
  fOutputContainer->Add(new TH1F("hVertexZSelectEvent","VertexZ SelectEvent;Zvtx (cm)",1000,-50.,50.));
  fOutputContainer->Add(new TH2F(Form("hCentrality%svsNContributor",fEstimator.Data()),Form("Centrality %s vs. Ncontributor;centrality (%%);N_{contributor}",fEstimator.Data()),101,0.,101,500,0,5000));
  fOutputContainer->Add(new TH2F("hNTPCclsvsNSDDSSDcls","N_{cls}^{TPC} vs. N_{cls}^{SDD+SSD};N_{cls}^{TPC};N_{cls}^{SDD+SSD}",300,0,6e+6,300,0,6e+4));//for pileup plot

	//simple track QA
	const Int_t Ndim    = 3;
	Int_t Nbin[Ndim]    = {100, 20,             36};
	Double_t xmin[Ndim] = {  0, -1,              0};
	Double_t xmax[Ndim] = { 10, +1, TMath::TwoPi()};

	THnSparseF *hs_PtEtaPhi = new THnSparseF("hs_PtEtaPhi","hs_PtEtaPhi;p_{T} (GeV/c);#eta;#varphi (rad.);",Ndim,Nbin,xmin,xmax);
	hs_PtEtaPhi->Sumw2();
  fOutputContainer->Add(hs_PtEtaPhi);

  fOutputContainer->Add(new TH2F("hTrackDCA","DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)",100,-5,+5,100,-5,+5));
  fOutputContainer->Add(new TH1F("hTrackNclsTPC"   ,"Number of clusters TPC;N_{cls}^{TPC}"                ,161,-0.5,160.5));
  fOutputContainer->Add(new TH1F("hTrackNclsPIDTPC","Number of clusters for PID TPC;N_{cls}^{TPC} for PID",161,-0.5,160.5));
  fOutputContainer->Add(new TH1F("hTrackNcrTPC"    ,"Number of crossed rows TPC;N_{cr}^{TPC}"            ,161,-0.5,160.5));
  fOutputContainer->Add(new TH1F("hTrackNfTPC"     ,"Number of findable clusters TPC;N_{f}^{TPC}"       ,161,-0.5,160.5));
  fOutputContainer->Add(new TH1F("hTrackRatioNcrtoNfTPC","ratio of N_{cr}^{TPC}/N_{f}^{TPC};N_{cr}^{TPC}/N_{f}^{TPC}",200,0.,2));
  fOutputContainer->Add(new TH1F("hTrackChi2TPC","chi2 TPC;#chi^{2}/N_{cls}^{TPC}",100,0,10));
  fOutputContainer->Add(new TH1F("hTrackGoldenChi2","golden chi2;golden #chi^{2}" ,100,0,100));

  fOutputContainer->Add(new TH1F("hTrackNclsITS","Number of clusters ITS;N_{cls}^{ITS}"      ,7,-0.5,6.5));
  fOutputContainer->Add(new TH1F("hTrackNscITS" ,"Number of shared clusters ITS;N_{sc}^{ITS}",7,-0.5,6.5));
  fOutputContainer->Add(new TH1F("hTrackChi2ITS","chi2 ITS;#chi^{2}/N_{cls}^{ITS}",100,0,10));

  fOutputContainer->Add(new TH2F("hTrackTPCdEdx","TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx (a.u.)",200,0.,10,200,0,200));
  fOutputContainer->Add(new TH2F("hTrackITSdEdx","ITS dE/dx vs. p_{in};p_{in} (GeV/c);ITS dE/dx (a.u.)",200,0.,10,400,0,400));
  fOutputContainer->Add(new TH2F("hTrackTOFbeta","TOF #beta vs. p_{in};p_{in} (GeV/c);TOF #beta"       ,200,0.,10,240,0,1.2));

  fOutputContainer->Add(new TH2F("hTrackNsigmaElTPCvsPin","n#sigma_{e}^{TPC} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TPC}",200,0.,10,100,-5,+5));
  fOutputContainer->Add(new TH2F("hTrackNsigmaElITSvsPin","n#sigma_{e}^{ITS} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{ITS}",200,0.,10,100,-5,+5));
  fOutputContainer->Add(new TH2F("hTrackNsigmaElTOFvsPin","n#sigma_{e}^{TOF} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TOF}",200,0.,10,100,-5,+5));

  fOutputContainer->Add(new TH2F("hTrackNsigmaElTPCvsEta","n#sigma_{e}^{TPC} vs. #eta;#eta;n#sigma_{e}^{TPC}",200,-1,+1,100,-5,+5));
  fOutputContainer->Add(new TH2F("hTrackNsigmaElITSvsEta","n#sigma_{e}^{ITS} vs. #eta;#eta;n#sigma_{e}^{ITS}",200,-1,+1,100,-5,+5));
  fOutputContainer->Add(new TH2F("hTrackNsigmaElTOFvsEta","n#sigma_{e}^{TOF} vs. #eta;#eta;n#sigma_{e}^{TOF}",200,-1,+1,100,-5,+5));

	//pair histograms
  fOutputContainer->Add(new TH2F("hPairMvsPhiV","m_{ee} vs.#varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})",100,0,TMath::Pi(),100,0,0.1));

	const Int_t Nmee = 150;
	Double_t mee[Nmee] = {};
	for(Int_t i=0  ;i<110 ;i++) mee[i] = 0.01 * (i-  0) + 0.0;//from 0 to 1.1 GeV/c2, every 0.01 GeV/c2
	for(Int_t i=110;i<Nmee;i++) mee[i] = 0.1  * (i-110) + 1.1;//from 1.1 to 5 GeV/c2, evety 0.1 GeV/c2

	const Int_t NpTe = 70;
	Double_t pTe[NpTe] = {};
	for(Int_t i=0  ;i<10  ;i++) pTe[i] = 0.01 * (i- 0) + 0.0;//from 0.0 to 0.1 GeV/c, every 0.01 GeV/c
	for(Int_t i=10 ;i<59;i++)   pTe[i] = 0.1  * (i-10) + 0.1;//from 0.1 to 5.0 GeV/c,evety 0.1 GeV/c
	for(Int_t i=59 ;i<NpTe;i++) pTe[i] = 0.5  * (i-59) + 5.0;//from 5.0 to 10 GeV/c, evety 0.5 GeV/c


	const TString probetype[3]  = {"Probe","PassingProbe"};
	const TString chargetype[3]  = {"ULS","LSpp","LSnn"};
	const TString eventtype[2] = {"same","mix"};

	for(Int_t ip=0;ip<2;ip++){
		for(Int_t ic=0;ic<3;ic++){
			for(Int_t ie=0;ie<2;ie++){
				TH2F *h2TAP = new TH2F(Form("h%s_%s_%s",probetype[ip].Data(),chargetype[ic].Data(),eventtype[ie].Data()),Form("h%s_%s_%s",probetype[ip].Data(),chargetype[ic].Data(),eventtype[ie].Data()),Nmee-1,mee,NpTe-1,pTe);
				h2TAP->SetXTitle("m_{ee} (GeV/c^{2})");
				h2TAP->SetYTitle("p_{T,e} (GeV/c)");
				h2TAP->Sumw2();
				fOutputContainer->Add(h2TAP);
			}
		}
	}



//	const Int_t Ndim_TAP        = 4;//mee, pT_e, eta_e, phi_e
//	Int_t Nbin_TAP[Ndim_TAP]    = {Nmee-1, NpTe-1, 20,             18};
//	Double_t xmin_TAP[Ndim_TAP] = {     0,      0, -1,              0};
//	Double_t xmax_TAP[Ndim_TAP] = {     5,     10, +1, TMath::TwoPi()};
//
//	THnSparseF *hsProbe_ULS_same         = new THnSparseF("hsProbe_ULS_same"        ,"hsProbe_ULS_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"        ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_ULS_same  = new THnSparseF("hsPassingProbe_ULS_same" ,"hsPassingProbe_ULS_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);" ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsProbe_LSpp_same        = new THnSparseF("hsProbe_LSpp_same"       ,"hsProbe_LSpp_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"       ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_LSpp_same = new THnSparseF("hsPassingProbe_LSpp_same","hsPassingProbe_LSpp_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);",Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsProbe_LSnn_same        = new THnSparseF("hsProbe_LSnn_same"       ,"hsProbe_LSnn_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"       ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_LSnn_same = new THnSparseF("hsPassingProbe_LSnn_same","hsPassingProbe_LSnn_same;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);",Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//
//	THnSparseF *hsProbe_ULS_mix         = new THnSparseF("hsProbe_ULS_mix"        ,"hsProbe_ULS_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"        ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_ULS_mix  = new THnSparseF("hsPassingProbe_ULS_mix" ,"hsPassingProbe_ULS_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);" ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsProbe_LSpp_mix        = new THnSparseF("hsProbe_LSpp_mix"       ,"hsProbe_LSpp_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"       ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_LSpp_mix = new THnSparseF("hsPassingProbe_LSpp_mix","hsPassingProbe_LSpp_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);",Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsProbe_LSnn_mix        = new THnSparseF("hsProbe_LSnn_mix"       ,"hsProbe_LSnn_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);"       ,Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//	THnSparseF *hsPassingProbe_LSnn_mix = new THnSparseF("hsPassingProbe_LSnn_mix","hsPassingProbe_LSnn_mix;m_{ee} (GeV/c^{2});p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);",Ndim_TAP,Nbin_TAP,xmin_TAP,xmax_TAP);
//
//	hsProbe_ULS_same        ->Sumw2();
//	hsPassingProbe_ULS_same ->Sumw2();
//	hsProbe_LSpp_same       ->Sumw2();
//	hsPassingProbe_LSpp_same->Sumw2();
//	hsProbe_LSnn_same       ->Sumw2();
//	hsPassingProbe_LSnn_same->Sumw2();
//	hsProbe_ULS_mix         ->Sumw2();
//	hsPassingProbe_ULS_mix  ->Sumw2();
//	hsProbe_LSpp_mix        ->Sumw2();
//	hsPassingProbe_LSpp_mix ->Sumw2();
//	hsProbe_LSnn_mix        ->Sumw2();
//	hsPassingProbe_LSnn_mix ->Sumw2();
//
//	hsProbe_ULS_same->SetBinEdges(0,mee);
//	hsProbe_ULS_same->SetBinEdges(1,pTe);
//	hsPassingProbe_ULS_same->SetBinEdges(0,mee);
//	hsPassingProbe_ULS_same->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_ULS_same);
//  fOutputContainer->Add(hsPassingProbe_ULS_same);
//
//	hsProbe_LSpp_same->SetBinEdges(0,mee);
//	hsProbe_LSpp_same->SetBinEdges(1,pTe);
//	hsPassingProbe_LSpp_same->SetBinEdges(0,mee);
//	hsPassingProbe_LSpp_same->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_LSpp_same);
//  fOutputContainer->Add(hsPassingProbe_LSpp_same);
//
//	hsProbe_LSnn_same->SetBinEdges(0,mee);
//	hsProbe_LSnn_same->SetBinEdges(1,pTe);
//	hsPassingProbe_LSnn_same->SetBinEdges(0,mee);
//	hsPassingProbe_LSnn_same->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_LSnn_same);
//  fOutputContainer->Add(hsPassingProbe_LSnn_same);
//
//	hsProbe_ULS_mix->SetBinEdges(0,mee);
//	hsProbe_ULS_mix->SetBinEdges(1,pTe);
//	hsPassingProbe_ULS_mix->SetBinEdges(0,mee);
//	hsPassingProbe_ULS_mix->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_ULS_mix);
//  fOutputContainer->Add(hsPassingProbe_ULS_mix);
//
//	hsProbe_LSpp_mix->SetBinEdges(0,mee);
//	hsProbe_LSpp_mix->SetBinEdges(1,pTe);
//	hsPassingProbe_LSpp_mix->SetBinEdges(0,mee);
//	hsPassingProbe_LSpp_mix->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_LSpp_mix);
//  fOutputContainer->Add(hsPassingProbe_LSpp_mix);
//
//	hsProbe_LSnn_mix->SetBinEdges(0,mee);
//	hsProbe_LSnn_mix->SetBinEdges(1,pTe);
//	hsPassingProbe_LSnn_mix->SetBinEdges(0,mee);
//	hsPassingProbe_LSnn_mix->SetBinEdges(1,pTe);
//  fOutputContainer->Add(hsProbe_LSnn_mix);
//  fOutputContainer->Add(hsPassingProbe_LSnn_mix);

  PostData(1,fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

	if(!fESDEvent && !fAODEvent){
		AliInfo("event class is neither ESD nor AOD. return.");
		return;
	}

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  FillHistogramTH1(fOutputContainer,"hEventSummary",1);//all
  FillHistogramTH1(fOutputContainer,"hVertexZ" ,fVertex[2]);
  Int_t Ncontributor  = vVertex->GetNContributors();

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isPhysSelOK = fSelectMask & fTriggerMask;
	if(!isPhysSelOK){
		AliInfo("event is rejected by physics selection");
		return;
	}
  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//after physics selection

	UInt_t selectedMask_ev = (1<<fEventFilter->GetCuts()->GetEntries())-1;
	UInt_t cutmask_ev = fEventFilter->IsSelected(InputEvent());
	if(cutmask_ev != selectedMask_ev){
		AliInfo("event is rejected by event filter. return.");
		return;
	}

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliFatal("fPIDResponse does not exist!");
    return;
  }

	AliDielectronVarManager::SetPIDResponse( fInputHandler->GetPIDResponse() );
	AliDielectronVarManager::SetFillMap(fUsedVars);
	AliDielectronVarManager::SetEvent( InputEvent() );


//	if(!fEventFilter->IsSelected(InputEvent())){
//		AliInfo("event is rejected by event filter. return.");
//		return;
//	}

	AliDielectronPID::SetPIDCalibinPU(fPIDCalibinPU);
	for(Int_t id=0;id<15;id++){//detector loop TPC/ITS/TOF
		for(Int_t ip=0;ip<15;ip++){//particle loop e/mu/pi/k/p
			if(fPostPIDCntrdCorrPU[id][ip]) AliDielectronPID::SetCentroidCorrFunctionPU(id,ip,fPostPIDCntrdCorrPU[id][ip]);
			if(fPostPIDWdthCorrPU[id][ip])  AliDielectronPID::SetWidthCorrFunctionPU(   id,ip,fPostPIDWdthCorrPU[id][ip] );
		}
	}

	fTagTrackArray->Clear();
	fProbeTrackArray->Clear();
	fPassingProbeTrackArray->Clear();

	//Double_t values[AliDielectronVarManager::kNMaxValues] = {0.};
	//if(fESDEvent)      AliDielectronVarManager::Fill(fESDEvent,values);
	//else if(fAODEvent) AliDielectronVarManager::Fill(fAODEvent,values);

  fZvtxBin = (Int_t)((fVertex[2]+10.)/2.);//it should be 0-9.
  if(fZvtxBin < 0) fZvtxBin = 0;//protection to avoid fZvtxBin = -1.
  if(fZvtxBin > 9) fZvtxBin = 9;//protection to avoid fZvtxBin = 10.

  //centrality estimation
	//Get Centrality
	fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
	if(!fMultSelection){
		//If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
		AliWarning("AliMultSelection object not found! centrality is set to 0.");
		fCentrality = 0;
		//return;
	}
	else{
		fCentrality = fMultSelection->GetMultiplicityPercentile(fEstimator);
	}
	//printf("centrality NEW = %f\n",values[AliDielectronVarManager::kCentralityNew]);
	//printf("centrality V0M = %f\n",fCentrality);

  FillHistogramTH1(fOutputContainer,"hEventSummary",3);//selected
  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNContributor",fEstimator.Data()),fCentrality,Ncontributor);

	//Int_t NclsSPD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(0);
	//Int_t NclsSPD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(1);
	Int_t NclsSDD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(2);
	Int_t NclsSDD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(3);
	Int_t NclsSSD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(4);
	Int_t NclsSSD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(5);

	Int_t NclsTPC  = 0;
	if(fESDEvent) NclsTPC = fESDEvent->GetNumberOfTPCClusters();
	else if(fAODEvent) NclsTPC = fAODEvent->GetNumberOfTPCClusters();

	FillHistogramTH2(fOutputContainer,"hNTPCclsvsNSDDSSDcls",NclsTPC,NclsSDD0+NclsSDD1+NclsSSD0+NclsSSD1);

	if(!fEventList[0][fZvtxBin]) fEventList[0][fZvtxBin] = new TList();//0 -> probe
	if(!fEventList[1][fZvtxBin]) fEventList[1][fZvtxBin] = new TList();//1 -> passing probe
	TList *prevEvent_p  = fEventList[0][fZvtxBin];
	TList *prevEvent_pp = fEventList[1][fZvtxBin];

	TrackQA();
  CutEfficiency();

  //Now we either add current events to stack or remove
  //If no electron in current event - no need to add it to mixed
  if(fProbeTrackArray->GetEntriesFast() > 0){
    prevEvent_p ->AddFirst(fProbeTrackArray->Clone());
    prevEvent_pp->AddFirst(fPassingProbeTrackArray->Clone());

    if(prevEvent_p->GetSize() > fNMixed){//Remove redundant events
      TObjArray * tmp_p = static_cast<TObjArray*>(prevEvent_p->Last());
      prevEvent_p->RemoveLast();
      delete tmp_p;
      tmp_p = NULL;

      TObjArray * tmp_pp = static_cast<TObjArray*>(prevEvent_pp->Last());
      prevEvent_pp->RemoveLast();
      delete tmp_pp;
      tmp_pp = NULL;
    }
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //In principle, this function is not needed...

  AliInfo(Form("%s is done.",GetName()));

}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::TrackQA() 
{
  const Int_t trackMult = fEvent->GetNumberOfTracks();

	Double_t vec_3D[3] = {0,0,0};
	UInt_t selectedMask_probe        = (1<<fProbeFilter->GetCuts()->GetEntries())-1;
	UInt_t selectedMask_passingprobe = (1<<fPassingProbeFilter->GetCuts()->GetEntries())-1;

	Double_t values[AliDielectronVarManager::kNMaxValues] = {0.};

	for(Int_t itrack=0;itrack<trackMult;itrack++){
		AliVParticle *particle = (AliVParticle*)fEvent->GetTrack(itrack);

		//reduce unnecessary IsSelected()
		if(particle->Pt() < 0.15) continue;
		if(TMath::Abs(particle->Eta()) > 0.9) continue;

		UInt_t cutmask_probe = fProbeFilter->IsSelected(particle);

		if(cutmask_probe != selectedMask_probe) continue;
		AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(particle);
		AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(particle);

		if(fESDEvent)      AliDielectronVarManager::Fill(esdtrack,values);
		else if(fAODEvent) AliDielectronVarManager::Fill(aodtrack,values);

		vec_3D[0] = values[AliDielectronVarManager::kPt];
		vec_3D[1] = values[AliDielectronVarManager::kEta];
		vec_3D[2] = values[AliDielectronVarManager::kPhi];
		FillSparse(fOutputContainer,"hs_PtEtaPhi",vec_3D);

		FillHistogramTH2(fOutputContainer,"hTrackTPCdEdx",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kTPCsignal]);
		FillHistogramTH2(fOutputContainer,"hTrackITSdEdx",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kITSsignal]);
		FillHistogramTH2(fOutputContainer,"hTrackTOFbeta",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kTOFbeta]);

		FillHistogramTH2(fOutputContainer,"hTrackDCA"            ,values[AliDielectronVarManager::kImpactParXY],values[AliDielectronVarManager::kImpactParZ]);
		FillHistogramTH1(fOutputContainer,"hTrackNclsTPC"        ,values[AliDielectronVarManager::kNclsTPC]);
		FillHistogramTH1(fOutputContainer,"hTrackNclsPIDTPC"     ,values[AliDielectronVarManager::kTPCsignalN]);
		FillHistogramTH1(fOutputContainer,"hTrackNcrTPC"         ,values[AliDielectronVarManager::kNclsCrTPC]);
		FillHistogramTH1(fOutputContainer,"hTrackNfTPC"          ,values[AliDielectronVarManager::kNFclsTPC]);
		FillHistogramTH1(fOutputContainer,"hTrackRatioNcrtoNfTPC",values[AliDielectronVarManager::kNFclsTPCfCross]);
		FillHistogramTH1(fOutputContainer,"hTrackChi2TPC"        ,values[AliDielectronVarManager::kTPCchi2Cl]);
		FillHistogramTH1(fOutputContainer,"hTrackGoldenChi2"     ,values[AliDielectronVarManager::kChi2TPCConstrainedVsGlobal]);
		FillHistogramTH1(fOutputContainer,"hTrackNclsITS"        ,values[AliDielectronVarManager::kNclsITS]);
		FillHistogramTH1(fOutputContainer,"hTrackNscITS"         ,values[AliDielectronVarManager::kNclsSITS]);
		FillHistogramTH1(fOutputContainer,"hTrackChi2ITS"        ,values[AliDielectronVarManager::kITSchi2Cl]);

		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTPCvsPin",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kTPCnSigmaEle]);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElITSvsPin",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kITSnSigmaEle]);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTOFvsPin",values[AliDielectronVarManager::kPIn],values[AliDielectronVarManager::kTOFnSigmaEle]);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTPCvsEta",values[AliDielectronVarManager::kEta],values[AliDielectronVarManager::kTPCnSigmaEle]);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElITSvsEta",values[AliDielectronVarManager::kEta],values[AliDielectronVarManager::kITSnSigmaEle]);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTOFvsEta",values[AliDielectronVarManager::kEta],values[AliDielectronVarManager::kTOFnSigmaEle]);
		fProbeTrackArray->Add(particle);//array of probe electrons for mixed event, // copy constructor?

		UInt_t cutmask_passingprobe = fPassingProbeFilter->IsSelected(particle);
		if(cutmask_passingprobe == selectedMask_passingprobe){
			fPassingProbeTrackArray->Add(particle);//array of passing probe electrons for mixed event
		}
	}//end of track loop

	//create tag track array
	UInt_t selectedMask_tag = (1<<fTagFilter->GetCuts()->GetEntries())-1;
	for(Int_t itrack=0;itrack<trackMult;itrack++){
		AliVParticle *particle = (AliVParticle*)fEvent->GetTrack(itrack);

		//reduce unnecessary IsSelected()
		if(particle->Pt() < 0.15) continue;
		if(TMath::Abs(particle->Eta()) > 0.9) continue;

		UInt_t cutmask_tag= fTagFilter->IsSelected(particle);
		if(cutmask_tag == selectedMask_tag) fTagTrackArray->Add(particle);
	}//end of track loop

}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::CutEfficiency()
{
	//tag and probe method.

	const Double_t Me = TDatabasePDG::Instance()->GetParticle(11)->Mass();

	//const Int_t trackMult = fEvent->GetNumberOfTracks();
	//UInt_t selectedMask_tag          = (1<<fTagFilter->GetCuts()->GetEntries())-1;
	//UInt_t selectedMask_probe        = (1<<fProbeFilter->GetCuts()->GetEntries())-1;
	//UInt_t selectedMask_passingprobe = (1<<fPassingProbeFilter->GetCuts()->GetEntries())-1;

	Double_t value[4] = {0,0,0,0};
	const Int_t trackMult_tag = fTagTrackArray->GetEntries();
	const Int_t trackMult_p   = fProbeTrackArray->GetEntries();
	const Int_t trackMult_pp  = fPassingProbeTrackArray->GetEntries();

	//same particle combination is rejected by pointer.

	//fill ULS -+ same
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() > 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t i2=0;i2<trackMult_p;i2++){
			AliVParticle *par2 = (AliVParticle*)fProbeTrackArray->At(i2);

			if(par2 == par1) continue;
			if(par2->Charge() < 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
			FillHistogramTH2(fOutputContainer,"hPairMvsPhiV",phiv,p12.M());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hProbe_ULS_same",p12.M(),pv2.Pt());

		}//end of par2

		for(Int_t i2=0;i2<trackMult_pp;i2++){
			AliVParticle *par2 = (AliVParticle*)fPassingProbeTrackArray->At(i2);

			if(par2 == par1) continue;
			if(par2->Charge() < 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hPassingProbe_ULS_same",p12.M(),pv2.Pt());
		}//end of par2

	}//end of par1



	//fill ULS +- same
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		//AliVParticle *par1 = (AliVParticle*)fEvent->GetTrack(i1);
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() < 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t i2=0;i2<trackMult_p;i2++){
			//if(i2==i1) continue;//reject same track combination
			//AliVParticle *par2 = (AliVParticle*)fEvent->GetTrack(i2);
			AliVParticle *par2 = (AliVParticle*)fProbeTrackArray->At(i2);
			if(par2 == par1) continue;
			if(par2->Charge() > 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
			FillHistogramTH2(fOutputContainer,"hPairMvsPhiV",phiv,p12.M());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hProbe_ULS_same",p12.M(),pv2.Pt());

		}//end of par2


		for(Int_t i2=0;i2<trackMult_pp;i2++){
			//if(i2==i1) continue;//reject same track combination
			//AliVParticle *par2 = (AliVParticle*)fEvent->GetTrack(i2);
			AliVParticle *par2 = (AliVParticle*)fPassingProbeTrackArray->At(i2);
			if(par2 == par1) continue;
			if(par2->Charge() > 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hPassingProbe_ULS_same",p12.M(),pv2.Pt());
		}//end of par2

	}//end of par1

	//fill LS ++ same
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		//AliVParticle *par1 = (AliVParticle*)fEvent->GetTrack(i1);
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() < 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t i2=0;i2<trackMult_p;i2++){
			//if(i2==i1) continue;//reject same track combination
			//AliVParticle *par2 = (AliVParticle*)fEvent->GetTrack(i2);
			AliVParticle *par2 = (AliVParticle*)fProbeTrackArray->At(i2);

			if(par2 == par1) continue;
			if(par2->Charge() < 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hProbe_LSpp_same",p12.M(),pv2.Pt());

		}//end of par2

		for(Int_t i2=0;i2<trackMult_pp;i2++){
			//if(i2==i1) continue;//reject same track combination
			//AliVParticle *par2 = (AliVParticle*)fEvent->GetTrack(i2);
			AliVParticle *par2 = (AliVParticle*)fPassingProbeTrackArray->At(i2);

			if(par2 == par1) continue;
			if(par2->Charge() < 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hPassingProbe_LSpp_same",p12.M(),pv2.Pt());
		}//end of par2

	}//end of par1


	//fill LS -- same
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() > 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t i2=0;i2<trackMult_p;i2++){
			AliVParticle *par2 = (AliVParticle*)fProbeTrackArray->At(i2);

			if(par2 == par1) continue;
			if(par2->Charge() > 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hProbe_LSnn_same",p12.M(),pv2.Pt());

		}//end of par2

		for(Int_t i2=0;i2<trackMult_pp;i2++){
			AliVParticle *par2 = (AliVParticle*)fPassingProbeTrackArray->At(i2);
			if(par2 == par1) continue;
			if(par2->Charge() > 0) continue;

			TLorentzVector pv2 = TLorentzVector();
			pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

			TLorentzVector p12 = pv1 + pv2;
			value[0] = p12.M();
			value[1] = pv2.Pt();
			value[2] = pv2.Eta();
			value[3] = pv2.Phi();
			if(value[3] < 0) value[3] += TMath::TwoPi();

			Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());

			if(
					(0 < p12.M() && p12.M() < fMmax)
					&& (fPhiVmin < phiv && phiv < TMath::Pi())
				) continue;

			FillHistogramTH2(fOutputContainer,"hPassingProbe_LSnn_same",p12.M(),pv2.Pt());
		}//end of par2

	}//end of par1

	//next mixed event
	TList *prevEvent_p  = fEventList[0][fZvtxBin];
	TList *prevEvent_pp = fEventList[1][fZvtxBin];

	//fill ULS -+ mix
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);
		if(par1->Charge() > 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t iev=0;iev<prevEvent_p->GetEntries();iev++){
			TObjArray *arrmix = static_cast<TObjArray*>(prevEvent_p->At(iev));

			for(Int_t i2=0;i2<arrmix->GetEntriesFast();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() < 0) continue;

				//UInt_t cutmask_probe = fProbeFilter->IsSelected(par2);
				//if(cutmask_probe != selectedMask_probe) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;
				FillHistogramTH2(fOutputContainer,"hProbe_ULS_mix",p12.M(),pv2.Pt());

			}//end of par2
		}//end of mix event loop

		for(Int_t iev=0;iev<prevEvent_pp->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_pp->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() < 0) continue;

				//UInt_t cutmask_probe = fProbeFilter->IsSelected(par2);
				//if(cutmask_probe != selectedMask_probe) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;

				FillHistogramTH2(fOutputContainer,"hPassingProbe_ULS_mix",p12.M(),pv2.Pt());
			}//end of par2
		}//end of mix event loop

	}//end of par1

	//fill ULS +- mix
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() < 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t iev=0;iev<prevEvent_p->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_p->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() > 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;
				FillHistogramTH2(fOutputContainer,"hProbe_ULS_mix",p12.M(),pv2.Pt());
			}//end of par2
		}//end of mix event loop

		for(Int_t iev=0;iev<prevEvent_pp->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_pp->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() > 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;

				FillHistogramTH2(fOutputContainer,"hPassingProbe_ULS_mix",p12.M(),pv2.Pt());
			}//end of par2
		}//end of mix event loop
	}//end of par1

	//fill ULS ++ mix
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);

		if(par1->Charge() < 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t iev=0;iev<prevEvent_p->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_p->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() < 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;
				FillHistogramTH2(fOutputContainer,"hProbe_LSpp_mix",p12.M(),pv2.Pt());

			}//end of par2
		}//end of mix event loop

		for(Int_t iev=0;iev<prevEvent_pp->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_pp->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() < 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;

				FillHistogramTH2(fOutputContainer,"hPassingProbe_LSpp_mix",p12.M(),pv2.Pt());
			}//end of par2
		}//end of mix event loop
	}//end of par1

	//fill ULS -- mix
	for(Int_t i1=0;i1<trackMult_tag;i1++){
		AliVParticle *par1 = (AliVParticle*)fTagTrackArray->At(i1);
		if(par1->Charge() > 0) continue;

		//apply tight cut to particle1
		TLorentzVector pv1 = TLorentzVector();
		pv1.SetPtEtaPhiM(par1->Pt(),par1->Eta(),par1->Phi(),Me);

		for(Int_t iev=0;iev<prevEvent_p->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_p->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() > 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;
				FillHistogramTH2(fOutputContainer,"hProbe_LSnn_mix",p12.M(),pv2.Pt());

			}//end of par2
		}//end of mix event loop

		for(Int_t iev=0;iev<prevEvent_pp->GetEntries();iev++){
			TObjArray *arrmix = (TObjArray*)prevEvent_pp->At(iev);

			for(Int_t i2=0;i2<arrmix->GetEntries();i2++){
				AliVParticle *par2 = (AliVParticle*)arrmix->At(i2);
				if(par2->Charge() > 0) continue;

				TLorentzVector pv2 = TLorentzVector();
				pv2.SetPtEtaPhiM(par2->Pt(),par2->Eta(),par2->Phi(),Me);

				TLorentzVector p12 = pv1 + pv2;
				value[0] = p12.M();
				value[1] = pv2.Pt();
				value[2] = pv2.Eta();
				value[3] = pv2.Phi();
				if(value[3] < 0) value[3] += TMath::TwoPi();
				Float_t phiv = PhivPair(fEvent->GetMagneticField(),par1->Charge(),par2->Charge(),pv1.Vect(),pv2.Vect());
				if(
						(0 < p12.M() && p12.M() < fMmax)
						&& (fPhiVmin < phiv && phiv < TMath::Pi())
					) continue;

				FillHistogramTH2(fOutputContainer,"hPassingProbe_LSnn_mix",p12.M(),pv2.Pt());
			}//end of par2
		}//end of mix event loop
	}//end of par1

}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillHistogramTH1(TList *list, const Char_t *hname, Double_t x, Double_t w, Option_t *opt) const
{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(list->FindObject(hname));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH1) <%s> ",hname));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      // use bin width as weight
      Int_t bin = hist->GetXaxis()->FindBin(x);
      // check if not overflow or underflow bin
      if(bin != 0 && bin != hist->GetXaxis()->GetNbins()){
        Double_t binwidth = hist->GetXaxis()->GetBinWidth(bin);
        myweight = w/binwidth;
      }
    }
    hist->Fill(x, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w, Option_t *opt) const
{
  TH2 * hist = dynamic_cast<TH2*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH2) <%s> ",name));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      Int_t binx = hist->GetXaxis()->FindBin(x);
      if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
    }
    if(optstring.Contains("wy")){
      Int_t biny = hist->GetYaxis()->FindBin(y);
      if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
    }
    hist->Fill(x, y, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w, Option_t *opt) const
{
  TH3 * hist = dynamic_cast<TH3*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance TH3) <%s> ",name));
    return;
  }
  else{
    TString optstring(opt);
    Double_t myweight = optstring.Contains("w") ? 1. : w;
    if(optstring.Contains("wx")){
      Int_t binx = hist->GetXaxis()->FindBin(x);
      if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
    }
    if(optstring.Contains("wy")){
      Int_t biny = hist->GetYaxis()->FindBin(y);
      if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
    }
    if(optstring.Contains("wz")){
      Int_t binz = hist->GetZaxis()->FindBin(z);
      if(binz != 0 && binz != hist->GetZaxis()->GetNbins()) myweight *= 1./hist->GetZaxis()->GetBinWidth(binz);
    }
    hist->Fill(x, y, z, myweight);
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillSparse(TList *list, const Char_t *name, Double_t *x, Double_t w) const
{
  THnSparse * hist = dynamic_cast<THnSparse*>(list->FindObject(name));
  if(!hist){
    AliError(Form("can not find histogram (of instance THnSparse) <%s> ",name));
    return;
  }
  else{
    hist->Fill(x,w);
    return;
  }

}
//_______________________________________________________________________________
Double_t AliAnalysisTaskTagAndProbe::PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2) //const
{
  /// Following the idea to use opening of collinear pairs in magnetic field from e.g. PHENIX
  /// to identify conversions. Angle between ee plane and magnetic field is calculated (0 to pi).
  /// Due to tracking to the primary vertex, conversions with no intrinsic opening angle
  /// always end up as pair in "cowboy" configuration. The function as defined here then
  /// returns values close to pi.
  /// Correlated Like Sign pairs (from double conversion / dalitz + conversion) may show up
  /// at pi or at 0 depending on which leg has the higher momentum. (not checked yet)
  /// This expected ambiguity is not seen due to sorting of track arrays in this framework.
  /// To reach the same result as for ULS (~pi), the legs are flipped for LS.
  /// from PWGDQ/dielectron/core/AliDielectronPair.cxx

  //Define local buffer variables for leg properties
  Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
  Double_t px2=-9999.,py2=-9999.,pz2=-9999.;

  TVector3 fD1=dau1;
  TVector3 fD2=dau2;
  Int_t    d1Q=charge1;
  //Int_t    d2Q=charge2;

  if (charge1*charge2 > 0.) { // Like Sign
    if(MagField<0){ // inverted behaviour
      if(d1Q>0){
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }
    }
  }
  else { // Unlike Sign
    if(MagField>0){ // regular behaviour
      if(d1Q>0){
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();

        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();

        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();

        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();

        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }
    }
  }

  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t pz = pz1+pz2;
  Double_t dppair = TMath::Sqrt(px*px+py*py+pz*pz);

  //unit vector of (pep+pem)
  Double_t pl = dppair;
  Double_t ux = px/pl;
  Double_t uy = py/pl;
  Double_t uz = pz/pl;
  Double_t ax = uy/TMath::Sqrt(ux*ux+uy*uy);
  Double_t ay = -ux/TMath::Sqrt(ux*ux+uy*uy);

  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  //Double_t ptep = iep->Px()*ax + iep->Py()*ay;
  //Double_t ptem = iem->Px()*ax + iem->Py()*ay;

  Double_t pxep = px1;
  Double_t pyep = py1;
  Double_t pzep = pz1;
  Double_t pxem = px2;
  Double_t pyem = py2;
  Double_t pzem = pz2;

  //vector product of pep X pem
  Double_t vpx = pyep*pzem - pzep*pyem;
  Double_t vpy = pzep*pxem - pxep*pzem;
  Double_t vpz = pxep*pyem - pyep*pxem;
  Double_t vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
  //Double_t thev = acos(vpz/vp);

  //unit vector of pep X pem
  Double_t vx = vpx/vp;
  Double_t vy = vpy/vp;
  Double_t vz = vpz/vp;

  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  Double_t wx = uy*vz - uz*vy;
  Double_t wy = uz*vx - ux*vz;
  //Double_t wz = ux*vy - uy*vx;
  //Double_t wl = sqrt(wx*wx+wy*wy+wz*wz);
  // by construction, (wx,wy,wz) must be a unit vector.
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them
  // should be small if the pair is conversion.
  // this function then returns values close to pi!
  Double_t cosPhiV = wx*ax + wy*ay;
  Double_t phiv = TMath::ACos(cosPhiV);

  return phiv;
}
//________________________________________________________________________
