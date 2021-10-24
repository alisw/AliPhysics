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
#include "AliESDv0KineCuts.h"

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
#include "AliAODv0KineCuts.h"

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

#include "AliAnalysisTaskTagAndProbe.h"

//Author: Daiki Sekihata (Center for Nuclear Study, the University of Tokyo)
//daiki.sekihata@cern.ch
//This analysis task is for data-driven efficiency for a single track with tag-and-probe method.
//Since this task aims to measure the efficiency for a single track,
//electrons from a primary vertex or electrons from gamma conversions do not matter.
//In dielectron analyses, a hit on 1st SPD layer is required for electrons.
//Thus, converted electrons are from the beam pipe.
//They are considered as electrons from the primary vertex which have the similar PID/tracking efficiency.

using namespace std;

ClassImp(AliAnalysisTaskTagAndProbe)

//________________________________________________________________________
AliAnalysisTaskTagAndProbe::AliAnalysisTaskTagAndProbe():
  AliAnalysisTaskSE(),
  fOutputContainer(0x0),
	fTriggerMask(0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCEvent(0x0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentrality(0),
  fCentralityMin(0.),
  fCentralityMax(101.),
  fNMixed(10),
  fZvtxBin(-1),
	fPIDResponse(0x0),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	fPIDCalibinPU(kTRUE),
  fPostPIDCntrdCorrTPC(0x0),
  fPostPIDWdthCorrTPC(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0),
	fTagTrackArray(0x0),
	fProbeTrackArray(0x0),
	fPassingProbeTrackArray(0x0),
	fEventFilter(0x0),
	fTagFilter(0x0),
	fProbeFilter(0x0),
	fPassingProbeFilter(0x0),
  fPIDFilter(0x0),
	fMmax(-1),
	fPhiVmin(3.2),
  fESDtrackCutsGlobalNoDCA(0x0),
  fESDv0KineCuts(0x0),
  fAODv0KineCuts(0x0),
  fMCArrayESD(0x0),
  fMCArrayAOD(0x0),
  fPIDCalibMode(kTRUE),
  fIsMC(kFALSE)
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
  fMCEvent(0x0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentrality(0),
  fCentralityMin(0.),
  fCentralityMax(101.),
  fNMixed(10),
  fZvtxBin(-1),
	fPIDResponse(0x0),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	fPIDCalibinPU(kTRUE),
  fPostPIDCntrdCorrTPC(0x0),
  fPostPIDWdthCorrTPC(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0),
	fTagTrackArray(0x0),
	fProbeTrackArray(0x0),
	fPassingProbeTrackArray(0x0),
	fEventFilter(0x0),
	fTagFilter(0x0),
	fProbeFilter(0x0),
	fPassingProbeFilter(0x0),
  fPIDFilter(0x0),
	fMmax(-1),
	fPhiVmin(3.2),
  fESDtrackCutsGlobalNoDCA(0x0),
  fESDv0KineCuts(0x0),
  fAODv0KineCuts(0x0),
  fMCArrayESD(0x0),
  fMCArrayAOD(0x0),
  fPIDCalibMode(kTRUE),
  fIsMC(kFALSE)
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
	fPIDFilter = new AliAnalysisFilter("fPIDFilter","fPIDFilter");

  fESDtrackCutsGlobalNoDCA = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexXY(2.4);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexZ(3.2);
  fESDtrackCutsGlobalNoDCA->SetDCAToVertex2D(kTRUE);

  fESDv0KineCuts = new AliESDv0KineCuts();
  fAODv0KineCuts = new AliAODv0KineCuts();
  fESDv0KineCuts->SetMode(AliESDv0KineCuts::kPurity,AliESDv0KineCuts::kPbPb);
  fAODv0KineCuts->SetMode(AliAODv0KineCuts::kPurity,AliAODv0KineCuts::kPbPb);
  //Float_t Rcut[2] = {0,90};
  //fESDv0KineCuts->SetGammaCutVertexR(Rcut);
  //fAODv0KineCuts->SetGammaCutVertexR(Rcut);
  //fAODv0KineCuts->SetGammaCutInvMass(0.1);//upper limit for mee

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
	delete fPIDFilter;
	delete fTagTrackArray;
	delete fProbeTrackArray;
	delete fPassingProbeTrackArray;
  delete fESDtrackCutsGlobalNoDCA;
  delete fESDv0KineCuts;
  delete fAODv0KineCuts;
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

	const Int_t NdimPU      = 3;
	Int_t NbinPU[NdimPU]    = { 100,  600,  200};
	Double_t xminPU[NdimPU] = {   0, -300,    0};
	Double_t xmaxPU[NdimPU] = {6e+4, +300, 2e+4};

  const TString sidename[3] = {"A","C",""};
  for(Int_t iside=0;iside<3;iside++){
    fOutputContainer->Add(new THnSparseF(Form("hsTPCpileup%s",sidename[iside].Data()),Form("TPC pileup %s;N_{cls}^{SDD+SSD};TPC pileup Z (cm);TPC pileup M;",sidename[iside].Data()),NdimPU,NbinPU,xminPU,xmaxPU));
  }

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

	const TString probetype[2]  = {"Probe","PassingProbe"};
	const TString chargetype[3]  = {"ULS","LSpp","LSnn"};
	const TString eventtype[2] = {"same","mix"};

  for(Int_t ip=0;ip<2;ip++){
    for(Int_t ic=0;ic<3;ic++){
      for(Int_t ie=0;ie<2;ie++){
        TH2F *h2TAP = new TH2F(Form("h%s_%s_%s",probetype[ip].Data(),chargetype[ic].Data(),eventtype[ie].Data()),Form("h%s_%s_%s",probetype[ip].Data(),chargetype[ic].Data(),eventtype[ie].Data()),500,0,5,100,0,10);
        h2TAP->SetXTitle("m_{ee} (GeV/c^{2})");
        h2TAP->SetYTitle("p_{T,e} (GeV/c)");
        h2TAP->Sumw2();
        fOutputContainer->Add(h2TAP);
      }
    }
  }

  fOutputContainer->Add(new TH1F("hV0CosPointingAngle","V0 cos pointing angle;cos(#theta_{point})",100,0,1));
  fOutputContainer->Add(new TH2F("hV0Lxy","V0 L_{xy} vs. m_{ee};V0 L_{xy} (cm);m_{ee} (GeV/c^{2})",600,0,60,50,0,0.05));
  fOutputContainer->Add(new TH2F("hV0Lxy_GammaConv","V0 L_{xy} vs. m_{ee};V0 L_{xy} (cm);m_{ee} (GeV/c^{2})",600,0,60,50,0,0.05));

  fOutputContainer->Add(new TH2F("hV0AP","AP plot",200,-1,+1,300,0,0.3));
  const TString V0name[4] = {"GammaConv","K0S","Lambda","AntiLambda"};
  for(Int_t i=0;i<4;i++) fOutputContainer->Add(new TH2F(Form("hV0AP_%s",V0name[i].Data()),Form("V0 AP plot %s",V0name[i].Data()),200,-1,+1,300,0,0.3));

  const Int_t Ndim_PID = 8;//NclsSDDSSD, puZ, puM, pin, eta, nsigmaTPC, nsigmaITS, nsigmaTOF
  Int_t Nbin_PID[Ndim_PID]    = {    4,      6,    4, 28, 20, 100, 100, 100};
  Double_t xmin_PID[Ndim_PID] = {    0,   -300,    0,  0, -1,  -5,  -5,  -5};
  Double_t xmax_PID[Ndim_PID] = {20000,   +300,10000, 10, +1,  +5,  +5,  +5};

  const TString parname[6] = {"El_online","El_offline","El","Pi","Ka","Pr"};
  const Double_t NSDDSSD[5]       = {0, 2.5e3, 1e+4, 1.6e+4, 1e+5};//clusters on SDD+SSD layers
  const Double_t TPCpileupZ[7]    = {-300,-75,-25,0,+25,+75,+300};//in cm (A+C)/2 average
  const Double_t TPCpileupMult[5] = {0,400,1200,3000,20000};//pileup contributors A+C sum
  const Double_t pinbin[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3,4,5,6,7,8,9,10};//pin for PID calib

  for(Int_t i=0;i<6;i++){
    //if(parname[i] == "El") continue;
    THnSparseF *hsPID = new THnSparseF(Form("hsPID_V0%s",parname[i].Data()),Form("hsPID %s;",parname[i].Data()),Ndim_PID,Nbin_PID,xmin_PID,xmax_PID);
    hsPID->SetBinEdges(0,NSDDSSD);
    hsPID->SetBinEdges(1,TPCpileupZ);
    hsPID->SetBinEdges(2,TPCpileupMult);
    hsPID->SetBinEdges(3,pinbin);
    hsPID->GetAxis(0)->SetTitle("N_{cls}^{SDD+SSD}");
    hsPID->GetAxis(1)->SetTitle("Z_{puv} (cm)");
    hsPID->GetAxis(2)->SetTitle("M_{puv}");
    hsPID->GetAxis(3)->SetTitle("p_{in} (GeV/c)");
    hsPID->GetAxis(4)->SetTitle("#eta");
    hsPID->GetAxis(5)->SetTitle("n #sigma_{TPC}");
    hsPID->GetAxis(6)->SetTitle("n #sigma_{ITS}");
    hsPID->GetAxis(7)->SetTitle("n #sigma_{TOF}");
    fOutputContainer->Add(hsPID);
  }

  THnSparseF *hsAll_El_TAP = new THnSparseF("hsAll_El_TAP","hs all e^{#pm} for PID eff.",Ndim,Nbin,xmin,xmax);
  //hsAll_El_TAP->SetBinEdges(0,pinbin);
  hsAll_El_TAP->GetAxis(0)->SetTitle("p_{T,e} (GeV/c)");
  hsAll_El_TAP->GetAxis(1)->SetTitle("#eta_{e}");
  hsAll_El_TAP->GetAxis(2)->SetTitle("#varphi_{e} (rad.)");
  fOutputContainer->Add(hsAll_El_TAP);

  THnSparseF *hsSel_El_TAP = new THnSparseF("hsSel_El_TAP","hs all e^{#pm} for PID eff.",Ndim,Nbin,xmin,xmax);
  //hsSel_El_TAP->SetBinEdges(0,pinbin);
  hsSel_El_TAP->GetAxis(0)->SetTitle("p_{T,e} (GeV/c)");
  hsSel_El_TAP->GetAxis(1)->SetTitle("#eta_{e}");
  hsSel_El_TAP->GetAxis(2)->SetTitle("#varphi_{e} (rad.)");
  fOutputContainer->Add(hsSel_El_TAP);

  if(fIsMC){
    for(Int_t i=0;i<6;i++){
      if(parname[i].Contains("El_")) continue;
      THnSparseF *hsPID = new THnSparseF(Form("hsPID_MC%s",parname[i].Data()),Form("hsPID %s;",parname[i].Data()),Ndim_PID,Nbin_PID,xmin_PID,xmax_PID);
      hsPID->SetBinEdges(0,NSDDSSD);
      hsPID->SetBinEdges(1,TPCpileupZ);
      hsPID->SetBinEdges(2,TPCpileupMult);
      hsPID->SetBinEdges(3,pinbin);
      hsPID->GetAxis(0)->SetTitle("N_{cls}^{SDD+SSD}");
      hsPID->GetAxis(1)->SetTitle("Z_{puv} (cm)");
      hsPID->GetAxis(2)->SetTitle("M_{puv}");
      hsPID->GetAxis(3)->SetTitle("p_{in} (GeV/c)");
      hsPID->GetAxis(4)->SetTitle("#eta");
      hsPID->GetAxis(5)->SetTitle("n #sigma_{TPC}");
      hsPID->GetAxis(6)->SetTitle("n #sigma_{ITS}");
      hsPID->GetAxis(7)->SetTitle("n #sigma_{TOF}");
      fOutputContainer->Add(hsPID);
    }

    //histograms for PID efficiency in MC
    TH2F *hMCElall = new TH2F("hMCElall","all electron in M.C.;p_{T,e} (GeV/c);#eta_{e}",100,0,10,20,-1,+1);
    hMCElall->Sumw2();
    fOutputContainer->Add(hMCElall);

    TH2F *hMCElselected = new TH2F("hMCElselected","selected electron in M.C.;p_{T,e} (GeV/c);#eta_{e}",100,0,10,20,-1,+1);
    hMCElselected->Sumw2();
    fOutputContainer->Add(hMCElselected);

  }

  AliDielectronPID::SetPIDCalibinPU(fPIDCalibinPU);
  if(fPIDCalibinPU){
    for(Int_t id=0;id<15;id++){//detector loop TPC/ITS/TOF
      for(Int_t ip=0;ip<15;ip++){//particle loop e/mu/pi/k/p
        if(fPostPIDCntrdCorrPU[id][ip]) AliDielectronPID::SetCentroidCorrFunctionPU(id,ip,fPostPIDCntrdCorrPU[id][ip]);
        if(fPostPIDWdthCorrPU[id][ip])  AliDielectronPID::SetWidthCorrFunctionPU(   id,ip,fPostPIDWdthCorrPU[id][ip] );
      }
    }
  }
  else{
    if(fPostPIDCntrdCorrTPC)  AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorrTPC);
    if(fPostPIDWdthCorrTPC)   AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorrTPC);
    if(fPostPIDCntrdCorrITS)  AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
    if(fPostPIDWdthCorrITS)   AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
    if(fPostPIDCntrdCorrTOF)  AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorrTOF);
    if(fPostPIDWdthCorrTOF)   AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDWdthCorrTOF);
  }

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

  AliInputEventHandler *eventHandler   = nullptr;
  AliInputEventHandler *eventHandlerMC = nullptr;

  if((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){//for AOD
    eventHandler   = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
  else if((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliESDInputHandler::Class()){//for ESD
    eventHandlerMC = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    eventHandler   = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  }
  fMCEvent = (AliMCEvent*)eventHandlerMC->MCEvent();

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


	fTagTrackArray->Clear();
	fProbeTrackArray->Clear();
	fPassingProbeTrackArray->Clear();

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
  Int_t NSDDSSD = NclsSDD0 + NclsSDD1 + NclsSSD0 + NclsSSD1;

	Int_t NclsTPC  = 0;
	if(fESDEvent) NclsTPC = fESDEvent->GetNumberOfTPCClusters();
	else if(fAODEvent) NclsTPC = fAODEvent->GetNumberOfTPCClusters();

	FillHistogramTH2(fOutputContainer,"hNTPCclsvsNSDDSSDcls",NclsTPC,NSDDSSD);

  Float_t TPCpileupZA = 0;
  Float_t TPCpileupZC = 0;
  Float_t TPCpileupZ  = 0;
  Float_t TPCpileupMA = 0;
  Float_t TPCpileupMC = 0;
  Float_t TPCpileupM  = 0;

  if(fESDEvent){
    TVectorF tpcVertexInfo(10);
    AliESDUtils::GetTPCPileupVertexInfo(fESDEvent, tpcVertexInfo);
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }
  else if (fAODEvent){
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    static TVectorF dummyVertexInfo(10); // to be used with old AODs w/o vertex info
    const TVectorF &tpcVertexInfo = header->GetTPCPileUpInfo() ? *header->GetTPCPileUpInfo() : dummyVertexInfo;
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }

  Double_t valuePU[3] = {0,0,0};
  valuePU[0] = NSDDSSD;
  //for A side
  valuePU[1] = TPCpileupZA; valuePU[2] = TPCpileupMA;
	FillSparse(fOutputContainer,"hsTPCpileupA",valuePU);
  //for C side
  valuePU[1] = TPCpileupZC; valuePU[2] = TPCpileupMC;
	FillSparse(fOutputContainer,"hsTPCpileupC",valuePU);
  //for AC average
  valuePU[1] = TPCpileupZ; valuePU[2] = TPCpileupM;
	FillSparse(fOutputContainer,"hsTPCpileup",valuePU);

	if(!fEventList[0][fZvtxBin]) fEventList[0][fZvtxBin] = new TList();//0 -> probe
	if(!fEventList[1][fZvtxBin]) fEventList[1][fZvtxBin] = new TList();//1 -> passing probe
	TList *prevEvent_p  = fEventList[0][fZvtxBin];
	TList *prevEvent_pp = fEventList[1][fZvtxBin];

	TrackQA();
  if(fPIDCalibMode){
    if(fESDEvent)      FillV0InfoESD();
    else if(fAODEvent) FillV0InfoAOD();
  }
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

  if(fIsMC){
    if(fESDEvent) GetMCInfoESD();
    else if(fAODEvent) GetMCInfoAOD();
    ProcessMC(option);
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
void AliAnalysisTaskTagAndProbe::ProcessMC(Option_t *option) 
{
	Int_t NclsSDD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(2);
	Int_t NclsSDD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(3);
	Int_t NclsSSD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(4);
	Int_t NclsSSD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(5);
  Int_t NSDDSSD = NclsSDD0 + NclsSDD1 + NclsSSD0 + NclsSSD1;

	Int_t NclsTPC  = 0;
	if(fESDEvent) NclsTPC = fESDEvent->GetNumberOfTPCClusters();
	else if(fAODEvent) NclsTPC = fAODEvent->GetNumberOfTPCClusters();

  Float_t TPCpileupZA = 0;
  Float_t TPCpileupZC = 0;
  Float_t TPCpileupZ  = 0;
  Float_t TPCpileupMA = 0;
  Float_t TPCpileupMC = 0;
  Float_t TPCpileupM  = 0;

  if(fESDEvent){
    TVectorF tpcVertexInfo(10);
    AliESDUtils::GetTPCPileupVertexInfo(fESDEvent, tpcVertexInfo);
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }
  else if (fAODEvent){
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    static TVectorF dummyVertexInfo(10); // to be used with old AODs w/o vertex info
    const TVectorF &tpcVertexInfo = header->GetTPCPileUpInfo() ? *header->GetTPCPileUpInfo() : dummyVertexInfo;
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }

  Double_t value[8] = {};
  for(Int_t i=0;i<8;i++) value[i] = 0.0;
  value[0] = NSDDSSD;
  value[1] = TPCpileupZ;
  value[2] = TPCpileupM;

  const Int_t trackMult = fEvent->GetNumberOfTracks();
	UInt_t selectedMask_probe        = (1<<fProbeFilter->GetCuts()->GetEntries())-1;
	UInt_t selectedMask_passingprobe = (1<<fPassingProbeFilter->GetCuts()->GetEntries())-1;

  for(Int_t itrack=0;itrack<trackMult;itrack++){
    AliVParticle *particle = (AliVParticle*)fEvent->GetTrack(itrack);

    //reduce unnecessary IsSelected()
    if(particle->Pt() < 0.15) continue;
    if(TMath::Abs(particle->Eta()) > 0.9) continue;

    AliVTrack *track = dynamic_cast<AliVTrack*>(particle);
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(label);
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(label, fMCEvent)) continue;//particles from pileup collision should NOT be used.
    if(!p->IsPhysicalPrimary()) continue;
    Int_t pdg = p->GetPdgCode();

    Float_t nsigma_El_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
    Float_t nsigma_El_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
    Float_t nsigma_El_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);

    Float_t nsigma_Pi_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
    Float_t nsigma_Pi_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
    Float_t nsigma_Pi_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);

    Float_t nsigma_Ka_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
    Float_t nsigma_Ka_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
    Float_t nsigma_Ka_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);

    Float_t nsigma_Pr_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
    Float_t nsigma_Pr_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);
    Float_t nsigma_Pr_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);

    value[3] = track->GetTPCmomentum();
    value[4] = track->Eta();
    if(TMath::Abs(pdg) == 11){
      value[5] = nsigma_El_TPC;
      value[6] = nsigma_El_ITS;
      value[7] = nsigma_El_TOF;
      FillSparse(fOutputContainer,"hsPID_MCEl",value);

      UInt_t cutmask_probe = fProbeFilter->IsSelected(particle);
      if(cutmask_probe == selectedMask_probe) FillHistogramTH2(fOutputContainer,"hMCElall",track->Pt(),track->Eta());

      UInt_t cutmask_passingprobe = fPassingProbeFilter->IsSelected(particle);
      if(cutmask_passingprobe == selectedMask_passingprobe) FillHistogramTH2(fOutputContainer,"hMCElselected",track->Pt(),track->Eta());

    }
    else if(TMath::Abs(pdg) == 211){
      value[5] = nsigma_Pi_TPC;
      value[6] = nsigma_Pi_ITS;
      value[7] = nsigma_Pi_TOF;
      FillSparse(fOutputContainer,"hsPID_MCPi",value);
    }
    else if(TMath::Abs(pdg) == 321){
      value[5] = nsigma_Ka_TPC;
      value[6] = nsigma_Ka_ITS;
      value[7] = nsigma_Ka_TOF;
      FillSparse(fOutputContainer,"hsPID_MCKa",value);
    }
    else if(TMath::Abs(pdg) == 2212){
      value[5] = nsigma_Pr_TPC;
      value[6] = nsigma_Pr_ITS;
      value[7] = nsigma_Pr_TOF;
      FillSparse(fOutputContainer,"hsPID_MCPr",value);
    }

  }//end of track loop
}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::TrackQA() 
{
  const Int_t trackMult = fEvent->GetNumberOfTracks();
  AliTOFHeader* tofH = 0x0;          // from v5-02-Rev10 on subtract the start time
  if(fEvent) tofH = (AliTOFHeader*)fEvent->GetTOFHeader();

	Double_t vec_3D[3] = {0,0,0};
	UInt_t selectedMask_probe        = (1<<fProbeFilter->GetCuts()->GetEntries())-1;
	UInt_t selectedMask_passingprobe = (1<<fPassingProbeFilter->GetCuts()->GetEntries())-1;

  Double_t pT=0, eta=0, phi=0, pin=0;
  Double_t TPCsignal=0, ITSsignal = 0, TOFbeta = -1;
  Double_t TPCchi2 = 999, ITSchi2 = 999, ratioCRtoF = 0;
  Int_t NclsTPC = 0, NclsITS = 0, NcrTPC = 0, NscITS = 0, NfTPC = 0, NclsPIDTPC = 0;
  Float_t DCAxy = -999, DCAz = -999;

  Double32_t expt[5] = {0};
  Double_t l = 0, t = 0 , v = 0;

  Float_t nsigma_El_TPC = -999;
  Float_t nsigma_El_ITS = -999;
  Float_t nsigma_El_TOF = -999;

	for(Int_t itrack=0;itrack<trackMult;itrack++){
		AliVParticle *particle = (AliVParticle*)fEvent->GetTrack(itrack);

    if(fIsMC){
      Int_t label = TMath::Abs(particle->GetLabel());
      if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(label, fMCEvent)) continue;//particles from pileup collision should NOT be used.
    }

		//reduce unnecessary IsSelected()
		if(particle->Pt() < 0.15) continue;
		if(TMath::Abs(particle->Eta()) > 0.9) continue;

		UInt_t cutmask_probe = fProbeFilter->IsSelected(particle);
		if(cutmask_probe != selectedMask_probe) continue;

    AliVTrack *track = dynamic_cast<AliVTrack*>(particle);
    pT  = track->Pt();
    eta = track->Eta();
    phi = track->Phi();
    if(phi < 0) phi += TMath::TwoPi();

    vec_3D[0] = pT;
    vec_3D[1] = eta;
    vec_3D[2] = phi;
    FillSparse(fOutputContainer,"hs_PtEtaPhi",vec_3D);
  
    pin = track->GetTPCmomentum();
    TPCsignal = track->GetTPCsignal();
    ITSsignal = track->GetITSsignal();

    ULong64_t status = track->GetStatus();
    Bool_t isTIME = status & AliVTrack::kTIME;
    Bool_t isTOFout = status & AliVTrack::kTOFout;
    Bool_t isTOFOK = isTIME & isTOFout;

    if(isTOFOK){
      for(Int_t i=0;i<5;i++){expt[i] = 0.0;}
      track->GetIntegratedTimes(expt,5);// ps
      l = TMath::C() * expt[0] * 1e-12;  // m
      t = track->GetTOFsignal();      // ps start time subtracted (until v5-02-Rev09)
      if(tofH) t -= fPIDResponse->GetTOFResponse().GetStartTime(track->P()); // ps
      if( (l < 360.e-2 || l > 800.e-2) || (t <= 0.) ) { TOFbeta  = -1; }
      else {
        t *= 1e-12; //ps -> s
        v = l / t;
        TOFbeta = v / TMath::C();
      }
    }
    else TOFbeta = -1;

		FillHistogramTH2(fOutputContainer,"hTrackTPCdEdx",pin,TPCsignal);
		FillHistogramTH2(fOutputContainer,"hTrackITSdEdx",pin,ITSsignal);
		FillHistogramTH2(fOutputContainer,"hTrackTOFbeta",pin,TOFbeta);

    NclsTPC = track->GetNcls(1);
    NclsITS = track->GetNcls(0);
    NcrTPC = track->GetTPCCrossedRows();
    NfTPC = track->GetTPCNclsF();
    NclsPIDTPC = track->GetTPCsignalN();
    ratioCRtoF = NfTPC > 0 ? (Float_t)NcrTPC / (Float_t)NfTPC : 1;

    ITSchi2 = NclsITS > 0 ? track->GetITSchi2() / NclsITS : 999;
    TPCchi2 = NclsTPC > 0 ? track->GetTPCchi2() / NclsTPC : 999;

    DCAxy = -999, DCAz = -999;
    track->GetImpactParameters(DCAxy,DCAz);

    NscITS = 0;
    for(Int_t il=0;il<6;il++){
      if(track->HasSharedPointOnITSLayer(il)) NscITS++;
    }

		FillHistogramTH2(fOutputContainer,"hTrackDCA"            ,DCAxy, DCAz);
		FillHistogramTH1(fOutputContainer,"hTrackNclsTPC"        ,NclsTPC);
		FillHistogramTH1(fOutputContainer,"hTrackNclsPIDTPC"     ,NclsPIDTPC);
		FillHistogramTH1(fOutputContainer,"hTrackNcrTPC"         ,NcrTPC);
		FillHistogramTH1(fOutputContainer,"hTrackNfTPC"          ,NfTPC);
		FillHistogramTH1(fOutputContainer,"hTrackRatioNcrtoNfTPC",ratioCRtoF);
		FillHistogramTH1(fOutputContainer,"hTrackChi2TPC"        ,TPCchi2);
		FillHistogramTH1(fOutputContainer,"hTrackNclsITS"        ,NclsITS);
		FillHistogramTH1(fOutputContainer,"hTrackNscITS"         ,NscITS);
		FillHistogramTH1(fOutputContainer,"hTrackChi2ITS"        ,ITSchi2);

    //printf("pin = %f GeV/c , eta = %f\n",pin, track->Eta());
    nsigma_El_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
    nsigma_El_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
    nsigma_El_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
    //printf("before nsigma_El_TPC = %f\n",nsigma_El_TPC);
    //printf("before nsigma_El_ITS = %f\n",nsigma_El_ITS);
    //printf("before nsigma_El_TOF = %f\n",nsigma_El_TOF);

    nsigma_El_TPC = (fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(track,AliPID::kElectron)) / AliDielectronPID::GetWdthCorr(track,AliPID::kElectron);
    nsigma_El_ITS = (fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(track,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrITS(track,AliPID::kElectron);
    nsigma_El_TOF = (fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrTOF(track,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrTOF(track,AliPID::kElectron);

    //printf("For TPC El : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorr(track,AliPID::kElectron),AliDielectronPID::GetWdthCorr(track,AliPID::kElectron));
    //printf("For ITS El : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorrITS(track,AliPID::kElectron),AliDielectronPID::GetWdthCorrITS(track,AliPID::kElectron));
    //printf("For TOF El : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorrTOF(track,AliPID::kElectron),AliDielectronPID::GetWdthCorrTOF(track,AliPID::kElectron));

    //printf("For TPC Pi : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorr(track,AliPID::kPion),AliDielectronPID::GetWdthCorr(track,AliPID::kPion));
    //printf("For TPC Ka : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorr(track,AliPID::kKaon),AliDielectronPID::GetWdthCorr(track,AliPID::kKaon));
    //printf("For TPC Pr : cntrd = %f , width = %f\n",AliDielectronPID::GetCntrdCorr(track,AliPID::kProton),AliDielectronPID::GetWdthCorr(track,AliPID::kProton));

    //printf("after nsigma_El_TPC = %f\n",nsigma_El_TPC);
    //printf("after nsigma_El_ITS = %f\n",nsigma_El_ITS);
    //printf("after nsigma_El_TOF = %f\n",nsigma_El_TOF);

		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTPCvsPin",pin,nsigma_El_TPC);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElITSvsPin",pin,nsigma_El_ITS);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTOFvsPin",pin,nsigma_El_TOF);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTPCvsEta",eta,nsigma_El_TPC);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElITSvsEta",eta,nsigma_El_ITS);
		FillHistogramTH2(fOutputContainer,"hTrackNsigmaElTOFvsEta",eta,nsigma_El_TOF);

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

  //only for Kaon PID calibration
  //As statistics of K is too large and merge will fail due to too large output size.
  //Thus, Kaon is randomly rejected.
  Float_t nsigma_Ka_TPC = -999;
  Float_t nsigma_Ka_ITS = -999;
  Float_t nsigma_Ka_TOF = -999;
  Double_t value8[8] = {0.};
	Int_t NclsSDD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(2);
	Int_t NclsSDD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(3);
	Int_t NclsSSD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(4);
	Int_t NclsSSD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(5);
  Int_t NSDDSSD = NclsSDD0 + NclsSDD1 + NclsSSD0 + NclsSSD1;

  Float_t TPCpileupZA = 0;
  Float_t TPCpileupZC = 0;
  Float_t TPCpileupZ  = 0;
  Float_t TPCpileupMA = 0;
  Float_t TPCpileupMC = 0;
  Float_t TPCpileupM  = 0;

  if(fESDEvent){
    TVectorF tpcVertexInfo(10);
    AliESDUtils::GetTPCPileupVertexInfo(fESDEvent, tpcVertexInfo);
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }
  else if (fAODEvent){
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    static TVectorF dummyVertexInfo(10); // to be used with old AODs w/o vertex info
    const TVectorF &tpcVertexInfo = header->GetTPCPileUpInfo() ? *header->GetTPCPileUpInfo() : dummyVertexInfo;
    TPCpileupZA = tpcVertexInfo[0];
    TPCpileupZC = tpcVertexInfo[1];
    TPCpileupZ  = tpcVertexInfo[2];
    TPCpileupMA = tpcVertexInfo[3];
    TPCpileupMC = tpcVertexInfo[4];
    TPCpileupM  = tpcVertexInfo[5];
  }

  value8[0] = NSDDSSD;
  value8[1] = TPCpileupZ;
  value8[2] = TPCpileupM;

  TRandom3 *r3 = new TRandom3(1);
	for(Int_t itrack=0;itrack<trackMult;itrack++){
		AliVTrack *track = (AliVTrack*)fEvent->GetTrack(itrack);

		if(track->Pt() < 0.15) continue;
		if(TMath::Abs(track->Eta()) > 0.9) continue;

    if(fESDEvent){
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track);
      if(!fESDtrackCutsGlobalNoDCA->AcceptTrack(esdtrack)) continue;//standard cuts with very loose DCA cut //bit4
    }
    else if(fAODEvent){
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut //bit4
    }

    nsigma_Ka_TPC = (fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(track,AliPID::kKaon)) / AliDielectronPID::GetWdthCorr(track,AliPID::kKaon);
    nsigma_Ka_ITS = (fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon) - AliDielectronPID::GetCntrdCorrITS(track,AliPID::kKaon)) / AliDielectronPID::GetWdthCorrITS(track,AliPID::kKaon);
    nsigma_Ka_TOF = (fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon) - AliDielectronPID::GetCntrdCorrTOF(track,AliPID::kKaon)) / AliDielectronPID::GetWdthCorrTOF(track,AliPID::kKaon);

    if(nsigma_Ka_ITS < -2 || +2 < nsigma_Ka_ITS) continue;
    if(track->GetTPCmomentum() > 0.4 && (nsigma_Ka_TOF < -2 || +2 < nsigma_Ka_TOF)) continue;
    if(r3->Rndm() > 0.001) continue;//accept only 0.1% of Kaons.

    value8[3] = track->GetTPCmomentum();
    value8[4] = track->Eta();
    value8[5] = nsigma_Ka_TPC;
    value8[6] = nsigma_Ka_ITS;
    value8[7] = nsigma_Ka_TOF;
    FillSparse(fOutputContainer,"hsPID_V0Ka",value8);

	}//end of track loop

  delete r3;
  r3 = 0x0;
}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillV0InfoESD()
{
	Int_t NclsSDD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(2);
	Int_t NclsSDD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(3);
	Int_t NclsSSD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(4);
	Int_t NclsSSD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(5);
  Int_t NSDDSSD = NclsSDD0 + NclsSDD1 + NclsSSD0 + NclsSSD1;

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*vVertex);
  Double_t secVtx[3] = {primaryVertexKF.GetX(), primaryVertexKF.GetY(), primaryVertexKF.GetZ()};
  fESDv0KineCuts->SetEvent(InputEvent());
  fESDv0KineCuts->SetPrimaryVertex(&primaryVertexKF);
  const Double_t Me  = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  const Double_t Mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t Mp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

  TVectorF tpcVertexInfo(10);
  AliESDUtils::GetTPCPileupVertexInfo(fESDEvent, tpcVertexInfo);
  Float_t TPCpileupZ  = tpcVertexInfo[2];
  Float_t TPCpileupM  = tpcVertexInfo[5];

	UInt_t selectedMask_pid = (1<<fPIDFilter->GetCuts()->GetEntries())-1;
  Double_t M1 = 0;
  Double_t M2 = 0;
  Double_t M12 = 0;
  Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;
  Float_t qT = 0, alpha = 0;
  Float_t Lxy = 0;

  Float_t nsigma_El_TPC = -999 , nsigma_Pi_TPC = -999 , nsigma_Pr_TPC = -999;
  Float_t nsigma_El_ITS = -999 , nsigma_Pi_ITS = -999 , nsigma_Pr_ITS = -999;
  Float_t nsigma_El_TOF = -999 , nsigma_Pi_TOF = -999 , nsigma_Pr_TOF = -999;

  Double_t value3D[3] = {0,0,0};

  Double_t value[8] = {};
  for(Int_t i=0;i<8;i++) value[i] = 0.0;
  value[0] = NSDDSSD;
  value[1] = TPCpileupZ;
  value[2] = TPCpileupM;
  const Int_t Nv0 = fEvent->GetNumberOfV0s();  

  for(Int_t iv0=0;iv0<Nv0;iv0++){
    AliESDv0 *v0 = (AliESDv0*)fESDEvent->GetV0(iv0);
    if(v0->GetRr() > 60.) continue;
    if(v0->GetDcaV0Daughters() > 0.25) continue;
    if(v0->GetV0CosineOfPointingAngle(secVtx[0],secVtx[1],secVtx[2]) < 0.98) continue;

    AliESDtrack* legPos = fESDEvent->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = fESDEvent->GetTrack(v0->GetNindex());
    if(legPos->Charge() * legNeg->Charge() > 0) continue;//reject same sign pair

    if(legPos->Pt() < 0.15) continue;
    if(legNeg->Pt() < 0.15) continue;
    if(TMath::Abs(legPos->Eta()) > 0.9) continue;
    if(TMath::Abs(legNeg->Eta()) > 0.9) continue;

    Float_t DCAxy_leg = -999, DCAz_leg = -999;
    legPos->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    DCAxy_leg = -999; DCAz_leg = -999;
    legNeg->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    if(legPos->GetKinkIndex(0) != 0) continue;
    if(legNeg->GetKinkIndex(0) != 0) continue;

    if(!(legPos->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(legPos->GetNcls(0) < 2.5) continue;//minimum number of ITS cluster 3
    if(legNeg->GetNcls(0) < 2.5) continue;//minimum number of ITS cluster 3
    if(legPos->GetITSchi2() / legPos->GetNcls(0) > 36.) continue;//maximum chi2 per cluster ITS
    if(legNeg->GetITSchi2() / legNeg->GetNcls(0) > 36.) continue;//maximum chi2 per cluster ITS

    if(!(legPos->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(legPos->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legNeg->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legPos->GetTPCchi2() / legPos->GetNcls(1) > 2.5) continue;//maximum chi2 per cluster TPC
    if(legNeg->GetTPCchi2() / legNeg->GetNcls(1) > 2.5) continue;//maximum chi2 per cluster TPC

    Float_t ratio_pos = legPos->GetTPCNclsF() > 0 ? (Float_t)legPos->GetTPCCrossedRows() / (Float_t)legPos->GetTPCNclsF() : 1.0;
    Float_t ratio_neg = legNeg->GetTPCNclsF() > 0 ? (Float_t)legNeg->GetTPCCrossedRows() / (Float_t)legNeg->GetTPCNclsF() : 1.0;
    if(ratio_pos < 0.8) continue;
    if(ratio_neg < 0.8) continue;

    if(!legPos->HasPointOnITSLayer(0) && !legPos->HasPointOnITSLayer(1)) continue;//accept SPDany
    if(!legNeg->HasPointOnITSLayer(0) && !legNeg->HasPointOnITSLayer(1)) continue;//accept SPDany

    Lxy = v0->GetRr();
    alpha = v0->AlphaV0();
    qT    = v0->PtArmV0();
    FillHistogramTH2(fOutputContainer,"hV0AP",alpha,qT);

    pdgV0 = 0; pdgP = 0; pdgN = 0;
    if(!fESDv0KineCuts->ProcessV0(v0,pdgV0,pdgP,pdgN)) continue;
    for(Int_t i=3;i<8;i++) value[i] = 0.0;

    if(pdgV0 == 22 && TMath::Abs(pdgP) == 11 && TMath::Abs(pdgN) == 11){//GammaConv
      M1 = Me;
      M2 = Me;
      M12 = v0->GetEffMassExplicit(M1,M2);
      FillHistogramTH2(fOutputContainer,"hV0Lxy",Lxy,M12);

      nsigma_El_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kElectron);
      nsigma_El_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kElectron);
      nsigma_El_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kElectron);
      value[3] = legPos->GetTPCmomentum();
      value[4] = legPos->Eta();
      value[5] = nsigma_El_TPC;
      value[6] = nsigma_El_ITS;
      value[7] = nsigma_El_TOF;

      if(HasConversionPointOnSPD(v0,legPos,legNeg)){
        FillHistogramTH2(fOutputContainer,"hV0Lxy_GammaConv",Lxy,M12);
        FillHistogramTH2(fOutputContainer,"hV0AP_GammaConv",alpha,qT);
        FillSparse(fOutputContainer,"hsPID_V0El" ,value);
        if(v0->GetOnFlyStatus()) FillSparse(fOutputContainer,"hsPID_V0El_online" ,value);
        else                     FillSparse(fOutputContainer,"hsPID_V0El_offline",value);
      }

      nsigma_El_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kElectron);
      nsigma_El_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kElectron);
      nsigma_El_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kElectron);
      value[3] = legNeg->GetTPCmomentum();
      value[4] = legNeg->Eta();
      value[5] = nsigma_El_TPC;
      value[6] = nsigma_El_ITS;
      value[7] = nsigma_El_TOF;

      if(HasConversionPointOnSPD(v0,legPos,legNeg)){
        FillSparse(fOutputContainer,"hsPID_V0El" ,value);
        if(v0->GetOnFlyStatus()) FillSparse(fOutputContainer,"hsPID_V0El_online" ,value);
        else                     FillSparse(fOutputContainer,"hsPID_V0El_offline",value);
      }

      if(TMath::Abs(nsigma_El_TPC) < 3.){//electron is pre-selected by loose 3 sigma.
        //for PID efficiency by DDA
        //fill denominator
        value3D[0] = legPos->Pt();
        value3D[1] = legPos->Eta();
        value3D[2] = legPos->Phi();
        FillSparse(fOutputContainer,"hsAll_El_TAP",value3D);
        value3D[0] = legNeg->Pt();
        value3D[1] = legNeg->Eta();
        value3D[2] = legNeg->Phi();
        FillSparse(fOutputContainer,"hsAll_El_TAP",value3D);

        //fill nominator
        UInt_t cutmask_pid = fPIDFilter->IsSelected(legPos);
        if(cutmask_pid == selectedMask_pid){
          value3D[0] = legPos->Pt();
          value3D[1] = legPos->Eta();
          value3D[2] = legPos->Phi();
          FillSparse(fOutputContainer,"hsSel_El_TAP",value3D);
        }
        cutmask_pid = 0;
        cutmask_pid = fPIDFilter->IsSelected(legNeg);
        if(cutmask_pid == selectedMask_pid){
          value3D[0] = legNeg->Pt();
          value3D[1] = legNeg->Eta();
          value3D[2] = legNeg->Phi();
          FillSparse(fOutputContainer,"hsSel_El_TAP",value3D);
        }

      }
    }
    else if(pdgV0 == 310 && TMath::Abs(pdgP) == 211 && TMath::Abs(pdgN) == 211){//K0S
      M1 = Mpi;
      M2 = Mpi;
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_K0S",alpha,qT);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kPion);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion);
        value[6] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion);
        value[7] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion);
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kPion);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion);
        value[6] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion);
        value[7] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion);
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);
      }
    }
    else if(pdgV0 == 3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){//Lambda
      M1 = Mp;
      M2 = Mpi;
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_Lambda",alpha,qT);

        if(pdgP == -211 && pdgN == 2212){//swapped (this does NOT mean AntiLambda)
          M1 = Mpi;
          M2 = Mp;
          legPos = fESDEvent->GetTrack(v0->GetNindex());//proton
          legNeg = fESDEvent->GetTrack(v0->GetPindex());//pi-
        }
        nsigma_Pr_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kProton) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kProton);
        nsigma_Pr_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kProton) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kProton);
        nsigma_Pr_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kProton) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kProton);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = nsigma_Pr_TPC;
        value[6] = nsigma_Pr_ITS;
        value[7] = nsigma_Pr_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pr",value);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kPion);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = nsigma_Pi_TPC;
        value[6] = nsigma_Pi_ITS;
        value[7] = nsigma_Pi_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);
      }
    }
    else if(pdgV0 == -3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){//Anti-Lambda
      M1 = Mpi;
      M2 = Mp;
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_AntiLambda",alpha,qT);

        if(pdgP == -2212 && pdgN == 211){//swapped (this does NOT mean Lambda)
          M1 = Mp;
          M2 = Mpi;
          legPos = fESDEvent->GetTrack(v0->GetNindex());//pi+
          legNeg = fESDEvent->GetTrack(v0->GetPindex());//anti-proton
        }

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kPion);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = nsigma_Pi_TPC;
        value[6] = nsigma_Pi_ITS;
        value[7] = nsigma_Pi_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);

        nsigma_Pr_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kProton) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kProton);
        nsigma_Pr_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kProton) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kProton);
        nsigma_Pr_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kProton) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kProton);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = nsigma_Pr_TPC;
        value[6] = nsigma_Pr_ITS;
        value[7] = nsigma_Pr_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pr",value);
      }
    }

  }

}
//________________________________________________________________________
void AliAnalysisTaskTagAndProbe::FillV0InfoAOD()
{
	Int_t NclsSDD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(2);
	Int_t NclsSDD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(3);
	Int_t NclsSSD0 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(4);
	Int_t NclsSSD1 = fEvent->GetMultiplicity()->GetNumberOfITSClusters(5);
  Int_t NSDDSSD = NclsSDD0 + NclsSDD1 + NclsSSD0 + NclsSSD1;

  AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  static TVectorF dummyVertexInfo(10); // to be used with old AODs w/o vertex info
  const TVectorF &tpcVertexInfo = header->GetTPCPileUpInfo() ? *header->GetTPCPileUpInfo() : dummyVertexInfo;
  Float_t TPCpileupZ  = tpcVertexInfo[2];
  Float_t TPCpileupM  = tpcVertexInfo[5];

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*vVertex);
  fAODv0KineCuts->SetEvent(InputEvent());
  fAODv0KineCuts->SetPrimaryVertex(&primaryVertexKF);

	UInt_t selectedMask_pid = (1<<fPIDFilter->GetCuts()->GetEntries())-1;
  Double_t M12 = 0;
  Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;
  Float_t qT = 0, alpha = 0;
  Float_t Lxy = 0;

  Float_t nsigma_El_TPC = -999 , nsigma_Pi_TPC = -999 , nsigma_Pr_TPC = -999;
  Float_t nsigma_El_ITS = -999 , nsigma_Pi_ITS = -999 , nsigma_Pr_ITS = -999;
  Float_t nsigma_El_TOF = -999 , nsigma_Pi_TOF = -999 , nsigma_Pr_TOF = -999;

  Double_t value3D[3] = {0,0,0};

  const Int_t Nv0 = fEvent->GetNumberOfV0s();  
  Double_t value[8] = {};
  for(Int_t i=0;i<8;i++) value[i] = 0.0;
  value[0] = NSDDSSD;
  value[1] = TPCpileupZ;
  value[2] = TPCpileupM;

  for(Int_t iv0=0;iv0<Nv0;iv0++){
    AliAODv0 *v0 = (AliAODv0*)fAODEvent->GetV0(iv0);
    if(v0->RadiusV0() > 60.) continue;
    if(v0->CosPointingAngle(dynamic_cast<AliAODVertex *>(fAODEvent->GetPrimaryVertex())) < 0.98) continue;

    AliAODTrack *legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
    AliAODTrack *legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));
    if(legPos->Charge() * legNeg->Charge() > 0) continue;//reject same sign pair

    if(legPos->Pt() < 0.15) continue;
    if(legNeg->Pt() < 0.15) continue;
    if(TMath::Abs(legPos->Eta()) > 0.9) continue;
    if(TMath::Abs(legNeg->Eta()) > 0.9) continue;

    if(!legPos->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut //bit4
    if(!legNeg->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut //bit4

    Float_t DCAxy_leg = -999, DCAz_leg = -999;
    legPos->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    DCAxy_leg = -999; DCAz_leg = -999;
    legNeg->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    AliAODVertex *avp = (AliAODVertex*)legPos->GetProdVertex();
    AliAODVertex *avn = (AliAODVertex*)legNeg->GetProdVertex();
    if(avp->GetType() == AliAODVertex::kKink) continue;//reject kink
    if(avn->GetType() == AliAODVertex::kKink) continue;//reject kink

    if(!(legPos->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(legPos->GetNcls(0) < 2.5) continue;//minimum number of ITS cluster 3
    if(legNeg->GetNcls(0) < 2.5) continue;//minimum number of ITS cluster 3
    if(legPos->GetITSchi2() / legPos->GetNcls(0) > 36.) continue;//maximum chi2 per cluster ITS
    if(legNeg->GetITSchi2() / legNeg->GetNcls(0) > 36.) continue;//maximum chi2 per cluster ITS

    if(!(legPos->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(legPos->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legNeg->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legPos->GetTPCchi2() / legPos->GetNcls(1) > 2.5) continue;//maximum chi2 per cluster TPC
    if(legNeg->GetTPCchi2() / legNeg->GetNcls(1) > 2.5) continue;//maximum chi2 per cluster TPC

    Float_t ratio_pos = legPos->GetTPCNclsF() > 0 ? (Float_t)legPos->GetTPCCrossedRows() / (Float_t)legPos->GetTPCNclsF() : 1.0;
    Float_t ratio_neg = legNeg->GetTPCNclsF() > 0 ? (Float_t)legNeg->GetTPCCrossedRows() / (Float_t)legNeg->GetTPCNclsF() : 1.0;
    if(ratio_pos < 0.8) continue;
    if(ratio_neg < 0.8) continue;

    if(!legPos->HasPointOnITSLayer(0) && !legPos->HasPointOnITSLayer(1)) continue;//accept SPDany
    if(!legNeg->HasPointOnITSLayer(0) && !legNeg->HasPointOnITSLayer(1)) continue;//accept SPDany

    Lxy = v0->RadiusV0();//in cm

    alpha = v0->AlphaV0();
    qT    = v0->PtArmV0();
    FillHistogramTH2(fOutputContainer,"hV0AP",alpha,qT);
    FillHistogramTH1(fOutputContainer,"hV0CosPointingAngle",v0->CosPointingAngle(dynamic_cast<AliAODVertex *>(fAODEvent->GetPrimaryVertex())) );

    //ULong64_t status1 = legPos->GetStatus();
    //ULong64_t status2 = legNeg->GetStatus();

    //Bool_t isTIME1 = status1 & AliVTrack::kTIME;
    //Bool_t isTIME2 = status2 & AliVTrack::kTIME;

    //Bool_t isTOFout1 = status1 & AliVTrack::kTOFout;
    //Bool_t isTOFout2 = status2 & AliVTrack::kTOFout;

    //Bool_t isTOFOK1 = isTIME1 & isTOFout1;
    //Bool_t isTOFOK2 = isTIME2 & isTOFout2;

    pdgV0 = 0; pdgP = 0; pdgN = 0;
    if(!fAODv0KineCuts->ProcessV0(v0,pdgV0,pdgP,pdgN)) continue;
    //fV0Mass.push_back(v0->InvMass2Prongs(0,1,TMath::Abs(pdgP),TMath::Abs(pdgN)));

    for(Int_t i=3;i<8;i++) value[i] = 0.0;

    if(pdgV0 == 22 && TMath::Abs(pdgP) == 11 && TMath::Abs(pdgN) == 11){//GammaConv
      M12 = v0->InvMass2Prongs(0,1,TMath::Abs(pdgP),TMath::Abs(pdgN));
      FillHistogramTH2(fOutputContainer,"hV0Lxy",Lxy,M12);

      nsigma_El_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kElectron);
      nsigma_El_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kElectron);
      nsigma_El_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kElectron);
      value[3] = legPos->GetTPCmomentum();
      value[4] = legPos->Eta();
      value[5] = nsigma_El_TPC;
      value[6] = nsigma_El_ITS;
      value[7] = nsigma_El_TOF;

      if(HasConversionPointOnSPD(v0,legPos,legNeg)){
        FillHistogramTH2(fOutputContainer,"hV0Lxy_GammaConv",Lxy,M12);
        FillHistogramTH2(fOutputContainer,"hV0AP_GammaConv",alpha,qT);
        FillSparse(fOutputContainer,"hsPID_V0El" ,value);
        if(v0->GetOnFlyStatus()) FillSparse(fOutputContainer,"hsPID_V0El_online" ,value);
        else                     FillSparse(fOutputContainer,"hsPID_V0El_offline",value);
      }

      nsigma_El_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kElectron);
      nsigma_El_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kElectron);
      nsigma_El_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kElectron)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kElectron);
      value[3] = legNeg->GetTPCmomentum();
      value[4] = legNeg->Eta();
      value[5] = nsigma_El_TPC;
      value[6] = nsigma_El_ITS;
      value[7] = nsigma_El_TOF;

      if(HasConversionPointOnSPD(v0,legPos,legNeg)){
        FillSparse(fOutputContainer,"hsPID_V0El" ,value);
        if(v0->GetOnFlyStatus()) FillSparse(fOutputContainer,"hsPID_V0El_online" ,value);
        else                     FillSparse(fOutputContainer,"hsPID_V0El_offline",value);
      }

      if(TMath::Abs(nsigma_El_TPC) < 3.){//electron is pre-selected by loose 3 sigma.
        //for PID efficiency by DDA
        //fill denominator
        value3D[0] = legPos->Pt();
        value3D[1] = legPos->Eta();
        value3D[2] = legPos->Phi();
        FillSparse(fOutputContainer,"hsAll_El_TAP",value3D);
        value3D[0] = legNeg->Pt();
        value3D[1] = legNeg->Eta();
        value3D[2] = legNeg->Phi();
        FillSparse(fOutputContainer,"hsAll_El_TAP",value3D);

        //fill nominator
        UInt_t cutmask_pid = fPIDFilter->IsSelected(legPos);
        if(cutmask_pid == selectedMask_pid){
          value3D[0] = legPos->Pt();
          value3D[1] = legPos->Eta();
          value3D[2] = legPos->Phi();
          FillSparse(fOutputContainer,"hsSel_El_TAP",value3D);
        }
        cutmask_pid = 0;
        cutmask_pid = fPIDFilter->IsSelected(legNeg);
        if(cutmask_pid == selectedMask_pid){
          value3D[0] = legNeg->Pt();
          value3D[1] = legNeg->Eta();
          value3D[2] = legNeg->Phi();
          FillSparse(fOutputContainer,"hsSel_El_TAP",value3D);
        }
      }
    }
    else if(pdgV0 == 310 && TMath::Abs(pdgP) == 211 && TMath::Abs(pdgN) == 211){//K0S
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_K0S",alpha,qT);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kPion);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion);
        value[6] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion);
        value[7] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion);
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kPion);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion);
        value[6] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion);
        value[7] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion);
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);
      }
    }
    else if(pdgV0 == 3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){//Lambda
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_Lambda",alpha,qT);

        if(pdgP == -211 && pdgN == 2212){//swapped (this does NOT mean AntiLambda)
          legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));//proton
          legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));//pi-
        }
        nsigma_Pr_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kProton) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kProton);
        nsigma_Pr_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kProton) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kProton);
        nsigma_Pr_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kProton) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kProton)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kProton);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = nsigma_Pr_TPC;
        value[6] = nsigma_Pr_ITS;
        value[7] = nsigma_Pr_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pr",value);

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kPion);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = nsigma_Pi_TPC;
        value[6] = nsigma_Pi_ITS;
        value[7] = nsigma_Pi_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);
      }
    }
    else if(pdgV0 == -3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){//Anti-Lambda
      if(!v0->GetOnFlyStatus()){
        FillHistogramTH2(fOutputContainer,"hV0AP_AntiLambda",alpha,qT);

        if(pdgP == -2212 && pdgN == 211){//swapped (this does NOT mean Lambda)
          legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));//pi+
          legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));//anti-proton
        }

        nsigma_Pi_TPC = (fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorr(legPos,AliPID::kPion);
        nsigma_Pi_ITS = (fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrITS(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrITS(legPos,AliPID::kPion);
        nsigma_Pi_TOF = (fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion) - AliDielectronPID::GetCntrdCorrTOF(legPos,AliPID::kPion)) / AliDielectronPID::GetWdthCorrTOF(legPos,AliPID::kPion);
        value[3] = legPos->GetTPCmomentum();
        value[4] = legPos->Eta();
        value[5] = nsigma_Pi_TPC;
        value[6] = nsigma_Pi_ITS;
        value[7] = nsigma_Pi_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pi",value);

        nsigma_Pr_TPC = (fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kProton) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorr(legNeg,AliPID::kProton);
        nsigma_Pr_ITS = (fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kProton) - AliDielectronPID::GetCntrdCorrITS(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorrITS(legNeg,AliPID::kProton);
        nsigma_Pr_TOF = (fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kProton) - AliDielectronPID::GetCntrdCorrTOF(legNeg,AliPID::kProton)) / AliDielectronPID::GetWdthCorrTOF(legNeg,AliPID::kProton);
        value[3] = legNeg->GetTPCmomentum();
        value[4] = legNeg->Eta();
        value[5] = nsigma_Pr_TPC;
        value[6] = nsigma_Pr_ITS;
        value[7] = nsigma_Pr_TOF;
        FillSparse(fOutputContainer,"hsPID_V0Pr",value);
      }
    }

  }//end of v0 loop

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
//_______________________________________________________________________________
void AliAnalysisTaskTagAndProbe::GetMCInfoESD() 
{
  fMCArrayESD = 0x0;
  AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  if(eventHandler){
    AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
    if(mcEventHandler) fMCArrayESD = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
  }

  if(!fMCArrayESD) AliError("Could not get MC Stack!");

}
//_______________________________________________________________________________
void AliAnalysisTaskTagAndProbe::GetMCInfoAOD() 
{
  fMCArrayAOD = 0x0;
  AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(aodHandler){
    AliAODEvent *aod=aodHandler->GetEvent();
    if(aod){
      fMCArrayAOD = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fMCArrayAOD) AliError("Could not retrieve MC array!");
    }
    else AliError("Could not retrieve AOD event!");
  }

}
//_______________________________________________________________________________
