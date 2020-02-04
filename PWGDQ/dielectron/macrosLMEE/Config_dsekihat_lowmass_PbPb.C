#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

void InitHistograms(AliDielectron *die, Int_t cutDefinition);

TString TrackCutnames[] = {
	"DefaultTrackCut_Nsc0"
 ,"DefaultTrackCut_Nsc01"
//,"LooseTrackCut"
//,"TightTrackCut"
//,"PIDCalibTrackCut"
};
const Int_t nTC = sizeof(TrackCutnames)/sizeof(TrackCutnames[0]);
Int_t GetNTC(){return nTC;}

TString PIDnames[] = {
  "DefaultPID"
 ,"ITSTPChadrejORTOFrec"
 ,"ITSTPChadrej"
 ,"ITSTOFrecover"
 ,"TPChadrejORTOFrec"
 ,"TPChadrej"
 ,"TOFrecover"
// ,"noPID"
};
const Int_t nPID = sizeof(PIDnames)/sizeof(PIDnames[0]);
Int_t GetNPID(){return nPID;}

TString PFnames[] = {
  "woPF"
// ,"wPF"
};
const Int_t nPF = sizeof(PFnames)/sizeof(PFnames[0]);
Int_t GetNPF(){return nPF;}

const Int_t nDie = nTC * nPID * nPF;
Int_t GetN(){return nDie;}

AliDielectron* Config_dsekihat_lowmass_PbPb(
		const Int_t cutDefinitionTC,
		const Int_t cutDefinitionPID,
		const Int_t cutDefinitionPF,
    const Float_t PtMin ,
    const Float_t PtMax ,
    const Float_t EtaMin,
    const Float_t EtaMax
		)
{
  // Setup the instance of AliDielectron
  // create the actual framework object
  TString name = Form("%s_%s_%s",TrackCutnames[cutDefinitionTC].Data(),PIDnames[cutDefinitionPID].Data(),PFnames[cutDefinitionPF].Data());

  Bool_t isAOD = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  printf("isAOD = %d\n",isAOD);

  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("cuts: %s",name.Data()));
  die->Print();

  // set track cuts
  LMEECutLib *lib = new LMEECutLib(name);
  die->GetTrackFilter().AddCuts(lib->SetupTrackCuts(PtMin,PtMax,EtaMin,EtaMax));//AliDielectronCutGroup
  die->GetTrackFilter().AddCuts(lib->SetupPIDCuts());//AliDielectronCutGroup
  if(!isAOD) die->GetTrackFilter().AddCuts(lib->SetupESDtrackCuts());

  if(name.Contains("wPF",TString::kIgnoreCase)) die->GetPairPreFilter().AddCuts(lib->SetupPhiVPreFilter());

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);

	if(name.Contains("PIDCalib",TString::kIgnoreCase) || name.Contains("noPID",TString::kIgnoreCase)) die->SetNoPairing(kTRUE);
	else                                                                                              die->SetNoPairing(kFALSE);

  InitHistograms(die,0);

  return die;
}
//______________________________________________________________________________________
//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //Setup histogram classes
  AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(), die->GetTitle());
  TString name = die->GetName();

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  //Event class
  histos->AddClass("Event");

	for (Int_t i=0; i<2; ++i) histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
	//--- Pair sub-classes
	//    same event ++, +-, --
	for (Int_t i=0; i<3; ++i) histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));

  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","hVertexZ","vertex Z;Z_{vtx} (cm)",40,-20,20,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","hCentralityV0M","centrality;centrality V0M (%)",101,0,101,AliDielectronVarManager::kCentralityNew);//V0M in Run2
  histos->UserHistogram("Event","hCentralityCL0","centrality;centrality CL0 (%)",101,0,101,AliDielectronVarManager::kCentralityCL0);
  histos->UserHistogram("Event","hCentralityCL1","centrality;centrality CL1 (%)",101,0,101,AliDielectronVarManager::kCentralityCL1);
  //histos->UserHistogram("Event","hCentralityV0A","centrality;centrality V0A (%)",101,0,101,AliDielectronVarManager::kCentralityV0A);
  //histos->UserHistogram("Event","hCentralityV0C","centrality;centrality V0C (%)",101,0,101,AliDielectronVarManager::kCentralityV0C);
  //histos->UserHistogram("Event","hCentralityZNA","centrality;centrality ZNA (%)",101,0,101,AliDielectronVarManager::kCentralityZNA);
  //histos->UserHistogram("Event","hCentralityZNC","centrality;centrality ZNC (%)",101,0,101,AliDielectronVarManager::kCentralityZNC);
  histos->UserHistogram("Event","hNclsITS1","Number of clusters on ITS1;N_{cls}^{ITS1}",200,0,2e+4,AliDielectronVarManager::kNclsITS1);
  histos->UserHistogram("Event","hNclsITS2","Number of clusters on ITS2;N_{cls}^{ITS2}",200,0,2e+4,AliDielectronVarManager::kNclsITS2);
  histos->UserHistogram("Event","hNclsITS3","Number of clusters on ITS3;N_{cls}^{ITS3}",200,0,2e+4,AliDielectronVarManager::kNclsITS3);
  histos->UserHistogram("Event","hNclsITS4","Number of clusters on ITS4;N_{cls}^{ITS4}",200,0,2e+4,AliDielectronVarManager::kNclsITS4);
  histos->UserHistogram("Event","hNclsITS5","Number of clusters on ITS5;N_{cls}^{ITS5}",200,0,2e+4,AliDielectronVarManager::kNclsITS5);
  histos->UserHistogram("Event","hNclsITS6","Number of clusters on ITS6;N_{cls}^{ITS6}",200,0,2e+4,AliDielectronVarManager::kNclsITS6);
  histos->UserHistogram("Event","hTPCpileupA","TPC pileup Z vs. M on A side;Z_{puv}^{A} (cm);N_{contrib}^{A}",500,-250,250,500,0,5000,AliDielectronVarManager::kTPCpileupZA,AliDielectronVarManager::kTPCpileupMA);
  histos->UserHistogram("Event","hTPCpileupC","TPC pileup Z vs. M on C side;Z_{puv}^{C} (cm);N_{contrib}^{C}",500,-250,250,500,0,5000,AliDielectronVarManager::kTPCpileupZC,AliDielectronVarManager::kTPCpileupMC);
  histos->UserHistogram("Event","hTPCpileup","TPC pileup Z vs. M;(Z_{puv}^{A} + Z_{puv}^{C})/2 (cm);N_{contrib}^{A}+N_{contrib}^{C}",500,-250,250,100,0,10000,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
  histos->UserHistogram("Event","hV0TPCcorr","Number of clusters of TPC vs. V0 multiplicity;N_{cls}^{TPC};total V0 multiplicity",500,0,5e+6,600,0,6e+4,AliDielectronVarManager::kNTPCclsEvent,AliDielectronVarManager::kMultV0);
  histos->UserHistogram("Event","hV0MvsNcontrib","V0 multiplicity vs. N_{contributors}^{PV};total V0 multiplicity;N_{contributors} to PV",120,0,6e+4,500,0,5000,AliDielectronVarManager::kMultV0,AliDielectronVarManager::kNVtxContrib);

  //add histograms to track class
	const Int_t Ndim = 3;
	Int_t Nbin[Ndim]    = {100, 20,             36};
	Double_t xmin[Ndim] = {  0, -1,              0};
	Double_t xmax[Ndim] = { 10, +1, TMath::TwoPi()};
	UInt_t var[Ndim]  = {AliDielectronVarManager::kPt, AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi};
  histos->UserSparse("Track","hsPtEtaPhi","p_{T} vs. #eta vs. #varphi;",Ndim,Nbin,xmin,xmax,var);
  histos->UserHistogram("Track","hDCAxyz","DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)",100,-5.,5,100,-5.,5,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","hNclsTPC","Number of clusters TPC; N_{cls}^{TPC};Number of tracks",161,-0.5,160.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","hNcrTPC","Number of crossed rows TPC;N_{crossed rows}^{TPC};Number of track",161,-0.5,160.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","hRatioNcrToNf","N_{cr}/N_{f} TPC;N_{cr}/N_{f} TPC;Number of tracks",200,0.,2.0,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","hChi2TPC","chi2/N_{cls}^{TPC};#chi^{2}/N_{cls}^{TPC};Number of tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","hNclsITS","Number of clusters ITS; N_{cls}^{ITS};Number of tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","hChi2ITS","chi2/N_{cls}^{ITS};#chi^{2}/N_{cls}^{ITS};Number of tracks",100,0.,10,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","hNscITS","Number of shared clusters ITS;N_{sc}^{ITS};Number of tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsSITS);

  histos->UserHistogram("Track","hTPCdEdxvsPin","TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx (a.u.)",200,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","hITSdEdxvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);ITS dE/dx (a.u.)",200,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","hTOFbetavsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);TOF #beta"       ,200,0.,10.,240,0.,1.2 ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  //histos->UserHistogram("Track","hTPCnSigmaElvsPin","TPC dE/dx vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TPC}",200,0.,10.,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  //histos->UserHistogram("Track","hITSnSigmaElvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{ITS}",200,0.,10.,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  //histos->UserHistogram("Track","hTOFnSigmaElvsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TOF}",200,0.,10.,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

//  if(name.Contains("PIDCalib",TString::kIgnoreCase)){
    const Int_t Ndim_PID = 8;
    Int_t Nbin_PID[Ndim_PID]    = {   20,     20,  60,100, 20, 100, 100, 100};
    Double_t xmin_PID[Ndim_PID] = {    0,   -250,   0,  0, -1,  -5,  -5,  -5};
    Double_t xmax_PID[Ndim_PID] = {20000,   +250,6000, 10, +1,  +5,  +5,  +5};
    UInt_t var_PID[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kITSnSigmaEle, AliDielectronVarManager::kTOFnSigmaEle};
    histos->UserSparse("Track","hsPID_V0El","hsPID_V0El",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID);

    //UInt_t var_PID_Pio[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPio, AliDielectronVarManager::kITSnSigmaPio, AliDielectronVarManager::kTOFnSigmaPio};
    //histos->UserSparse("Track","hsPID_V0Pi","hsPID_V0Pi",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pio);

    //UInt_t var_PID_Pro[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPro, AliDielectronVarManager::kITSnSigmaPro, AliDielectronVarManager::kTOFnSigmaPro};
    //histos->UserSparse("Track","hsPID_V0Pr","hsPID_V0Pr",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pro);

    //histos->UserHistogram("Track","hAPplot","AP plot;#alpha;q_{T} (GeV/c)",100,-1.,+1,300,0.,0.3,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);
//  }

  const Int_t Nmee = 150;
  Double_t mee[Nmee] = {};
  for(Int_t i=0  ;i<110 ;i++) mee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 1.1 GeV/c2, every 0.01 GeV/c2
  for(Int_t i=110;i<Nmee;i++) mee[i] = 0.1  * (i-110) +  1.1;//from 1.1 to 5 GeV/c2, evety 0.1 GeV/c2

  const Int_t NpTee = 111;
  Double_t pTee[NpTee] = {};
  for(Int_t i=0  ;i<10   ;i++) pTee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 0.09 GeV/c, every 0.01 GeV/c
  for(Int_t i=10 ;i<100  ;i++) pTee[i] = 0.1  * (i- 10) +  0.1;//from 0.1 to 10 GeV/c, evety 0.1 GeV/c
  for(Int_t i=100;i<NpTee;i++) pTee[i] = 1.0  * (i-100) + 10.0;//from 10 to 20 GeV/c, evety 1.0 GeV/c

  //add histograms to pair classes
  //histos->UserHistogram("Pair","hMvsPtvsPhiV"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 500,0.,5.,100,0.,20,90,0.,TMath::Pi(), AliDielectronVarManager::kM, AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);

  //histos->UserHistogram("Pair","hMvsPtvsPhiV"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)",v_mee, v_pTee,v_phiv, AliDielectronVarManager::kM, AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);
  //histos->UserHistogram("Pair","hEtaPhi" ,"pair #eta vs #varphi;#varphi (rad.);#eta",72,0.,TMath::TwoPi(), 200,-1.,1.,AliDielectronVarManager::kPhi, AliDielectronVarManager::kEta);

  const Int_t Ndim_pair = 5;//mee, pT,eta,phi,phiV
  Int_t Nbin_pair[Ndim_pair]    = {149, 110, 20,          72,          72};
  Double_t xmin_pair[Ndim_pair] = {  0,   0, -1,           0,           0};
  Double_t xmax_pair[Ndim_pair] = {  5,  20, +1, TMath::Pi(), TMath::Pi()};
  UInt_t var_pair[Ndim_pair] = {AliDielectronVarManager::kM,AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi, AliDielectronVarManager::kPhivPair};

  histos->UserSparse("Pair","hsPair","ee pairs",Ndim_pair,Nbin_pair,xmin_pair,xmax_pair,var_pair);

  for (Int_t i=0; i<8; ++i){
    THnSparseD *hs_pair = dynamic_cast<THnSparseD*>(histos->GetHist(Form("Pair_%s",AliDielectron::PairClassName(i)),"hsPair"));
    if(!hs_pair) continue;
    hs_pair->SetBinEdges(0,mee);
    hs_pair->SetBinEdges(1,pTee);
    hs_pair->GetAxis(0)->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    hs_pair->GetAxis(1)->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    hs_pair->GetAxis(2)->SetTitle("#it{#eta}_{ee}");
    hs_pair->GetAxis(3)->SetTitle("#it{#varphi}_{ee} (rad.)");
    hs_pair->GetAxis(4)->SetTitle("#it{#varphi}_{V} (rad.)");
  }

  die->SetHistogramManager(histos);

}
//___________________________________________________________________________________________________

