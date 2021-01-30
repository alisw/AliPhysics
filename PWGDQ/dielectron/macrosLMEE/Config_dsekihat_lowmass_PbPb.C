#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

void InitHistograms(AliDielectron *die, Int_t cutDefinition,Bool_t isMC);
void SetupMCSignals(AliDielectron *die);

TString TrackCutnames[] = {
//	"DefaultTrackCut_Nsc0"
// "DefaultTrackCut_Nsc01",
//  "DefaultTrackCut_Nsc01"
 "DefaultTrackCut_Nsc01_woPU"
// ,"DefaultTrackCut_Nsc01_onlyPU"
//,"LooseTrackCut"
//,"TightTrackCut"
//,"PIDCalibTrackCut"
};

const Int_t nTC = sizeof(TrackCutnames)/sizeof(TrackCutnames[0]);
Int_t GetNTC(){return nTC;}

TString PIDnames[] = {
  "DefaultPID",
//  ,"ITSTPChadrejORTOFrec"
// ,"ITSTPChadrej"
// ,"ITSTOFrecover"
 "TPChadrejORTOFrec"
 ,"TPChadrejORTOFrec_daiki"
// ,"TPChadrej"
// ,"TOFrecover"
// ,"TightTPCTOF"
// ,"TightTPC"
 ,"noPID"
};
const Int_t nPID = sizeof(PIDnames)/sizeof(PIDnames[0]);
Int_t GetNPID(){return nPID;}

TString PFnames[] = {
  "woPF"
// ,"wPF"
};
const Int_t nPF = sizeof(PFnames)/sizeof(PFnames[0]);
//Int_t GetNPF(){return nPF;}
Int_t GetNPF(){return 1;}

const Int_t nDie = nTC * nPID * nPF;
Int_t GetN(){return nDie;}

AliDielectron* Config_dsekihat_lowmass_PbPb(
		const Int_t cutDefinitionTC,
		const Int_t cutDefinitionPID,
		const Int_t cutDefinitionPF,
		const Bool_t applyPairCut,
		const Bool_t isMC
		)
{
	// Setup the instance of AliDielectron
	// create the actual framework object
	TString name = Form("%s_%s_PairCut%d",TrackCutnames[cutDefinitionTC].Data(),PIDnames[cutDefinitionPID].Data(),applyPairCut);

	Bool_t isAOD = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
	printf("isAOD = %d\n",isAOD);

	AliDielectron *die = new AliDielectron(Form("anaFilter_%s",name.Data()), Form("cuts: %s",name.Data()));
	die->Print();

	// set track cuts
	LMEECutLib *lib = new LMEECutLib(name);
	die->GetTrackFilter().AddCuts(lib->SetupTrackCuts());//AliDielectronCutGroup
	die->GetTrackFilter().AddCuts(lib->SetupPIDCuts());//AliDielectronCutGroup
	if(!isAOD) die->GetTrackFilter().AddCuts(lib->SetupESDtrackCuts());

	//if(name.Contains("wPF",TString::kIgnoreCase)) die->GetPairPreFilter().AddCuts(lib->SetupPhiVPreFilter());
	if(applyPairCut) die->GetPairFilter().AddCuts(lib->SetupPairCuts());

	//pairing with TLorentzVector
	die->SetUseKF(kFALSE);
  if(!isMC){
    if(name.Contains("woPU")){
      TF1 *f1min = new TF1("f1min","pol2(0)",0,1e+8);
      f1min->SetNpx(1000);
      f1min->FixParameter(0,-3000);
      f1min->FixParameter(1,0.0099);
      f1min->FixParameter(2,9.42e-10);
      AliDielectronEventCuts*  pileupcuts = new AliDielectronEventCuts("pileupcuts","pileupcuts");
      pileupcuts->SetMinCorrCutFunction(f1min, AliDielectronVarManager::kNTPCclsEvent, AliDielectronVarManager::kNSDDSSDclsEvent);
      pileupcuts->Print();
      die->GetEventFilter().AddCuts(pileupcuts);
    }
    else if(name.Contains("onlyPU")){
      TF1 *f1min = new TF1("f1min","pol2(0)",0,1e+8);
      f1min->SetNpx(1000);
      f1min->FixParameter(0,-3000);
      f1min->FixParameter(1,0.0099);
      f1min->FixParameter(2,9.42e-10);
      AliDielectronEventCuts*  pileupcuts = new AliDielectronEventCuts("pileupcuts","pileupcuts");
      pileupcuts->SetMaxCorrCutFunction(f1min, AliDielectronVarManager::kNTPCclsEvent, AliDielectronVarManager::kNSDDSSDclsEvent);
      pileupcuts->Print();
      die->GetEventFilter().AddCuts(pileupcuts);
    }

  }

	if(name.Contains("PIDCalib",TString::kIgnoreCase)) die->SetNoPairing(kTRUE);
	else                                               die->SetNoPairing(kFALSE);

	if(isMC) SetupMCSignals(die);
	InitHistograms(die,0,isMC);

	return die;
}
//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition, Bool_t isMC)
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

	if(die->GetMCSignals()) {
		for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
			histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
			histos->AddClass(Form("Track_%s_%s",AliDielectron::PairClassName(1),die->GetMCSignals()->At(i)->GetName()));
			//histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
			//histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
		}
	}


	//add histograms to event class
	histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
	histos->UserHistogram("Event","hVertexZ","vertex Z;Z_{vtx} (cm)",400,-20,20,AliDielectronVarManager::kZvPrim);
	histos->UserHistogram("Event","hCentralityV0M","centrality;centrality V0M (%)",101,0,101,AliDielectronVarManager::kCentralityNew);//V0M in Run2
	histos->UserHistogram("Event","hCentralityCL0","centrality;centrality CL0 (%)",101,0,101,AliDielectronVarManager::kCentralityCL0);
	histos->UserHistogram("Event","hCentralityCL1","centrality;centrality CL1 (%)",101,0,101,AliDielectronVarManager::kCentralityCL1);
	//histos->UserHistogram("Event","hCentralityV0A","centrality;centrality V0A (%)",101,0,101,AliDielectronVarManager::kCentralityV0A);
	//histos->UserHistogram("Event","hCentralityV0C","centrality;centrality V0C (%)",101,0,101,AliDielectronVarManager::kCentralityV0C);
	//histos->UserHistogram("Event","hCentralityZNA","centrality;centrality ZNA (%)",101,0,101,AliDielectronVarManager::kCentralityZNA);
	//histos->UserHistogram("Event","hCentralityZNC","centrality;centrality ZNC (%)",101,0,101,AliDielectronVarManager::kCentralityZNC);
	//histos->UserHistogram("Event","hNclsITS1","Number of clusters on ITS1;N_{cls}^{ITS1}",200,0,2e+4,AliDielectronVarManager::kNclsITS1);
	//histos->UserHistogram("Event","hNclsITS2","Number of clusters on ITS2;N_{cls}^{ITS2}",200,0,2e+4,AliDielectronVarManager::kNclsITS2);
	//histos->UserHistogram("Event","hNclsITS3","Number of clusters on ITS3;N_{cls}^{ITS3}",200,0,2e+4,AliDielectronVarManager::kNclsITS3);
	//histos->UserHistogram("Event","hNclsITS4","Number of clusters on ITS4;N_{cls}^{ITS4}",200,0,2e+4,AliDielectronVarManager::kNclsITS4);
	//histos->UserHistogram("Event","hNclsITS5","Number of clusters on ITS5;N_{cls}^{ITS5}",200,0,2e+4,AliDielectronVarManager::kNclsITS5);
	//histos->UserHistogram("Event","hNclsITS6","Number of clusters on ITS6;N_{cls}^{ITS6}",200,0,2e+4,AliDielectronVarManager::kNclsITS6);
	histos->UserHistogram("Event","hTPCpileupA","TPC pileup Z vs. M on A side;Z_{puv}^{A} (cm);N_{contrib}^{A}",500,-250,250,500,0,5000,AliDielectronVarManager::kTPCpileupZA,AliDielectronVarManager::kTPCpileupMA);
	histos->UserHistogram("Event","hTPCpileupC","TPC pileup Z vs. M on C side;Z_{puv}^{C} (cm);N_{contrib}^{C}",500,-250,250,500,0,5000,AliDielectronVarManager::kTPCpileupZC,AliDielectronVarManager::kTPCpileupMC);
	histos->UserHistogram("Event","hTPCpileup","TPC pileup Z vs. M;(Z_{puv}^{A} + Z_{puv}^{C})/2 (cm);N_{contrib}^{A}+N_{contrib}^{C}",500,-250,250,100,0,10000,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
	histos->UserHistogram("Event","hV0TPCcorr","Number of clusters of TPC vs. V0 multiplicity;N_{cls}^{TPC};total V0 multiplicity",500,0,5e+6,600,0,6e+4,AliDielectronVarManager::kNTPCclsEvent,AliDielectronVarManager::kMultV0);
	histos->UserHistogram("Event","hV0MvsNcontrib","V0 multiplicity vs. N_{contributors}^{PV};total V0 multiplicity;N_{contributors} to PV",120,0,6e+4,500,0,5000,AliDielectronVarManager::kMultV0,AliDielectronVarManager::kNVtxContrib);
	histos->UserHistogram("Event","hNTPCclsvsNSDDSSDcls","Number of clusters of TPC vs. SDD+SSD cls;N_{cls}^{TPC};N_{cls}^{SDD+SSD}",300,0,6e+6,300,0,6e+4,AliDielectronVarManager::kNTPCclsEvent,AliDielectronVarManager::kNSDDSSDclsEvent);

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
	histos->UserHistogram("Track","hITSdEdxvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);ITS dE/dx (a.u.)",200,0.,10.,400,0.,400.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
	histos->UserHistogram("Track","hTOFbetavsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);TOF #beta"       ,200,0.,10.,240,0.,1.2 ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);

	histos->UserHistogram("Track","hTPCnSigmaElvsPin","TPC n#sigma_{e}^{TPC} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TPC}",200,0,+10.0,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
	histos->UserHistogram("Track","hITSnSigmaElvsPin","ITS n#sigma_{e}^{ITS} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{ITS}",200,0,+10.0,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
	histos->UserHistogram("Track","hTOFnSigmaElvsPin","TOF n#sigma_{e}^{TOF} vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TOF}",200,0,+10.0,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
	histos->UserHistogram("Track","hTPCnSigmaElvsEta","TPC n#sigma_{e}^{TPC} vs. #eta;#eta;n#sigma_{e}^{TPC}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
	histos->UserHistogram("Track","hITSnSigmaElvsEta","ITS n#sigma_{e}^{ITS} vs. #eta;#eta;n#sigma_{e}^{ITS}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
	histos->UserHistogram("Track","hTOFnSigmaElvsEta","TOF n#sigma_{e}^{TOF} vs. #eta;#eta;n#sigma_{e}^{TOF}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);

	histos->UserHistogram("Track","hTPCnSigmaPivsPin","TPC n#sigma_{#pi}^{TPC} vs. p_{in};p_{in} (GeV/c);n#sigma_{#pi}^{TPC}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
	histos->UserHistogram("Track","hTPCnSigmaPivsEta","TPC n#sigma_{#pi}^{TPC} vs. #eta;#eta;n#sigma_{#pi}^{TPC}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio);
	//  //histos->UserHistogram("Track","hITSnSigmaPivsEta","ITS n#sigma_{#pi}^{ITS} vs. #eta;#eta;n#sigma_{#pi}^{ITS}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaPio);
	//  //histos->UserHistogram("Track","hTOFnSigmaPivsEta","TOF n#sigma_{#pi}^{TOF} vs. #eta;#eta;n#sigma_{#pi}^{TOF}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPio);
	//  //histos->UserHistogram("Track","hITSnSigmaPivsPin","ITS n#sigma_{#pi}^{ITS} vs. p_{in};p_{in} (GeV/c);n#sigma_{#pi}^{ITS}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaPio);
	//  //histos->UserHistogram("Track","hTOFnSigmaPivsPin","TOF n#sigma_{#pi}^{TOF} vs. p_{in};p_{in} (GeV/c);n#sigma_{#pi}^{TOF}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);

	histos->UserHistogram("Track","hTPCnSigmaKavsPin","TPC n#sigma_{K}^{TPC} vs. p_{in};p_{in} (GeV/c);n#sigma_{K}^{TPC}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
	histos->UserHistogram("Track","hTPCnSigmaKavsEta","TPC n#sigma_{K}^{TPC} vs. #eta;#eta;n#sigma_{K}^{TPC}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
	//  //histos->UserHistogram("Track","hITSnSigmaKavsEta","ITS n#sigma_{K}^{ITS} vs. #eta;#eta;n#sigma_{K}^{ITS}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaKao);
	//  //histos->UserHistogram("Track","hTOFnSigmaKavsEta","TOF n#sigma_{K}^{TOF} vs. #eta;#eta;n#sigma_{K}^{TOF}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaKao);
	//  //histos->UserHistogram("Track","hITSnSigmaKavsPin","ITS n#sigma_{K}^{ITS} vs. p_{in};p_{in} (GeV/c);n#sigma_{K}^{ITS}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaKao);
	//  //histos->UserHistogram("Track","hTOFnSigmaKavsPin","TOF n#sigma_{K}^{TOF} vs. p_{in};p_{in} (GeV/c);n#sigma_{K}^{TOF}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);

	histos->UserHistogram("Track","hTPCnSigmaPrvsPin","TPC n#sigma_{p}^{TPC} vs. p_{in};p_{in} (GeV/c);n#sigma_{p}^{TPC}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
	histos->UserHistogram("Track","hTPCnSigmaPrvsEta","TPC n#sigma_{p}^{TPC} vs. #eta;#eta;n#sigma_{p}^{TPC}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPro);
	//  //histos->UserHistogram("Track","hITSnSigmaPrvsEta","ITS n#sigma_{p}^{ITS} vs. #eta;#eta;n#sigma_{p}^{ITS}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaPro);
	//  //histos->UserHistogram("Track","hTOFnSigmaPrvsEta","TOF n#sigma_{p}^{TOF} vs. #eta;#eta;n#sigma_{p}^{TOF}",20,-1.0,+1.0,100,-5,+5,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPro);
	//  //histos->UserHistogram("Track","hITSnSigmaPrvsPin","ITS n#sigma_{p}^{ITS} vs. p_{in};p_{in} (GeV/c);n#sigma_{p}^{ITS}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaPro);
	//  //histos->UserHistogram("Track","hTOFnSigmaPrvsPin","TOF n#sigma_{p}^{TOF} vs. p_{in};p_{in} (GeV/c);n#sigma_{p}^{TOF}",200,0,10,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);


	if(name.Contains("PIDCalib",TString::kIgnoreCase) || isMC){
		const Int_t Ndim_PID = 8;
		Int_t Nbin_PID[Ndim_PID]    = {    4,      8,    4,100, 20, 100, 100, 100};
		Double_t xmin_PID[Ndim_PID] = {    0,   -250,    0,  0, -1,  -5,  -5,  -5};
		Double_t xmax_PID[Ndim_PID] = {20000,   +250,10000, 10, +1,  +5,  +5,  +5};
		UInt_t var_PID[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kITSnSigmaEle, AliDielectronVarManager::kTOFnSigmaEle};
		histos->UserSparse("Track","hsPID_V0El","hsPID_V0El",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID);

		//UInt_t var_PID_Pio[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPio, AliDielectronVarManager::kITSnSigmaPio, AliDielectronVarManager::kTOFnSigmaPio};
		//histos->UserSparse("Track","hsPID_V0Pi","hsPID_V0Pi",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pio);

		//UInt_t var_PID_Pro[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPro, AliDielectronVarManager::kITSnSigmaPro, AliDielectronVarManager::kTOFnSigmaPro};
		//histos->UserSparse("Track","hsPID_V0Pr","hsPID_V0Pr",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pro);

		//histos->UserHistogram("Track","hAPplot","AP plot;#alpha;q_{T} (GeV/c)",100,-1.,+1,300,0.,0.3,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);

		const Double_t NSPD0[5]         = {0,1000,3000,5000,20000};//clusters on a first SPD layer
		const Double_t TPCpileupZ[9]    = {-300,-100,-50,-25,0,+25,+50,+100,+300};//in cm (A+C)/2 average
		const Double_t TPCpileupMult[5] = {0,400,800,1600,10000};//pileup contributors A+C sum

		const TString parnames[1] = {"El"};
		for(Int_t ipar=0;ipar<1;ipar++){
			for(Int_t i=0; i<2; ++i){
				THnSparseD *hs_PID = dynamic_cast<THnSparseD*>(histos->GetHist(Form("Track_%s",AliDielectron::TrackClassName(i)),Form("hsPID_V0%s",parnames[ipar].Data())));
				if(!hs_PID) continue;
				hs_PID->SetBinEdges(0,NSPD0);
				hs_PID->SetBinEdges(1,TPCpileupZ);
				hs_PID->SetBinEdges(2,TPCpileupMult);
				hs_PID->GetAxis(0)->SetTitle("N_{cls}^{SPD0}");
				hs_PID->GetAxis(1)->SetTitle("Z_{puv} (cm)");
				hs_PID->GetAxis(2)->SetTitle("M_{puv}");
			}
		}
	}

	const Int_t Nmee = 150;
	Double_t mee[Nmee] = {};
	for(Int_t i=0  ;i<110 ;i++) mee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 1.1 GeV/c2, every 0.01 GeV/c2
	for(Int_t i=110;i<Nmee;i++) mee[i] = 0.1  * (i-110) +  1.1;//from 1.1 to 5 GeV/c2, evety 0.1 GeV/c2

	const Int_t NpTee = 121;
	Double_t pTee[NpTee] = {};
	for(Int_t i=0  ;i<10   ;i++) pTee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 0.09 GeV/c, every 0.01 GeV/c
	for(Int_t i=10 ;i<110  ;i++) pTee[i] = 0.1  * (i- 10) +  0.1;//from 0.1 to 10 GeV/c, evety 0.1 GeV/c
	for(Int_t i=110;i<NpTee;i++) pTee[i] = 1.0  * (i-110) + 10.0;//from 10 to 20 GeV/c, evety 1.0 GeV/c

	TVectorD *v_mee = new TVectorD(Nmee);
	for(Int_t i=0;i<Nmee;i++) (*v_mee)[i] = mee[i];

	TVectorD *v_pTee = new TVectorD(NpTee);
	for(Int_t i=0;i<NpTee;i++) (*v_pTee)[i] = pTee[i];

	TVectorD *v_phiv = AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
	TVectorD *v_lmee = AliDielectronHelper::MakeLinBinning( 20, 0., 0.2);

	//histos->UserHistogram("Pair","hMvsPt"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", (TVectorD*)v_mee->Clone(), (TVectorD*)v_pTee->Clone(), AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
	//histos->UserHistogram("Pair","hMvsPt"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 500,0.,5.,200,0.,20, AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
	//histos->UserHistogram("Pair","hMvsPtvsPhiV","m_{ee} vs. p_{T,ee} vs. #varphi_{V};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);#varphi_{V} (rad.)", (TVectorD*)v_lmee->Clone(), (TVectorD*)v_pTee->Clone(), (TVectorD*)v_phiv->Clone(), AliDielectronVarManager::kM, AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);
	//histos->UserHistogram("Pair","hMvsPhiV"    ,"m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", (TVectorD*)v_phiv->Clone(), (TVectorD*)v_lmee->Clone(), AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kM);

	//add histograms to pair classes
	//histos->UserHistogram("Pair","hMvsPtvsPhiV"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 500,0.,5.,100,0.,20,90,0.,TMath::Pi(), AliDielectronVarManager::kM, AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);

	//histos->UserHistogram("Pair","hMvsPtvsPhiV"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)",v_mee, v_pTee,v_phiv, AliDielectronVarManager::kM, AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);
	//histos->UserHistogram("Pair","hEtaPhi" ,"pair #eta vs #varphi;#varphi (rad.);#eta",72,0.,TMath::TwoPi(), 200,-1.,1.,AliDielectronVarManager::kPhi, AliDielectronVarManager::kEta);

	//const Int_t Ndim_pair = 5;//mee, pT,eta,phi,phiV
	//Int_t Nbin_pair[Ndim_pair]    = {149, 110, 20,          36,          72};
	//Double_t xmin_pair[Ndim_pair] = {  0,   0, -1,           0,           0};
	//Double_t xmax_pair[Ndim_pair] = {  5,  20, +1, TMath::Pi(), TMath::Pi()};
	//UInt_t var_pair[Ndim_pair] = {AliDielectronVarManager::kM,AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi, AliDielectronVarManager::kPhivPair};
	const Int_t Ndim_pair = 3;//mee, pT, phiV
	Int_t Nbin_pair[Ndim_pair]    = {149, 120,         100};
	Double_t xmin_pair[Ndim_pair] = {  0,   0,           0};
	Double_t xmax_pair[Ndim_pair] = {  5,  20, TMath::Pi()};
	UInt_t var_pair[Ndim_pair] = {AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair};

	histos->UserSparse("Pair","hsPair","ee pairs",Ndim_pair,Nbin_pair,xmin_pair,xmax_pair,var_pair);

	for (Int_t i=0; i<8; ++i){
		THnSparseD *hs_pair = dynamic_cast<THnSparseD*>(histos->GetHist(Form("Pair_%s",AliDielectron::PairClassName(i)),"hsPair"));
		if(!hs_pair) continue;
		hs_pair->SetBinEdges(0,mee);
		hs_pair->SetBinEdges(1,pTee);
		hs_pair->GetAxis(0)->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
		hs_pair->GetAxis(1)->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
		hs_pair->GetAxis(2)->SetTitle("#it{#varphi}_{V} (rad.)");
		//hs_pair->GetAxis(2)->SetTitle("#it{#eta}_{ee}");
		//hs_pair->GetAxis(3)->SetTitle("#it{#varphi}_{ee} (rad.)");
		//hs_pair->GetAxis(4)->SetTitle("#it{#varphi}_{V} (rad.)");
	}

	if(die->GetMCSignals()) {
		for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
			THnSparseD *hs_pair = dynamic_cast<THnSparseD*>(histos->GetHist(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()),"hsPair"));
			if(!hs_pair) continue;
			hs_pair->SetBinEdges(0,mee);
			hs_pair->SetBinEdges(1,pTee);
			hs_pair->GetAxis(0)->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
			hs_pair->GetAxis(1)->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
			hs_pair->GetAxis(2)->SetTitle("#it{#varphi}_{V} (rad.)");
		}
	}

	die->SetHistogramManager(histos);

}
//___________________________________________________________________________________________________
void SetupMCSignals(AliDielectron *die)
{

	AliDielectronSignalMC* LFSig = new AliDielectronSignalMC("LF", "LFSignal"); ///all LF+J/psi
  LFSig->SetLegPDGs(11,-11);
  LFSig->SetMotherPDGs(600,600);
  LFSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  LFSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  LFSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  LFSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  LFSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  LFSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(LFSig);

//	AliDielectronSignalMC* pi0Sig = new AliDielectronSignalMC("pi0", "pi0Signal"); ///pi0 dalitz pairs
//  pi0Sig->SetLegPDGs(11,-11);
//  pi0Sig->SetMotherPDGs(111,111);
//  pi0Sig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  pi0Sig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  pi0Sig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  pi0Sig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  pi0Sig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  pi0Sig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(pi0Sig);
//
//  //AliDielectronSignalMC* pi0All = new AliDielectronSignalMC("pi0", "pi0All"); ///pi0 dalitz pairs (also from secondary)
//  //pi0All->SetLegPDGs(11,-11);
//  //pi0All->SetMotherPDGs(111,111);
//  //pi0All->SetMothersRelation(AliDielectronSignalMC::kSame);
//  //pi0All->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  //pi0All->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  //pi0All->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  //pi0All->SetFillPureMCStep(kTRUE);
//  //die->AddSignalMC(pi0All);
//
//  AliDielectronSignalMC* etaSig = new AliDielectronSignalMC("Eta", "etaSignal"); ///eta dalitz pairs
//  etaSig->SetLegPDGs(11,-11);
//  etaSig->SetMotherPDGs(221,221);
//  etaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  etaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  etaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  etaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  etaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  etaSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(etaSig);
//
//
//  AliDielectronSignalMC* etaprimeSig = new AliDielectronSignalMC("Etaprime", "etaprimeSignal"); ///etaprime pairs
//  etaprimeSig->SetLegPDGs(11,-11);
//  etaprimeSig->SetMotherPDGs(331,331);
//  etaprimeSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  etaprimeSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  etaprimeSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  etaprimeSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  etaprimeSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  etaprimeSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(etaprimeSig);
//
//
//  AliDielectronSignalMC* rhoSig = new AliDielectronSignalMC("Rho", "rhoSignal"); ///rho pairs
//  rhoSig->SetLegPDGs(11,-11);
//  rhoSig->SetMotherPDGs(113,113);
//  rhoSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  rhoSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  rhoSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  rhoSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  rhoSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  rhoSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(rhoSig);
//
//  AliDielectronSignalMC* omegaSig = new AliDielectronSignalMC("Omega", "omegaSignal"); ///omega pairs
//  omegaSig->SetLegPDGs(11,-11);
//  omegaSig->SetMotherPDGs(223,223);
//  omegaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  omegaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  omegaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  omegaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  omegaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  omegaSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(omegaSig);
//
//  AliDielectronSignalMC* phiSig = new AliDielectronSignalMC("Phi", "phiSignal"); ///phi pairs
//  phiSig->SetLegPDGs(11,-11);
//  phiSig->SetMotherPDGs(333,333);
//  phiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  phiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  phiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  phiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  phiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  phiSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(phiSig);
//
//  AliDielectronSignalMC* jpsiSig = new AliDielectronSignalMC("Jpsi", "jpsiSignal"); ///jpsi pairs
//  jpsiSig->SetLegPDGs(11,-11);
//  jpsiSig->SetMotherPDGs(433,433);
//  jpsiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
//  jpsiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  jpsiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
//  jpsiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  jpsiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  jpsiSig->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(jpsiSig);

  //AliDielectronSignalMC* pcmSig = new AliDielectronSignalMC("GammaConv", "pcmSignal"); ///gamma conversion
  //pcmSig->SetLegPDGs(11,-11);
  //pcmSig->SetMotherPDGs(22,22);
  //pcmSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  //pcmSig->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  //pcmSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  //pcmSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //pcmSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  //pcmSig->SetFillPureMCStep(kTRUE);
  //die->AddSignalMC(pcmSig);


  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("chadrons","di-electrons from open charm hadrons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharm->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenCharm->SetCheckStackForPDG(kTRUE);
  diEleOpenCharm->SetPDGforStack(503);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharm);

	//#################### D-Mesons
//  AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
//  diEleOpenCharmCharged->SetLegPDGs(11,-11);
//  diEleOpenCharmCharged->SetMotherPDGs(401,401);
//  diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
//  diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
//  diEleOpenCharmCharged->SetPDGforStack(503);
//  diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenCharmCharged);
//
//  AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
//  diEleOpenCharmNeutral->SetLegPDGs(11,-11);
//  diEleOpenCharmNeutral->SetMotherPDGs(405,405);
//  diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
//  diEleOpenCharmNeutral->SetPDGforStack(503);
//  diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenCharmNeutral);


  //b hadrons
  AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("bhadrons","di-electrons from open beauty hadrons");  // dielectrons originating from beauty hadrons
  diEleOpenB->SetLegPDGs(11,-11);
  diEleOpenB->SetMotherPDGs(502,502);
  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenB);

//  //B meson (3)
//  AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
//  diEleOneOpenB->SetLegPDGs(11,-11);
//  diEleOneOpenB->SetMotherPDGs(401,501);
//  diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOneOpenB->SetGrandMotherPDGs(501,0);
//  diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOneOpenB);
//
//  //B meson (1)(1)
//  AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
//  diEleOpenB->SetLegPDGs(11,-11);
//  diEleOpenB->SetMotherPDGs(501,501);
//  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenB);
//
//  //B meson (2)(2)
//  AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
//  diEleOpenBtoD->SetLegPDGs(11,-11);
//  diEleOpenBtoD->SetMotherPDGs(401,401);
//  diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenBtoD->SetGrandMotherPDGs(501,501);
//  diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
//  diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  die->AddSignalMC(diEleOpenBtoD);
//
//  //B meson (1)(2)
//  AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
//  diEleOpenBandBtoD->SetLegPDGs(11,11);
//  diEleOpenBandBtoD->SetMotherPDGs(401,501);
//  diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
//  diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//  diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
//  diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
//  //do i need this?
//  diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
//  die->AddSignalMC(diEleOpenBandBtoD);

}
//___________________________________________________________________________________________________
