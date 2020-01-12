void InitHistograms(AliDielectron *die, Int_t cutDefinition);

//TString arrNames[] = {
//  "DefaultPID_wpf"
// ,"TPCHadRejTOFrecPID_wpf"
// ,"TPCHadRejPID_wpf"
// ,"TOFrecoverPID_wpf"
// ,"noPID_wpf"
//};

TString TrackCutnames[] = {
  "DefaultTrackCut"
// ,"LooseTrackCut"
// ,"TightTrackCut"
};
const Int_t nTC = sizeof(TrackCutnames)/sizeof(TrackCutnames[0]);
Int_t GetNTC(){return nTC;}

TString PIDnames[] = {
  "DefaultPID"
 ,"ITSTPChadrejORTOFrec"
 ,"ITSTPChadrej"
 ,"ITSTOFrecover"
 ,"noPID"
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
    const Bool_t isAOD,
    const Float_t PtMin ,
    const Float_t PtMax ,
    const Float_t EtaMin,
    const Float_t EtaMax
		)
{
  // Setup the instance of AliDielectron
  // create the actual framework object
  TString name = Form("%s_%s_%s",TrackCutnames[cutDefinitionTC].Data(),PIDnames[cutDefinitionPID].Data(),PFnames[cutDefinitionPF].Data());

  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("cuts: %s",name.Data()));
  die->Print();

  // set track cuts
  LMEECutLib *lib = new LMEECutLib(name,isAOD);
  die->GetTrackFilter().AddCuts(lib->SetupTrackCuts(PtMin,PtMax,EtaMin,EtaMax));//AliDielectronCutGroup
  die->GetTrackFilter().AddCuts(lib->SetupPIDCuts());//AliDielectronCutGroup
  if(!isAOD) die->GetTrackFilter().AddCuts(lib->SetupESDtrackCuts());

  if(name.Contains("wPF",TString::kIgnoreCase)) die->GetPairPreFilter().AddCuts(lib->SetupPhiVPreFilter());

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);

	if(name.Contains("postPIDCalib",TString::kIgnoreCase) || name.Contains("noPID",TString::kIgnoreCase)) die->SetNoPairing(kTRUE);
	else                                                                                                  die->SetNoPairing(kFALSE);

	const Int_t Ndim    = 5;//pin,eta,NSPD0,puvZ,puzM
	Int_t Nbin[Ndim]    = {100,  20,   4,   8,   4};
	Double_t xmin[Ndim] = {  0,-1.0,   0,-250,   0};
	Double_t xmax[Ndim] = { 10,+1.0,2e+4,+250,1e+4};
	THnF *hs_mean  = new THnF("hs_mean" ,"hs_mean" ,Ndim,Nbin,xmin,xmax);//bin index for THn starts from 0!
	THnF *hs_width = new THnF("hs_width","hs_width",Ndim,Nbin,xmin,xmax);//bin index for THn starts from 0!
	const Int_t Ntot = 100 * 20 * 4 * 8 * 4;
  const Int_t Nall = hs_mean->GetNbins();//include over/underflow bins (i.e. +2 for each axis)
  printf("Ntot = %d , Nall = %d\n",Ntot,Nall);
	for(Int_t i=0;i<Nall;i++){
		hs_mean ->SetBinContent((Long64_t)i,-2);
		hs_width->SetBinContent((Long64_t)i,1.5);
	}

//	die->SetPIDCaibinPU(kTRUE);
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kPion,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kPion,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kPion,hs_mean ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
//	die->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kPion,hs_width,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
  //this should be in AddTask.C to avoid many calling open file.

  InitHistograms(die,0);

  return die;
}
//______________________________________________________________________________________
//void SetupCuts(AliDielectron *die, Int_t cutDefinition)
//{
//  // Setup the track cuts
//
//  die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
//  die->GetTrackFilter().AddCuts(SetPIDCuts(cutDefinition));
//  die->GetTrackFilter().AddCuts(SetupTrackCuts(cutDefinition));
//
//	TString cutname = arrNames[cutDefinition];
//	if(cutname.Contains("wpf",TString::kIgnoreCase)){
//		printf("PhiV pre-filter will be used!\n");
//		AliDielectronVarCuts *phiVcut = new AliDielectronVarCuts("phiVcut","phiVcut");
//		phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.5, 3.2);
//		phiVcut->AddCut(AliDielectronVarManager::kM       , 0.0,0.05);
//		die->GetPairPreFilter().AddCuts(phiVcut);
//	}
//	else if(cutname.Contains("wopf",TString::kIgnoreCase)){
//		printf("No PhiV pre-filter will be used!\n");
//	}
//	else {
//		printf("No PhiV pre-filter will be used!\n");
//	}
//
//}
//______________________________________________________________________________________
//AliAnalysisCuts *SetupTrackCuts(Int_t cutDefinition)
//{
//
//	AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
//	TString cutname = arrNames[cutDefinition];
//	if(!cutname.Contains("postPIDCalib",TString::kIgnoreCase)){//normal analysis
//
//		//AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
//
//		/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
//		AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
//		trkFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
//		//trkFilter->SetRequireITSRefit(kTRUE);
//		//trkFilter->SetRequireTPCRefit(kTRUE);
//
//		/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
//		AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
//
//		if(cutname.Contains("Pt200MeV"))      varCuts->AddCut(AliDielectronVarManager::kPt,0.20,10.);
//		else if(cutname.Contains("Pt150MeV")) varCuts->AddCut(AliDielectronVarManager::kPt,0.15,10.);
//		else                                  varCuts->AddCut(AliDielectronVarManager::kPt,0.20,10.);
//
//		if(cutname.Contains("Eta08"))      varCuts->AddCut(AliDielectronVarManager::kEta,-0.8,+0.8);
//		else if(cutname.Contains("Eta09")) varCuts->AddCut(AliDielectronVarManager::kEta,-0.9,+0.9);
//		else                               varCuts->AddCut(AliDielectronVarManager::kEta,-0.8,+0.8);
//
//		varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
//		varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,3.0);
//
//		varCuts->AddCut(AliDielectronVarManager::kNclsITS,  4.0,6.0);
//		varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);
//
//		//varCuts->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//should not be used in 2018 PbPb analyses
//		varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
//		varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
//		varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);
//
//		//ITS shared cluster cut
//		varCuts->AddCut(AliDielectronVarManager::kNclsSITS,-0.1,0.1);//only accept 0
//
//		cuts->AddCut(trkFilter);
//		cuts->AddCut(varCuts);
//
//	}
//	else{
//		AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
//		gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
//		// gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,  0.25, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,  0.05, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,  0.05, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,  0.05, kFALSE);
//		gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35, 0.35, kFALSE); // should increase purity...
//		gammaV0Cuts->SetExcludeTracks(kFALSE);
//
//		/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
//		AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
//
//		varCuts->AddCut(AliDielectronVarManager::kPt,0.15,10.);
//		varCuts->AddCut(AliDielectronVarManager::kEta,-0.9,+0.9);
//
//		varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
//		varCuts->AddCut(AliDielectronVarManager::kImpactParZ, -3.0,3.0);
//
//		//varCuts->AddCut(AliDielectronVarManager::kNclsITS,  4.0,6.0);
//		//varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,0.0,5.0);
//
//		varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);//crossed rows
//		varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,    2.);
//		varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);
//
//		//ITS shared cluster cut
//		varCuts->AddCut(AliDielectronVarManager::kNclsSFracITS,0.2,1.0);
//
//		cuts->AddCut(gammaV0Cuts);
//		cuts->AddCut(varCuts);
//	}
//
//  AliAnalysisCuts* TrackCuts=0x0;
//  TrackCuts = cuts;
//  return TrackCuts;
//}
////______________________________________________________________________________________
//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //Setup histogram classes
  AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(), die->GetTitle());

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

  histos->UserHistogram("Track","hTPCdEdxvsPin","TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx (a.u.)",1000,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","hITSdEdxvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);ITS dE/dx (a.u.)",1000,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","hTOFbetavsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);TOF #beta"       ,1000,0.,10.,240,0.,1.2 ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","hTPCnSigmaElvsPin","TPC dE/dx vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TPC}",1000,0.,10.,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","hITSnSigmaElvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{ITS}",1000,0.,10.,100,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","hTOFnSigmaElvsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);n#sigma_{e}^{TOF}",1000,0.,10.,140,-5,+5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

	//const Int_t Ndim_PID = 8;
	//Int_t Nbin_PID[Ndim_PID]    = {   20,     20,  100,100, 20, 100, 100, 100};
	//Double_t xmin_PID[Ndim_PID] = {    0,   -250,    0,  0, -1,  -5,  -5,  -5};
	//Double_t xmax_PID[Ndim_PID] = {20000,   +250,10000, 10, +1,  +5,  +5,  +5};
	//UInt_t var_PID[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kITSnSigmaEle, AliDielectronVarManager::kTOFnSigmaEle};
  //histos->UserSparse("Track","hsPID_V0El","hsPID_V0El",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID);
	//UInt_t var_PID_Pio[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPio, AliDielectronVarManager::kITSnSigmaPio, AliDielectronVarManager::kTOFnSigmaPio};
  //histos->UserSparse("Track","hsPID_V0Pi","hsPID_V0Pi",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pio);
	//UInt_t var_PID_Pro[Ndim_PID] = {AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaPro, AliDielectronVarManager::kITSnSigmaPro, AliDielectronVarManager::kTOFnSigmaPro};
  //histos->UserSparse("Track","hsPID_V0Pi","hsPID_V0Pi",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID_Pro);

  //histos->UserHistogram("Track","hAPplot","AP plot;#alpha;q_{T} (GeV/c)",100,-1.,+1,300,0.,0.3,AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);

  //add histograms to pair classes
  histos->UserHistogram("Pair","hMvsPt"  ,"m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 500,0.,5.,200,0.,20, AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","hEtaPhi" ,"pair #eta vs #varphi;#varphi (rad.);#eta",100,0.,TMath::TwoPi(), 200,-1.,1.,AliDielectronVarManager::kPhi, AliDielectronVarManager::kEta);
  histos->UserHistogram("Pair","hMvsPhiV","m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})",100,0.,TMath::Pi(),100,0.,0.1,AliDielectronVarManager::kPhivPair,AliDielectronVarManager::kM);
	histos->UserHistogram("Pair","hMvsOA"  ,"m_{ee} vs. opening angle;opening angle;m_{ee} (GeV/c^{2})",100,0.,TMath::Pi(),100,0.,0.1 ,AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);

  // histos->UserHistogram("Pair","DCA_lin_sqr","#it{DCA}_{ee,lin} vs. #it{DCA}_{ee,sqr}",100,0.,10.,100,0.,10., AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);

  //3d histos for pair dca
  //histos->UserHistogram("Pair","InvMass_DCAsigma_pPt","", 200,0.,4, 50,0.,20., 80,0.,8., AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPt);

  die->SetHistogramManager(histos);

}
//___________________________________________________________________________________________________
//const AliDielectronEventCuts *GetEventCuts(){
//  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Any && |vtxZ|<10 && ncontrib>0");
//  eventCuts->SetRequireVertex(kTRUE);
//  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
//  eventCuts->SetVertexZ(-10.,10.);
//  eventCuts->SetMinVtxContributors(1);
//	eventCuts->SetTimeRangeCut(kTRUE);
//  return eventCuts;
//}
////___________________________________________________________________________________________________
