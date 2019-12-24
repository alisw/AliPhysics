//class LMEECutlib{
//
//
//	ClassDef();
//};

void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void SetupCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition);
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliAnalysisCuts *SetupTrackCuts(Int_t cutDefinition);  //Setup track cuts independent of AOD/ESD
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
const AliDielectronEventCuts *GetEventCuts();


Bool_t kMix = kTRUE;
Bool_t IsESD = kTRUE;
Bool_t randomizeDau = kFALSE;

TString names("cut06_pf_pt400");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

Int_t GetN(){return nDie;}

AliDielectron* Config_dsekihat_lowmass_PbPb_ESD(
		const Int_t cutDefinition=1,
		const Int_t CenMin=0,
		const Int_t CenMax=101,
		const Int_t Nmix=10
		)
{
  // ESD handler?
  // IsESD = (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  //
  // Setup the instance of AliDielectron
  //
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("Track cuts: %s",name.Data()));

  if(kMix){
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -8., -6., -4., -2. , 0., 2., 4., 6., 8. , 10.");
    mix->AddVariable(AliDielectronVarManager::kCentralityNew,"0,10,30,50,70,90,101");
    mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
    mix->SetDepth(10);
    die->SetMixingHandler(mix);
  }//kMix

  // set track cuts
  SetupCuts(die,cutDefinition);

  //pairing with TLorentzVector
  die->SetUseKF(kFALSE);

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled

  AliDielectronVarCuts*  centCuts = new AliDielectronVarCuts("centCuts",Form("kPbPb%02d%02d",CenMin,CenMax));
  centCuts->AddCut(AliDielectronVarManager::kCentralityNew, (Float_t)CenMin, (Float_t)CenMax);
  die->GetEventFilter().AddCuts(centCuts);

  InitHistograms(die,cutDefinition);
  //die->SetNoPairing();

  return die;
}
//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
  // Setup the track cuts
  //AliDielectronVarCuts *nSharedClsITS = new AliDielectronVarCuts("nSharedClsITS","nSharedClsITS");
  //nSharedClsITS->AddCut(AliDielectronVarManager::kNclsSITS,               1.0,   6.0, kTRUE);

  //This definition is bogus


  die->GetTrackFilter().AddCuts(SetupPreFilterESDtrackCuts(cutDefinition));
  die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
 
  AliDielectronVarCuts *phiVcut = new AliDielectronVarCuts("phiVcut","phiVcut");
  //phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.0, 3.2,kTRUE);
  //phiVcut->AddCut(AliDielectronVarManager::kM       , 0.0,0.05,kTRUE);
  phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.0, 3.2);
  phiVcut->AddCut(AliDielectronVarManager::kM       , 0.0,0.05);
  //phiVcut->SetCutType(AliDielectronVarCuts::kAny);
  //die->GetPairPreFilter().AddCuts(phiVcut);
 
// AliDielectronVarCuts *OpAngle = new AliDielectronVarCuts("OpAngle","OpAngle");
//  OpAngle->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050);
//  die->GetPairPreFilter().AddCuts(OpAngle);
}
//______________________________________________________________________________________
AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  AliDielectronPID *pid = new AliDielectronPID();
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
  pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -2. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  return pid;
}
//______________________________________________________________________________________
AliAnalysisCuts *SetupTrackCuts(Int_t cutDefinition)
{
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?
  //
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2,   100.);

  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);

  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);

  varCuts->AddCut(AliDielectronVarManager::kNclsITS,        5.0,   100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.50);
  //trkFilter->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  trkFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);


  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,        100.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 160.0);  //crossed rows?
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.5,   100.);  // fraction crossed rows/findable clusters in the TPC, as done in AliESDtrackCuts
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0,   4.0);

  cuts->AddCut(trkFilter);
  cuts->AddCut(trkCuts);
  cuts->AddCut(varCuts);

  AliAnalysisCuts* TrackCuts=0x0;
  TrackCuts = cuts;

  return TrackCuts;

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition){
  AliESDtrackCuts *fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);//bit4
  fesdTrackCuts->SetMaxDCAToVertexXY(2.4);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.2);
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  return fesdTrackCuts;
}
//______________________________________________________________________________________
AliESDtrackCuts *SetupPreFilterESDtrackCuts(Int_t cutDefinition){
  AliESDtrackCuts *fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);//bit4
  fesdTrackCuts->SetMaxDCAToVertexXY(2.4);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.2);
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  return fesdTrackCuts;
}
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
	for (Int_t i=0; i<3; ++i){
		histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
	}
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));

  //if(die->GetMCSignals()) for(Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));

  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",100,-50,50,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","hCentralityV0M","centrality;centrality V0M (%)",101,0,101,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","hNclsITS0","Number of clusters on ITS0;N_{cls}^{ITS0}",1000,0,1e+4,AliDielectronVarManager::kNclsITS0);
  histos->UserHistogram("Event","hNclsITS1","Number of clusters on ITS1;N_{cls}^{ITS1}",1000,0,1e+4,AliDielectronVarManager::kNclsITS1);
  histos->UserHistogram("Event","hNclsITS2","Number of clusters on ITS2;N_{cls}^{ITS2}",1000,0,1e+4,AliDielectronVarManager::kNclsITS2);
  histos->UserHistogram("Event","hNclsITS3","Number of clusters on ITS3;N_{cls}^{ITS3}",1000,0,1e+4,AliDielectronVarManager::kNclsITS3);
  histos->UserHistogram("Event","hNclsITS4","Number of clusters on ITS4;N_{cls}^{ITS4}",1000,0,1e+4,AliDielectronVarManager::kNclsITS4);
  histos->UserHistogram("Event","hNclsITS5","Number of clusters on ITS5;N_{cls}^{ITS5}",1000,0,1e+4,AliDielectronVarManager::kNclsITS5);

  //add histograms to track class
	const Int_t Ndim = 3;
	Int_t Nbin[Ndim]    = {100, 20,            100};
	Double_t xmin[Ndim] = {  0, -1,              0};
	Double_t xmax[Ndim] = { 10, +1, TMath::TwoPi()};
	UInt_t var[Ndim]  = {AliDielectronVarManager::kPt, AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi};
  //histos->UserSparse("Track",Ndim,Nbin,xmin,xmax,var);
  histos->UserSparse("Track","hsPtEtaPhi","p_{T} vs. #eta vs. #varphi;",Ndim,Nbin,xmin,xmax,var);
  histos->UserHistogram("Track","hDCAxyz","DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)",100,-5.,5,100,-5.,5,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","hNclsTPC","Number of clusters TPC; N_{cls}^{TPC};Number of tracks",161,-0.5,160.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","hNcrTPC","Number of crossed rows TPC;N_{crossed rows}^{TPC};Number of track",161,-0.5,160.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","hRatioNcrToNf","N_{cr}/N_{f} TPC;N_{cr}/N_{f} TPC;Number of tracks",200,0.,2.0,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","hChi2TPC","#chi2/N_{cls}^{TPC};#chi^{2}/N_{cls}^{TPC};Number of tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","hNclsITS","Number of clusters ITS; N_{cls}^{ITS};Number of tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","hChi2ITS","chi2/N_{cls}^{ITS};#chi^{2}/N_{cls}^{ITS};Number of tracks",100,0.,10,AliDielectronVarManager::kITSchi2Cl);

  histos->UserHistogram("Track","hTPCdEdxvsPin","TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx (a.u.)",1000,0.,10.,500,0.,500.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","hITSdEdxvsPin","ITS dE/dx vs. p_{in};p_{in} (GeV/c);ITS dE/dx (a.u.)",1000,0.,10.,500,0.,500.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","hTOFbetavsPin","TOF #beta vs. p_{in};p_{in} (GeV/c);TOF #beta"       ,1000,0.,10.,120,0.,1.2 ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);

	const Int_t Ndim_PID = 5;
	Int_t Nbin_PID[Ndim_PID]    = {100, 20, 100, 100, 100};
	Double_t xmin_PID[Ndim_PID] = {  0, -1,  -5,  -5,  -5};
	Double_t xmax_PID[Ndim_PID] = { 10, +1,  +5,  +5,  +5};
	UInt_t var_PID[Ndim_PID] = {AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kITSnSigmaEle, AliDielectronVarManager::kTOFnSigmaEle};
  //histos->UserSparse("Track",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID);
  histos->UserSparse("Track","hsPID","hsPID",Ndim_PID,Nbin_PID,xmin_PID,xmax_PID,var_PID);

  //histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;p;TPCsigmaEle",1000,0.,10.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  //histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;p;ITSsigmaEle",1000,0.,10.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
	//histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;p;TOFsigmaEle",1000,0.,10.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

  // histos->UserHistogram("Track","TPCnSigma_EtaEle","TPC number of sigmas Electrons vs Eta;Eta;TPCsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","ITSnSigma_EtaEle","ITS number of sigmas Electrons vs Eta;Eta;ITSsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  // histos->UserHistogram("Track","TOFnSigma_EtaEle","TOF number of sigmas Electrons vs Eta;Eta;TOFsigmaEle", 200,-2.,2.,400,-20., 20. ,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);

  // histos->UserHistogram("Track","TPCnSigma_PhiEle","TPC number of sigmas Electrons vs Phi;Phi;TPCsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","ITSnSigma_PhiEle","ITS number of sigmas Electrons vs Phi;Phi;ITSsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
  // histos->UserHistogram("Track","TOFnSigma_PhiEle","TOF number of sigmas Electrons vs Phi;Phi;TOFsigmaEle", 200,0.,7.,400,-20., 20. ,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);

  //add histograms to pair classes
  histos->UserHistogram("Pair", "hMvsPt","m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", 500,0.,5.,200,0.,20, AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair", "hEtaPhi","pair #eta vs #varphi;#varphi (rad.);#eta",100,0.,TMath::TwoPi(), 200,-1.,1.,AliDielectronVarManager::kPhi, AliDielectronVarManager::kEta);
  histos->UserHistogram("Pair", "hMvsPhiV","m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})",100,0.,TMath::Pi(),100,0.,0.1,AliDielectronVarManager::kPhivPair,AliDielectronVarManager::kM);

  // histos->UserHistogram("Pair",
  // 	                   	"InvMass_OpAngle","InvMass_OpAngle;Invariant Mass;Opening angle",
  // 	                   	100, 0., 0.5 ,160, 0., 3.2,
  // 	                   	AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  // histos->UserHistogram("Pair",
  //                       "Y","Y;counts;Y",
  //                       60, -1.2 , 1.2,
  //                       AliDielectronVarManager::kY);

  // histos->UserHistogram("Pair","DCA_lin_sqr","#it{DCA}_{ee,lin} vs. #it{DCA}_{ee,sqr}",100,0.,10.,100,0.,10., AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);

  //3d histos for pair dca
  //histos->UserHistogram("Pair","InvMass_DCAsigma_pPt","", 200,0.,4, 50,0.,20., 80,0.,8., AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPt);

  die->SetHistogramManager(histos);

}
//___________________________________________________________________________________________________
const AliDielectronEventCuts *GetEventCuts(){
  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);
  return eventCuts;
}
//___________________________________________________________________________________________________
