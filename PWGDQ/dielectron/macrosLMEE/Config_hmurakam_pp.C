void InitHistograms(AliDielectron *die,Bool_t isMix);
void SetSignalsMC (AliDielectron *die);
void InitCF(AliDielectron *die, Int_t cutDefinition);
void SetupCuts(AliDielectron *die, Int_t cutDefinition);
AliAnalysisCuts *SetuptrackCuts(Int_t cutDefinition);
AliAnalysisCuts *SetPIDcuts(Int_t cutDefinition);

TString names=("cut1;cut2;cut3;cut4");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Int_t GetN(){return nDie;}

const AliDielectronEventCuts* GetEventCuts(Float_t centmin, Float_t centmax);

AliDielectron* Config_hmurakam_pp(Int_t cutDefinition=0, Bool_t isMC=kFALSE, Bool_t isMix=kTRUE)
{
  //
  // Setup the instance of AliDielectron
  //

  //=== create the actual framework object ========================
  TString name=arrNames->At(cutDefinition)->GetName();

  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));

  // set track cuts
  SetupCuts(die,cutDefinition);

  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //

  //  if(isMC) SetSignalsMC(die);
  InitHistograms(die,isMix);
  die->SetNoPairing(kFALSE); //do paring always

  return die;
}


//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
 
  //pairing with TLorentzVector
  die->SetUseKF(kFALSE); 
  die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
  die->GetTrackFilter().AddCuts(SetuptrackCuts(cutDefinition));

}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliAnalysisCuts *SetPIDcuts(Int_t cutDefinition){

  //PID1
  AliDielectronPID *pidTPCTOFreq = new AliDielectronPID("pidTPCTOFreq","pidTPCTOFreq");//<--> same as recoverTOF
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,     -3. , 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. , 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,      -3., 3., 0.4, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kP);// TOF required
  //PID2
  AliDielectronPID *pidTPCHadRej = new AliDielectronPID("pidTPCHadRej","pidTPCHadRej");//pure TPC pid
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,  -3., 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kPion,    -100., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,      -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kProton,    -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC

  AliDielectronCutGroup* combinedPIDcuts = new AliDielectronCutGroup("combinedPIDcuts","combinedPIDcuts",AliDielectronCutGroup::kCompOR);
  combinedPIDcuts->AddCut(pidTPCTOFreq);//PID1
  combinedPIDcuts->AddCut(pidTPCHadRej);//PID2

  AliAnalysisCuts* fancyCut=NULL;
  fancyCut = reinterpret_cast<AliAnalysisCuts*>(combinedPIDcuts);
  return fancyCut;

}

//______________________________________________________________________________________
//-----------------------------------track cuts-----------------------------------------
AliAnalysisCuts *SetuptrackCuts(Int_t cutDefinition){

  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;
  
  AliDielectronTrackCuts* trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  // pT and eta
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,   1e30);
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,   0.8);

  //TPC
  trackCutsDiel->SetRequireTPCRefit(kTRUE);

  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);//(1)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,  0.8,   1.5);//(2)
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,       0.0,   4.0);//(3)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//(4)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,    0.0,   0.4);//(5)

  //ITS
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,         3.0, 100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,       0.0,   4.5);

  //primary selection
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,   3.0);

  printf("Add shared cluster cut\n");
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, 1.0, 6.0, kTRUE);//accept no shared cls hit (default)

  //Pair cut Phiv rejection
  AliDielectronVarCuts* PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
  if(cutDefinition==0){//off
    printf("No phiv rejection\n");
    //    PhiV->SetCutType(1);//-->crash
  }else if(cutDefinition==1){//default
    printf("Current default 40MeV \n");
    PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.0 , 3.2,kTRUE);
    PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.100,kTRUE);
    //    PhiV->SetCutType(1);//-->crash
    PhiV->SetCutType(AliDielectronVarCuts::CutType::kAny);
  }else if(cutDefinition==2){//loose
    printf("loose 2.3 \n");
    PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.3 , 3.2,kTRUE);
    PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.100,kTRUE);
    //    PhiV->SetCutType(1);//-->crash
    PhiV->SetCutType(AliDielectronVarCuts::CutType::kAny);
  }else if(cutDefinition==3){//tight
    printf("loose TMath::Pi()/2.0 \n");
    PhiV->AddCut(AliDielectronVarManager::kPhivPair, TMath::Pi()/2.0 , 3.2,kTRUE);
    PhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.100,kTRUE);
    //    PhiV->SetCutType(1);//-->crash
    PhiV->SetCutType(AliDielectronVarCuts::CutType::kAny);
  }


  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);
  trackCuts->AddCut(PhiV);

  AliAnalysisCuts* returnCut = NULL;
  returnCut = reinterpret_cast<AliAnalysisCuts*>(trackCuts);

  return returnCut;
}

void InitHistograms(AliDielectron *die,Bool_t isMix)
{
  //Setup histogram classes
  AliDielectronHistos *histos= new AliDielectronHistos(die->GetName(),
						       die->GetTitle());

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");

  //Event class
  histos->AddClass("Event");

  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }

  //Pair classes
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    // Legs of final Pairs. Both charges together. No duplicate entries.
  }

  //ME and track rot
  if (isMix) {
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
  }

  //=== add histograms to event class ===============================
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","VtxZ","Vertex Z;Vertex Z [cm];N of events",150,-15.,15.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;N of Vertex Contributors;N of events",200,-0.5,199.5,AliDielectronVarManager::kNVtxContrib);
  // --- V0 info ---------------------------------------------
  histos->UserHistogram("Event","MultV0","Multiplicity V0;V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kMultV0);
  histos->UserHistogram("Event","EqMultV0","Equalized Multiplicity V0;Equalized V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kEqMultV0);
  histos->UserHistogram("Event","ChMultV0","Charged Multiplicity V0;Charged V0M amplitude",1000,-0.5,999.5,AliDielectronVarManager::kVZEROchMult);
  histos->UserHistogram("Event","CentralityV0M","Centrality V0;V0M percentile",300,-50,250,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","CentralityV0Mzoomed","Centrality V0 zoomed;V0M percentile",200,0,2,AliDielectronVarManager::kCentralityNew);
  //--- track checks (ITS) ---------------------------------------------
  histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;N of ITS clusters;N of tracks",10,-0.5,9.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS #chi^{2}/N of ITS clusters;N of tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","ITSSharedClusters","N of ITS shared clusters;N of shared clusters ITS;N of tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsSITS);
  //--- PID information ---------------------------------------------
  /* //--- ITS --------------------------------------------------------- */
  /* histos->UserHistogram("Track","ITS_dEdx_P","ITS dEdx;#it{p} (GeV/#it{c});ITS signal (arb units)",400,0.,20.,700,0.,700.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal); */
  /* histos->UserHistogram("Track","ITSnSigmaEle_P","ITS number of sigmas Electrons;#it{p} (GeV/#it{c});ITS number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle); */
  /* histos->UserHistogram("Track","ITSnSigmaPio_P","ITS number of sigmas Pions;#it{p} (GeV/#it{c});ITS number of sigmas Pions",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio); */
  /* histos->UserHistogram("Track","ITSnSigmaKao_P","ITS number of sigmas Kaons;#it{p} (GeV/#it{c});ITS number of sigmas Kaons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao); */
  /* histos->UserHistogram("Track","ITSnSigmaPro_P","ITS number of sigmas Protons;#it{p} (GeV/#it{c});ITS number of sigmas Protons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro); */
  /* histos->UserHistogram("Track","ITSnSigmaEle_Eta_P_lin","ITS number of sigmas Electrons vs Eta and P;#eta;n#sigma_{ele}^{ITS};#it{p} (GeV/#it{c})",200, -1., 1., 200, -10., 10., 100, 0., 5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP); */
  /* //--- TPC ---------------------------------------------------------- */
  /* histos->UserHistogram("Track","TPC_dEdx_P","TPC dEdx;#it{p} (GeV/#it{c});TPC signal (arb units)",400,0.,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal); */
  /* histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;#it{p} (GeV/#it{c});TPC number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle); */
  /* histos->UserHistogram("Track","TPCnSigmaEleRaw_P","raw TPC number of sigmas Electrons;#it{p} (GeV/#it{c});raw TPC number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEleRaw); */
  /* histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;#it{p} (GeV/#it{c});TPC number of sigmas Pions",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio); */
  /* histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;#it{p} (GeV/#it{c});TPC number of sigmas Kaons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao); */
  /* histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;#it{p} (GeV/#it{c});TPC number of sigmas Protons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro); */
  /* histos->UserHistogram("Track","TPCnSigmaEle_Eta_P_lin",";Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",200, -1., 1., 200, -10., 10., 100, 0., 5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kP); */
  /* histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_P_lin",";Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})",200, -1., 1., 200, -10., 10., 100, 0., 5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw,AliDielectronVarManager::kP); */
  /* //--- TOF ---------------------------------------------------------- */
  /* histos->UserHistogram("Track","TOFbeta","TOF beta;#it{p} (GeV/#it{c});TOF #beta",400,0.,20.,300,0.,3.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta); */
  /* histos->UserHistogram("Track","TOFbeta_Pt","TOF beta;#it{p}_{T} (GeV/#it{c});TOF #beta",400,0.,20.,300,0.,3.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFbeta); */
  /* histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;#it{p} (GeV/#it{c});TOF number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle); */
  /* histos->UserHistogram("Track","TOFnSigmaEleRaw_P","raw TOF number of sigmas Electrons;#it{p} (GeV/#it{c});raw TOF number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEleRaw); */
  /* histos->UserHistogram("Track","TOFnSigmaEle_Pt","TOF number of sigmas Electrons;#it{p} (GeV/#it{c});TOF number of sigmas Electrons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTOFnSigmaEle); */
  /* histos->UserHistogram("Track","TOFnSigmaPio_P","TOF number of sigmas Pions;#it{p} (GeV/#it{c});TOF number of sigmas Pions",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio); */
  /* histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;#it{p} (GeV/#it{c});TOF number of sigmas Kaons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao); */
  /* histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;#it{p} (GeV/#it{c});TOF number of sigmas Protons",400,0.,20.,600,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro); */
  /* histos->UserHistogram("Track","TOFnSigmaEle_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p} (GeV/#it{c})",200, -1., 1., 200, -10., 10., 100, 0., 5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle,AliDielectronVarManager::kP); */
  /* histos->UserHistogram("Track","TOFnSigmaEleRaw_Eta_P_lin","raw TOF number of sigmas Electrons vs Eta and P;Eta;n#sigma_{ele}^{TOF};#it{p} (GeV/#it{c})",200, -1., 1., 200, -10., 10., 100, 0., 5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEleRaw,AliDielectronVarManager::kP); */
  
  //=== add histograms to Pair classes ==============================
  
  //Version 4.0 2022/01/20 more finner binning
  //mee
  /* const Int_t Nmee = 801; */
  /* Double_t mee[Nmee] = {}; */
  /* for(Int_t i=0  ;i<Nmee ;i++) mee[i] = 0.005 * (i - 0) + 0.0;//from 0 to 4 GeV/c2, every 5 MeV/c2 */

  //Version 5.0 2022/02/20
  // mee
  /* const Int_t Nmee = 9; */
  /* Double_t mee[Nmee] = {0.00,  0.04, 0.08, 0.14, 0.35, 1.03, 2.80, 3.10, 4.00}; */
  /* TVectorD *v_mee = new TVectorD(Nmee); */
  /* for(Int_t i=0;i<Nmee;i++) (*v_mee)[i] = mee[i]; */

  //Version 5.0 2022/09/09 resolution study
  /* const Int_t Nmee = 4001; */
  /* Double_t mee[Nmee] = {}; */
  /* for(Int_t i=0  ;i<Nmee ;i++) mee[i] = 0.001 * (i - 0) + 0.0;//from 0 to 4 GeV/c2, every 1 MeV/c2 */
  /* TVectorD *v_mee = new TVectorD(Nmee); */
  /* for(Int_t i=0;i<Nmee;i++) (*v_mee)[i] = mee[i]; */

  /* //Version 6.0 2022/10/29 */
  /* // mee */
  /* const Int_t Nmee = 45; */
  /* Double_t mee[Nmee] = {0,00, */
  /* 			0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, */
  /* 			0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, */
  /* 			0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, */
  /* 			0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, */
  /* 			0.41, 0.40, 0.45, 0.52}; */
  /* TVectorD *v_mee = new TVectorD(Nmee); */
  /* for(Int_t i=0;i<Nmee;i++) (*v_mee)[i] = mee[i]; */

  /* // ptee */
  /* const int Nptee = 17; */
  /* Double_t ptee[Nptee] = {0.0, */
  /* 			  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, */
  /* 			  1.5, 2.0, 3.0, 4.0, 5.0, 6.0}; */
  /* TVectorD *v_ptee = new TVectorD(Nptee); */
  /* for(Int_t i=0;i<Nptee;i++) (*v_ptee)[i] = ptee[i]; */
  /* // phiv */
  /* TVectorD *v_phiv = AliDielectronHelper::MakeLinBinning(1, 0., TMath::Pi()); */
  /* histos->UserHistogram("Pair","InvMass_PtRebinned_PhivPair","InvMass:Pt:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#phi_{V}",v_mee,v_ptee,v_phiv,AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);//913 --> modified and use */

  //Version 4.0 2022/01/20 Harry suggestion --> test more finner binning
  //mee
  const Int_t Nmee = 801;
  Double_t mee[Nmee] = {};
  for(Int_t i=0  ;i<Nmee ;i++) mee[i] = 0.005 * (i - 0) + 0.0;//from 0 to 4 GeV/c2, every 5 MeV/c2
  TVectorD *v_mee = new TVectorD(Nmee);
  for(Int_t i=0;i<Nmee;i++) (*v_mee)[i] = mee[i];
  // ptee
  const int Nptee = 55;
  Double_t ptee[Nptee] = {0.0,
			  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
			  1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
			  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
			  3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
			  5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0,
			  7.0, 8.0, 9.0, 10.0};

  TVectorD *v_ptee = new TVectorD(Nptee);
  for(Int_t i=0;i<Nptee;i++) (*v_ptee)[i] = ptee[i];        
  // phiv
  TVectorD *v_phiv = AliDielectronHelper::MakeLinBinning(90, 0., TMath::Pi());
  histos->UserHistogram("Pair","InvMass_PtRebinned_PhivPair","InvMass:Pt:PhivPair;#it{m}_{ee} (GeV/#it{c}^{2});Pair #it{p}_{T} (GeV/#it{c});#phi_{V}",v_mee,v_ptee,v_phiv,AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);//913 --> modified and use

  //=== add histograms relevant for mixing ============================
  /* AliDielectronMixingHandler *mix = die->GetMixingHandler(); */
  /* if (kMixing && kPairing){ */
  /*   const Int_t nMix=mix->GetNumberOfBins(); */
  /*   histos->UserHistogram("Event","","", nMix,0.,(Double_t)nMix,AliDielectronVarManager::kMixingBin); */
  /*   histos->UserHistogram("Pair","","", nMix, 0, nMix, AliDielectronVarManager::kMixingBin); */
  /* } */

  
  //=== add histograms relevant for mixing ============================ --> crash 2023/01/04
  /* AliDielectronMixingHandler *mix=die->GetMixingHandler(); */
  /* if (isMix) {     */
  /*   const Int_t nMix= mix->GetNumberOfBins(); */
  /*   histos->UserHistogram("Event","","", nMix,0.,(Double_t)nMix,AliDielectronVarManager::kMixingBin); */
  /*   histos->UserHistogram("Pair","","", nMix, 0, nMix, AliDielectronVarManager::kMixingBin); */
  /* } */

  die->SetHistogramManager(histos);
    
}
//-----------------------------------not used------------------------------------------------
/* void InitCF(AliDielectron *die, Int_t cutDefinition){ */

/*   AliDielectronCF *cf = new AliDielectronCF (die->GetName(), die->GetTitle()); */

/*   //Add pair variables */
/*   cf->AddVariable(AliDielectronVarManager::kM,200,0.,4.);     // Inv. Mass */
/*   cf->AddVariable(AliDielectronVarManager::kPt,60,0.,6.);    // Pair pT */
/*   cf->AddVariable(AliDielectronVarManager::kPairDCAsigXY,50,0.,20.);     // Pair DCA (sigma) square sum */
/*   cf->AddVariable(AliDielectronVarManager::kPairLinDCAsigXY,100,0.,20.);  // Pair DCA (sigma) linear sum */
/*   cf->AddVariable(AliDielectronVarManager::kPairType,11,0.,11.);    // PairType */
/*   //Add leg variables */
/*   cf->AddVariable(AliDielectronVarManager::kPt,200, 0.,10.,kTRUE); */
/*   //  cf->AddVariable(AliDielectronVarManager::kPdgCode,50, -20.,20.,kTRUE); */
/*   //  cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,50,-600.,600.,kTRUE); */
/*   //  cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,50,-600.,600.,kTRUE); */

/*   cf->SetStepsForMCtruthOnly(kTRUE); */
/*   die->SetCFManagerPair(cf); */
/* } */
//AliDielectronEventCuts* GetEventCuts()
const AliDielectronEventCuts* GetEventCuts(Float_t centmin, Float_t centmax)
{

  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Any && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex(kTRUE);
  //  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny); // used in data analysis (i)Any or (ii)SPD which is default ?
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD ? 
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //  eventCuts->SetTimeRangeCut(kTRUE); // For 2018 only
  if(-1 < centmin && -1 < centmax){
    eventCuts->SetCentralityRange(centmin,centmax,kTRUE);
  }  
  return eventCuts;  
}
//_____________________________________________________________________________________________
void SetSignalsMC (AliDielectron *die){}
