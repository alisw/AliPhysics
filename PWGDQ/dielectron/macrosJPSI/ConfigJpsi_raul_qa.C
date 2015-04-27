void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
//QAtask

void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Add(AliDielectron *die, Int_t cutDefinition);
TVectorD *GetRunNumbers();


void SetupMCsignals(AliDielectron *die);

TVectorD *GetRunNumbers() {
  // returns a vector with the runnumber used in the period                                                                                                                                                                                                                                                            
  Double_t first=0;
  Double_t last =1;
	
	
  switch(iPeriod) {
  case k10b: first=114737; last=117223; break;
  case k10c: first=117777; last=121417; break;
  case k10d: first=121692; last=126438; break;
  case k10e: first=127102; last=130851; break;
  case k10f: first=130931; last=135031; break;
  case k10h: first=136831; last=139517; break;
  case k11a: first=141052; last=146974; break;
  case k11d: first=155838; last=159649; break;
  case k11h: first=165772; last=170718; break;
  case k12h: first=188720; last=192738; break;
  }
  // printf("iPeriod: %d \t %.0f-%.0f \n",iPeriod,first,last);                                                                                                                                                                                                                                                        
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
}


	enum ConfDef {kDefaultCuts,kHF,kLmee,kDefault_activevolume,kDefault_conversions,kDefault_conversions_wPID};
	TString names=("JPsi;kHFe;Lmee;Jpsi_activevolume;Jpsi_conversions;Jpsi_conversions_wPID");
	TObjArray *arrNames=names.Tokenize(";");
	const Int_t nDie=arrNames->GetEntries();



//______________________________________________________________________________________
//______________________________________________________________________________________
//______________________________________________________________________________________
//
// Here the configuration part starts
//
AliDielectron* ConfigJpsi_raul_qa(Int_t cutDefinition)
{
	//
	// Setup the instance of AliDielectron
	//
	
	
	// create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =
    new AliDielectron(Form("%s",name.Data()),
		      Form("Track cuts: %s",name.Data()));
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  if (cutDefinition ==ConfDef::kDefaultCuts){
    SetupPairCuts(die,cutDefinition);
  }
  
  if (cutDefinition ==ConfDef::kDefaultCuts || cutDefinition ==ConfDef::kLmee || cutDefinition ==ConfDef::kDefault_activevolume){
    SetupV0Cuts(die,cutDefinition);
  }
  
  //V0s to have a pure e+e- sample to check the TPC nsigma	
  if (cutDefinition ==ConfDef::kDefault_conversions || cutDefinition ==ConfDef::kDefault_conversions_wPID){
    SetupV0add(die,cutDefinition);
  }
  
	
  if (cutDefinition !=ConfDef::kDefaultCuts){
    die->SetNoPairing();
  }	
  
  //SetupV0Cuts(die,cutDefinition); 
  
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // Monte Carlo Signals
  if (hasMC) SetupMCsignals(die);
  // prefilter settings
  // die->SetPreFilterUnlikeOnly();//  die->SetNoPairing();//  die->SetPreFilterAllSigns();
  // cut QA
  // die->SetCutQA();
  
	/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  InitHistograms(die,cutDefinition);
  //No CF container for now
  // if (cutDefinition ==0){
  // InitCF(die,cutDefinition);
  // }
  //   AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  //   mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-8,-5,0,5,8,10");
  //   mix->AddVariable(AliDielectronVarManager::kPhi ,"0,3.14,6.3");
  //   mix->SetDepth(10);
  //  die->SetMixingHandler(mix);
  //
  
  return die;
}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //
  
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  //eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);
  
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
	//
	// Setup the track cuts
	//
	
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  
  //default quality cuts
  AliDielectronTrackCuts *refit=new AliDielectronTrackCuts("refit","refit");
  
  if (cutDefinition==ConfDef::kDefaultCuts){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if(cutDefinition==ConfDef::kHF){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if(cutDefinition==ConfDef::kLmee){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if (cutDefinition ==ConfDef::kDefault_activevolume){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if (cutDefinition ==ConfDef::kDefault_conversions){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    //    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if(cutDefinition ==ConfDef::kDefault_conversions_wPID){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    //  refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		
  }
  cuts->AddCut(refit);
  
  
  //pt and kink mother
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  
  if (cutDefinition==kDefaultCuts){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  else if (cutDefinition==kHF) {
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    //to be checked
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-0.5.,3.);
    if (AliDielectronVarManager::kTOFPIDBit >0.8){
      pt->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.);
    }
    //TPC PID ClusteR
    pt->AddCut(AliDielectronVarManager::kTPCsignalN,80.,160.);
		//rjim NTPCcluster cut
    pt->AddCut(AliDielectronVarManager::kNclsTPC,120.,160.);  
    pt->AddCut(AliDielectronVarManager::kNFclsTPCfCross,0.6,1.1);    
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
		//    pt->AddCut(AliDielectronVarManager::kNClsITS,4.,200); 
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-2.,2.);  
  }
  else if (cutDefinition==kLmee) {
    pt->AddCut(AliDielectronVarManager::kPt,0.2,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kTOFnSigmaEle,-3.,3.); 
    //rjim Nitsclusters cut
    //  pt->AddCut(AliDielectronVarManager::kNclsITS,4.,160.);
    //add ncrossed rows tpc instead cluster instead Ncluster
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kNFclsTPCfCross,0.8,1.0);    
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.); 
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if(cutDefinition==kDefault_activevolume){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    //    pt->AddCut(AliDielectronVarManager::kTPCactvol,120.,200.);
    // NTPCclusters
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    //impact parameter
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
    //  pt->AddCut(AliDielectronVarManager::kTOFPIDBit,0.8,2.0);
  }else if(cutDefinition==kDefault_conversions){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    //    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }else if(cutDefinition==kDefault_conversions_wPID){
    pt->AddCut(AliDielectronVarManager::kPt,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-3.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }
  cuts->AddCut(pt);
  
  
	
  /*
    
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv /
  AliDielectronVarCuts *varAccCuts   = new AliDielectronVarCuts("acc","acc");
  varAccCuts->AddCut(AliDielectronVarManager::kPt,           0.8, 1e30);
  varAccCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,  0.9);
  die->GetTrackFilter().AddCuts(varAccCuts);
  varAccCuts->Print();
  
  
  AliDielectronVarCuts   *varRecCuts = new AliDielectronVarCuts("VarRecCuts","VarRecCuts");
  varRecCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.,   160.);
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varRecCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varRecCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  //  varRecCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,     0.  ,   36.     ); // not defined in AOD
  varRecCuts->AddCut(AliDielectronVarManager::kKinkIndex0,.000001,1e30,kTRUE);
  
  AliDielectronTrackCuts *trkRecCuts = new AliDielectronTrackCuts("TrkRecCuts","TrkRecCuts");
  trkRecCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkRecCuts->SetRequireITSRefit(kTRUE);
  trkRecCuts->SetRequireTPCRefit(kTRUE);
  
  AliDielectronCutGroup  *grpRecCuts = new AliDielectronCutGroup("rec","rec",AliDielectronCutGroup::kCompAND);
  grpRecCuts->AddCut(trkRecCuts);
  grpRecCuts->AddCut(varRecCuts);
  die->GetTrackFilter().AddCuts(grpRecCuts);
  grpRecCuts->Print();
	 
  
  */ 
  
  //exclude conversion electrons selected by the tender
  //   AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  //   noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //   cuts->AddCut(noconv);
  //	cuts->Print();
}
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //
  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);//ktrue excludes tracks v0s, 
  //kfalse 
  // gammaV0Cuts->Print();
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
  //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
  die->GetTrackFilter().AddCuts(gammaV0Cuts);
  //gammaV0Cuts->Print();	// 
}

void SetupV0add(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //
  AliDielectronV0Cuts *gammaV0Add = new AliDielectronV0Cuts("IsGamma2","IsGamma2");
  gammaV0Add->SetPdgCodes(22,11,11);
  gammaV0Add->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Add->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Add->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
  // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35,0.35, kFALSE); // not sure if it works as expected
	
  gammaV0Add->SetExcludeTracks(kFALSE);//ktrue excludes tracks v0s, 
  
  //kfalse 
  // gammaV0Cuts->Print();
  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
  //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; && const Double_t cutQTG2 < 0.04;
  die->GetTrackFilter().AddCuts(gammaV0Add);
  //gammaV0Add->Print();	// 
}
//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
  // add conversion rejection
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  die->GetPairPreFilter().AddCuts(gammaCut);
  
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
			    die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  if(cutDefinition==0){
    //Pair classes
    // to fill also mixed event histograms loop until 10
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    }
  }
  //add MC signal histograms to track and pair class
  if(die->GetMCSignals()) {
    for(Int_t isig=0; isig<die->GetMCSignals()->GetEntriesFast(); isig++) {
      TString sigMCname = die->GetMCSignals()->At(isig)->GetName(); 
      
      // mc truth
      histos->AddClass(Form("Pair_%s_MCtruth",       sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s_MCtruth", sigMCname.Data())); 
      // mc reconstructed
      histos->AddClass(Form("Pair_%s",               sigMCname.Data()));
      histos->AddClass(Form("Track_Legs_%s",         sigMCname.Data())); 
    }
  }
  if(cutDefinition==0){
    
    //add histograms to event class
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",700,0.,700.,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",100,0.,10.,AliDielectronVarManager::kPairs);
    histos->UserHistogram("Event","","",4,0.,4.,AliDielectronVarManager::kNevents);    
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    // histos->UserHistogram("Event","","",AliDielectronHelper::MakeArbitraryBinning(run_numbers),AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","",AliDielectronVarManager::kNTrk, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","",AliDielectronVarManager::kPairs, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(700,0.,700),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(10000,0.,10),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPairs);    
  }
  
  //add histograms to Track classes
  //run number dependence
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kITSnSigmaEle);  
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(20,-1.0,1.0),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(63,0.,6.32),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPhi);
  //nsigma
  histos->UserHistogram("Track","","",
			200,0.2,20.,AliDielectronVarManager::kPt); 
  histos->UserHistogram("Track","","",
			80,-1.0,1.0,AliDielectronVarManager::kEta); 
  histos->UserHistogram("Track","","",
			400,0.0,4.0,AliDielectronVarManager::kPhi);  
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);    	
  histos->UserHistogram("Track","","",
			200,0.2,20.,100,0.,1.2,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,20.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);   
  histos->UserHistogram("Track","","",
			100,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal); 
  histos->UserHistogram("Track","","",
			144,0.0,6.285,200,0.0,200,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
			40,-1.0,1.0,100,0.0,200,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","","",
			100,-5,5.,100,-5.,5.,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
			100,-5,5.,100,-5.,5.,AliDielectronVarManager::kTOFnSigmaEle,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",
			160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",
			100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",
			160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","","",
			160,0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  histos->UserHistogram("Track","","",200,0.,200.,AliDielectronVarManager::kNclsTPC, 200,0.,20., AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","","",
			150,-15,15,160,-0.5,159.5,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCsignalN); 
  //  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(300,0.5,300.5),AliDielectronVarManager::kNclsTRD);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kNaccTrckltsEsd10Corr);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(5,-0.5,2.0),AliDielectronVarManager::kTOFPIDBit);
  //rjim findable cluster vs pt
  histos->UserHistogram("Track","","",200,0.,20.,161,-0.5,161.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  //rjim frac find vs pt
  histos->UserHistogram("Track","","",			
			200,0.0,20.0,160,0,1.1,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCrFrac);
  // histos jpsi 2nd part
  //check tofbit
  //  histos->UserHistogram("Track","","",3.,-0.5,2.5, 200,0.,20.,AliDielectronVarManager::kTOFPIDBit,AliDielectronVarManager::kPt);
  
  //jpsi nsigma vseta
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10,10,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10,10,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,-1.0,1.0,100,-10,10,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle); 
  //jpsi nsigma vs phi
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10,10,AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10,10,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			20,0.0,7.0,100,-10,10,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);
  
  //for now SPD tracklets but multiplicity should be implemented
  histos->UserHistogram("Track","","",
			100,0.0,100.,100,-10.,10.,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			100,0.0,100.,200,-20,20,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);      
  //
  
  histos->UserHistogram("Track","","",
			200,-20.,20.,100,-10,10,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kTPCnSigmaEle);
  
  histos->UserHistogram("Track","","",
			200,-20.,20.,200,0.2,20.,AliDielectronVarManager::kYsignedIn,AliDielectronVarManager::kPIn);
  
  //inner and outer read out TPC clusters
  histos->UserHistogram("Track","","",
			70,0.0,7.0,160,0.0,160.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCclsIRO);
  histos->UserHistogram("Track","","",
			70,0.0,7.0,160,0.0,160.0,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCclsORO);
  histos->UserHistogram("Track","","",
			10000,0.0,1.0,AliDielectronVarManager::kM);
  //  histos->UserHistogram("Track","","",
  //	200,0.0,20.0,10000,0.0,1.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kM);	
  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  
  //add histograms to Pair classes
  
  if(cutDefinition ==0){
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","","",
			  100,-1.,1.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",
			  100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,200,0.2,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
    
  }
  die->SetHistogramManager(histos);
  
}

//______________________________________________________________________________________
void InitCF(AliDielectron* die, Int_t cutDefinition)
{
	//
	// Setupd the CF Manager if needed
	//
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
	
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 1.0, 1.3, 2.0, 3.0, 5., 7.0, 10.0, 100.0");
  
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.3,0,0.3,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.8, 1.0, 1.1, 1.2, 1.3, 100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);
  
  //cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0,6,kTRUE);
  
  if (hasMC && 0){ //ATTENTION SWITCHED OFF
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);
  }
  
  if(hasMC) {
    //only steps for efficiencies
    cf->SetStepsForMCtruthOnly();
		//only in this case write MC truth info
    if (cutDefinition==0){
      cf->SetStepForMCtruth();
    }
  }
  // cf->SetStepsForSignal();
	
  die->SetCFManagerPair(cf);
}

//______________________________________________________________________________________
void SetupMCsignals(AliDielectron *die){
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
  
  AliDielectronSignalMC* promptJpsiNonRad = new AliDielectronSignalMC("promptJpsiNonRad","Prompt J/psi non-Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiNonRad->SetLegPDGs(11,-11);
  promptJpsiNonRad->SetMotherPDGs(443,443);
  promptJpsiNonRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiNonRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiNonRad->SetFillPureMCStep(kTRUE);
  promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
  die->AddSignalMC(promptJpsiNonRad);
  
  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiRad->SetLegPDGs(11,-11);
  promptJpsiRad->SetMotherPDGs(443,443);
  promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiRad->SetFillPureMCStep(kTRUE);
  promptJpsiRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
  die->AddSignalMC(promptJpsiRad);
	
}

