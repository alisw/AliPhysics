void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
//void AddHistsEleEff(AliDielectron *die);
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die,  Int_t cutDefinition);
//QAtask
void SetupV0cuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Add(AliDielectron *die, Int_t cutDefinition);
//Eta correction
void SetEtaCorrection(AliDielectron *die);
TVectorD *GetRunNumbers();

void SetupMCsignals(AliDielectron *die);

Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

Bool_t isAOD=kTRUE;

TVectorD *GetRunNumbers() {
  Double_t first=0;
  Double_t last =1;
		
  switch(iPeriod) {   
  case k10h: first=136831; last=139517; break;
  case k11h: first=165772; last=170718; break;
  }
  return (AliDielectronHelper::MakeLinBinning(last-first, first, last));
}

enum ConfDef {kDefaultAny,kDefaultFirst,kITSNcluster3,kITSNcluster4}; 
TString names=("JPsi_Any; JPsi_First; JPsi_ITS3; JPsi_ITS4");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

//______________________________________________________________________________________
//______________________________________________________________________________________
//______________________________________________________________________________________
//
// Here the configuration part starts
//
AliDielectron* ConfigJpsi_zzhou_lowpt(Int_t cutDefinition) 
{
  //
  // Setup the instance of AliDielectron
  //
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);     
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
 	
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()) {
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die =new AliDielectron(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));
 
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  //### 4 different cuts (1)cutDefinition ==ConfDef::kDefaultAny  (2)ConfDef::kDefaultFirst (3)ConfDef::kITSNcluster3  (4)ConfDef::kITSNcluster4
  SetupEventCuts(die,cutDefinition);     // ##### Event cuts setup
  if (hasMC)SetupMCsignals(die);
  SetupTrackCuts(die,cutDefinition);     // ##### Track cuts setup
  SetupPairCuts(die,cutDefinition);      // ##### Pair cuts setup
   //does hfe apply v0 cuts like jpsi? check! 
  SetupV0Cuts(die,cutDefinition);        // ##### V0 cuts setup

  //###### post pid correction #####-----------------------------
  SetEtaCorrection(die);
  //-------------------------------------------------------------


    
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // Monte Carlo Signals
  //  if (hasMC) SetupMCsignals(die);
  // prefilter settings
  // die->SetPreFilterUnlikeOnly();//  die->SetNoPairing();//  die->SetPreFilterAllSigns();
  // cut QA
  // die->SetCutQA();
  
	/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  InitHistograms(die,cutDefinition);
  
  //CF container for efficiencies
  InitCF(die,cutDefinition);


  //  #### Event Mixing Handler #####
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
  eventCuts->SetCentralityRange(40.,90.);
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
  
  //Jpsi part
  if (cutDefinition==ConfDef::kDefaultAny){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  }else if(cutDefinition==ConfDef::kDefaultFirst){
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
    refit->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  }else if(cutDefinition==ConfDef::kITSNcluster3){ 
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
  }else if(cutDefinition==ConfDef::kITSNcluster4){ 
    refit->SetRequireITSRefit(kTRUE);
    refit->SetRequireTPCRefit(kTRUE);
  }


  cuts->AddCut(refit);
  
    //pt and kink mother
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("ptCut","pt cut");
  
  if (cutDefinition==ConfDef::kDefaultAny){
    //pt->AddCut(AliDielectronVarManager::kP,1.,1.e30);
    pt->AddCut(AliDielectronVarManager::kPt,0.85,8.); // same as Ionut's cuts
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8); // same as Ionut's cuts
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.); //3.5, 1000.
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.); //3.0, 1000.
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);  
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }  else if(cutDefinition==ConfDef::kDefaultFirst) {
    pt->AddCut(AliDielectronVarManager::kPt,0.85,8.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);  
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  } else if (cutDefinition==ConfDef::kITSNcluster3){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,8.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kNclsITS,3.,6.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  } else if (cutDefinition==ConfDef::kITSNcluster4){
    pt->AddCut(AliDielectronVarManager::kPt,0.85,8.);
    pt->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
    pt->AddCut(AliDielectronVarManager::kKinkIndex0,0.);	
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaEle,-2.,3.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPio,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kTPCnSigmaPro,3.5,1000.);
    pt->AddCut(AliDielectronVarManager::kNclsTPC,70.,160.);
    pt->AddCut(AliDielectronVarManager::kNclsITS,4.,6.);
    pt->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
    pt->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
    pt->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  }

  cuts->AddCut(pt);
  
}  
  
  //exclude conversion electrons selected by the tender
  //   AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  //   noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //   cuts->AddCut(noconv);
  //	cuts->Print();
  //}

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
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);
  gammaV0Cuts->SetExcludeTracks(kTRUE);//ktrue excludes tracks v0s, 
  die->GetTrackFilter().AddCuts(gammaV0Cuts);
  //gammaV0Cuts->Print();	// 
}

void SetupV0add(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  AliDielectronV0Cuts *gammaV0Add = new AliDielectronV0Cuts("IsGamma2","IsGamma2");
  gammaV0Add->SetPdgCodes(22,11,11);
  gammaV0Add->SetDefaultPID(13); // TPC+-10, TOF+-3 if available
  gammaV0Add->AddCut(AliDielectronVarManager::kCosPointingAngle,TMath::Cos(0.02),1.0, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kChi2NDF, 0.0,  10.0, kFALSE);//to be checked, if properly filled
  gammaV0Add->AddCut(AliDielectronVarManager::kLegDist, 0.0,   0.25, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kR, 3.0,  90.0, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kPsiPair, 0.0,   0.05, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kM, 0.0,   0.05, kFALSE);
  gammaV0Add->AddCut(AliDielectronVarManager::kArmPt, 0.0,   0.05, kFALSE);	
  gammaV0Add->SetExcludeTracks(kFALSE);//ktrue excludes tracks v0s, 

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
void InitHistograms(AliDielectron *die, Int_t cutDefinition){
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
  
    //Pair classes
    // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));  

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
  }
    
    //add histograms to event class
    histos->AddClass("Event");
    histos->UserHistogram("Event","","",700,0.,700.,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",100,0.,10.,AliDielectronVarManager::kPairs);
    histos->UserHistogram("Event","","",4,0.,4.,AliDielectronVarManager::kNevents);    
    histos->UserHistogram("Event","","",300,-15.,15.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","","",100,-2.,2.,100,-2.,2.,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserProfile("Event","","",AliDielectronVarManager::kNTrk, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","",AliDielectronVarManager::kPairs, GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(700,0.,700),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(10000,0.,10),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPairs);    
    //  }
  
  //add histograms to Track classes
  //run number dependence

  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(100,-5.,+5.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kITSnSigmaEle);  
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(200,0.,200.),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kLeg1DCAabsXY);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(20,-1.0,1.0),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","","",GetRunNumbers(),AliDielectronHelper::MakeLinBinning(63,0.,6.32),AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","","",200,0,2000,40,-1.0,1.0, AliDielectronVarManager::kNTrk, AliDielectronVarManager::kEta);

  //nsigma
  histos->UserHistogram("Track","","",200,0.2,8.,AliDielectronVarManager::kPt); 
  histos->UserHistogram("Track","","",200,0.,8.,AliDielectronVarManager::kP);
  histos->UserHistogram("Track","","",80,-1.0,1.0,AliDielectronVarManager::kEta); 
  histos->UserHistogram("Track","","",400,0.0,4.0,AliDielectronVarManager::kPhi);  
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,-15.,15.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE); 
  histos->UserHistogram("Track","","",
			200,0.2,8.,100,0.,1.2,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta,kTRUE);
  histos->UserHistogram("Track","","",
			200,0.2,8.,200,-20.,20.,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);   
  histos->UserHistogram("Track","","",
			200,0.2,8.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal); 
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
                        160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  histos->UserHistogram("Track","","",
			100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(150,0.,150.),AliDielectronVarManager::kNaccTrckltsEsd10Corr);
  histos->UserHistogram("Track","","",AliDielectronHelper::MakeLinBinning(5,-0.5,2.0),AliDielectronVarManager::kTOFPIDBit);
  //rjim findable cluster vs pt
  histos->UserHistogram("Track","","",200,0.,8.,161,-0.5,161.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  //rjim frac find vs pt
  histos->UserHistogram("Track","","",			
			200,0.0,8.0,160,0,1.1,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCrFrac);
  // histos jpsi 2nd part
  //check tofbit
  //  histos->UserHistogram("Track","","",3.,-0.5,2.5, 200,0.,20.,AliDielectronVarManager::kTOFPIDBit,AliDielectronVarManager::kPt);

  //jpsi leg  nsigma vs centrality 
  histos->UserHistogram("Track","","",
			12,40.,100., 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kITSnSigmaEle);

  histos->UserHistogram("Track","","",
			12,40.,100., 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCnSigmaEle);

  histos->UserHistogram("Track","","",
			12,40.,100., 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTOFnSigmaEle);

  //Jpsi leg  nsigma(TPC) vs centrality vs Pin (3D)
  histos->UserHistogram("Track","","",
			12,40.,100., 200,0.2,8., 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kTPCnSigmaEle);

  //Jpsi leg  nsigma(TPC) vs centrality vs Eta (3D)
  histos->UserHistogram("Track","","",
			12,40.,100., 40,-1.0,1.0, 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kEta,  AliDielectronVarManager::kTPCnSigmaEle);

  //Jpsi leg  nsigma(TPC) vs centrality vs Phi (3D)
  histos->UserHistogram("Track","","",
			12,40.,100., 200,0.2,8., 100,-10,10,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPhi,  AliDielectronVarManager::kTPCnSigmaEle);

  
  //Jpsi leg  de/dx(TPC) vs centrality vs Pin (3D)
  histos->UserHistogram("Track","","",
			12,40.,100., 40,-1.0,1.0,200,0.,200.,AliDielectronVarManager::kCentrality, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCsignal);
  
  

  //jpsi nsigma vs eta
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
  
  histos->UserHistogram("Track","","",
			100,0.0,100.,100,-10,10,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
			100,0.0,100.,200,-20,20,AliDielectronVarManager::kNaccTrckltsEsd10Corr,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);      

  histos->UserHistogram("Track","","",
                        500,0.0,5000.,100,-10,10,AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","","",
                        500,0.0,5000.,100,-10,10,AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);

//
  
  histos->UserHistogram("Track","","",10000,0.0,1.0,AliDielectronVarManager::kM);
  //  histos->UserHistogram("Track","","",
  //	200,0.0,20.0,10000,0.0,1.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kM);	
  histos->UserHistogram("Track","","",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","","",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
   
  //add histograms to Pair classes
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","","",
			  100,-1.,1.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",
			  100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",
			  301,-.01,6.01,500,0.,5.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

    /*
    //add MC signal histograms to track class
    if(die->GetMCSignals()) {

      TString className="";
      for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
	TString sigMCname = die->GetMCSignals()->At(i)->GetName();

	// mc truth
	histos->AddClass(Form("Pair_%s_MCtruth",      sigMCname.Data()));
	histos->AddClass(Form("Track_Legs_%s_MCtruth",sigMCname.Data()));
	// mc reconstructed
	histos->AddClass(Form("Pair_%s",              sigMCname.Data()));
	histos->AddClass(Form("Track_Legs_%s",        sigMCname.Data()));
	// tracks
	histos->AddClass(Form("Track_%s_%s",AliDielectron::PairClassName(AliDielectron::kEv1PM),sigMCname.Data()));
	histos->AddClass(Form("Track_%s_%s_MCtruth",AliDielectron::PairClassName(AliDielectron::kEv1PM),sigMCname.Data()));
      } //end: loop signals
    } //end: has signals
    */
    
    // add single electron histograms
    //AddHistsEleEff(die);

  die->SetHistogramManager(histos);
  
}
//______________________________________________________________________________________
void InitCF(AliDielectron* die, Int_t cutDefinition)
{
	//
	// Setupd the CF Manager if needed
	//
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  //Event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,"40.,50.,60.,70.,80.,90.,100.");
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 2.0, 3.0, 5.0");//, 7.0, 10.0
  
  cf->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
  cf->AddVariable(AliDielectronVarManager::kY,"-1,-0.9,-0.8,-0.3,0,0.3,0.9,1.0");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-4.0,-3.0,-2.5,-2.25, -2.0,-1.75,-1.5, -1.25, -1.0, 0.0, 1.0, 3.0,5.0");
  
  // cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"4.5, 5.0, 5.5, 6.0, 1000.");
  // cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"4.5, 5.0, 5.5, 6.0, 1000.");
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,"0.0, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 2.0, 5.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,"-1.,-0.9,-0.8,0,0.8,0.9,1.0",kTRUE);

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  if (hasMC){
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeGrandMother,10000,-5000.5,4999.5,kTRUE);

    //only steps for efficiencies
    //    cf->SetStepsForMCtruthOnly();
  }

  //only in this case write MC truth info
  //  if (cutDefinition==0){
  //      cf->SetStepForMCtruth();
      //      histos->UserHistogram("Track","PdgCodeMother",";mother PDG code;#tracks",10000,-5000.5,4999.5,AliDielectronVarManager::kPdgCodeMother);
      // histos->UserHistogram("Track","PdgCode",";tracks PDG code;#tracks",10000,-5000.5,4999.5,AliDielectronVarManager::kPdgCode);                        

      // }
    
  //   cf->SetStepsForSignal();
	
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
  /*  
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
  */
}


//####### use functions to correct pid ###############################
void SetEtaCorrection(AliDielectron *die) {
  Bool_t hasMC=die->GetHasMC();
  Bool_t hasTuneOnData=kFALSE;

  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv DATA vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // either data or MC with tune on data option
  if( !hasMC  ) {
    // 2-dimensional eta correction for the centroid of electron sigmas
    //Fit centroid  pol5 for eta and linear for centrality
    fCntrdCorr = new TF2("fCntrdCorr", 
    			 "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*x",
    			 40.0, 90.0, -0.9, +0.9);
    fCntrdCorr->SetParameters(0.402308+0.213, +0.0662174 , -1.56518, -0.451915, +3.54694, +0.442897,  +0.00161659); // 0.213 is the offset to the centroid correction
        
    // 2-dimensional eta correction for the width of electron sigmas
    //Fit width  pol6 for eta and linear for centrality
    fWdthCorr = new TF2("fWdthCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6)  + [7]*x", 40.0, 90.0, -0.9, +0.9);
    fWdthCorr->SetParameters(1.03565, -0.0225206, -0.678966, +0.0813398, +1.59697, -0.0679809, -1.13399, -0.00121459); // pol6*linear

    die->SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kEta);
    die->SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kEta);

    //printf(" DATA PID correction loaded!!!\n");
  }


}
