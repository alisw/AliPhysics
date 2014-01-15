void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

//void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron*die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron*die, Int_t cutDefinition);


Bool_t kRot = 0;
Bool_t kMix = 1;

//start systematics
TString names= ("Phiv;Pt10;Open;Mass");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Config_shin_pPbTRD(Int_t cutDefinition=1)
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
  
  
  if(kRot){
    AliDielectronTrackRotator *rot = new AliDielectronTrackRotator;
    rot->SetConeAnglePhi(TMath::Pi());
    rot->SetIterations(10);
    die->SetTrackRotator(rot);
  }//kRot
  
  
  if(kMix){
    AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
    mix->SetMixType(AliDielectronMixingHandler::kAll);
    mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
    mix->SetDepth(10);
    die->SetMixingHandler(mix);
  }//kMix
  
  
  // set track cuts
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  InitHistograms(die,cutDefinition);
  //  if(cutDefinition==0){
  //  InitCF(die,cutDefinition);
  //}

  return die;
  
}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //
  //pairing with TLorentzVector
  // die->SetUseKF(kFALSE);
  //options
  die->SetPreFilterAllSigns();
  
  

  AliDielectronVarCuts *varkinCuts = new AliDielectronVarCuts("trackkine","trackkine");
  varkinCuts->AddCut(AliDielectronVarManager::kPt,0.4,100.);
  if(cutDefinition==2){
    varkinCuts->AddCut(AliDielectronVarManager::kPt,1.,100.);
  }
  //  if(cutDefinition==3){
  //  varkinCuts->AddCut(AliDielectronVarManager::kPt,1.,100.);
  // }
  varkinCuts->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
  die->GetTrackFilter().AddCuts(varkinCuts);
  
  AliDielectronTrackCuts *TrackCuts = new AliDielectronTrackCuts("StandardCuts","StandardCut");
  //refit
  TrackCuts->SetRequireTPCRefit(kTRUE);
  TrackCuts->SetRequireITSRefit(kTRUE);
  //SPD require
  // if(cutDefinition == 0 )
  TrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  die->GetTrackFilter().AddCuts(TrackCuts);
  

  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("trackkineandTPCQ","trackkine_and_TPC");
  varCuts->AddCut(AliDielectronVarManager::kPt,0.2,100.);
  varCuts->AddCut(AliDielectronVarManager::kEta,-0.8,0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.,3.);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.,1.);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,0.);
  
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,3.,7.);
  
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,80.,500.);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,100.,500.); // or NFclsTPCr?
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross,0.5,500.);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.,4.);
  die->GetTrackFilter().AddCuts(varCuts);
    
  //  if(cutDefinition > 2){
    AliDielectronPID *pid = new AliDielectronPID("pid","pid"); 
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.2,100.,kFALSE);
      pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.,0.2,100.,kTRUE);
      pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.4, 5., kFALSE);
      die->GetTrackFilter().AddCuts(pid);
      // }
  
  // }
      //  if(cutDefinition>0){
    
    AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("V0","V0");
    gammaV0Cuts->SetPdgCodes(22,11,11);
    gammaV0Cuts->SetDefaultPID(16);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
    //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
    gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
    //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
    gammaV0Cuts->SetExcludeTracks(kTRUE);
    gammaV0Cuts->Print();
    
    //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
    //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;
    
    // if(cuts)
    //   ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
    // else
    die->GetTrackFilter().AddCuts(gammaV0Cuts);
    
    
    //  }
    
}
//----------------------------------Pair Cut-------------------------------------------
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //  if(cutDefinition > ){
  AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
  PhiV->AddCut(AliDielectronVarManager::kM, 0. , 0.05);
  PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.5 , 3.2);
  die->GetPairPreFilter().AddCuts(PhiV);
  // }
  if(cutDefinition==2){
    AliDielectronVarCuts *Open = new AliDielectronVarCuts("Open","Open");//mass and Phiv together
    Open->AddCut(AliDielectronVarManager::kM, 0. , 0.1);
    Open->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.1);
    die->GetPairPreFilter().AddCuts(Open);
    
  }
  if(cutDefinition==3){
    AliDielectronVarCuts *Mcut = new AliDielectronVarCuts("Mcut","Mcut");//mass and Phiv together
    Mcut->AddCut(AliDielectronVarManager::kM, 0. , 0.1);
    //    Open->AddCut(AliDielectronVarManager::kOpeningAngle, 0. , 0.1);
    die->GetPairPreFilter().AddCuts(Mcut);
    
  }
  
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram classes
  AliDielectronHistos *histos=
    new AliDielectronHistos(die->GetName(),
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
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }

  //  histos->AddClass(Form("Pair_%s","lowMassDiele"));

  if(kMix){
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3))); //ME ++
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));//ME -+
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));//ME +-
  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7))); // ME --
  }

  if(kRot)histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));//Rot


  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",480,-12.,12.,AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","nESDTracks","ESD tracks;ESD tracks;Number events",100,0,200,AliDielectronVarManager::kTracks);
  histos->UserHistogram("Event","Nacc","Number of accepted tracks;Number events",100,0,200,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;NVtx;Number events",100,0,100,AliDielectronVarManager::kNVtxContrib);


  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",1000,0.,10.,AliDielectronVarManager::kPt);
  // histos->UserHistogram("Track","Pt_bin2","Pt;Pt [GeV];#tracks",
  //		"0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.,5.,6.,8.,10"
  //		,AliDielectronVarManager::kPt);
  //  histos->UserHistogram("Track","Eta","Eta; Eta ;#tracks",100,-1.,1.,AliDielectronVarManager::kEta);
  //  histos->UserHistogram("Track","Phi","Phi; Phi ;#tracks",640,0.,6.4,AliDielectronVarManager::kPhi);

  //  histos->UserHistogram("Track","Pt_Eta_Phi","Pt",100,0.,10.,100,-1,1,320,0,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  //histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  //histos->UserHistogram("Track","nClsoverfindablecluster","Number of found Clusters TPC over findably ;TPC number cluster over findable;#tracks",160,0.0,1.1,AliDielectronVarManager::kNFclsTPCrFrac);
  //histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",100,-1.,1.,320,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_pt","Eta vs Pt;Eta;Pt",100,-1.,1.,500,0.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Phi_pt","Eta vs Pt;Phi;Pt",640,0.,6.4,500,0.,10.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta_Phi","Eta vs Phi;Phi;Pt",640,0.,6.4,200,-1.,1.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle",500,0.,5.,700,-30,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal", 1000,0,10.,800,0,200,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
   histos->UserHistogram("Track","TPCnSigma_P","TPC number of sigmas Electrons vs Eta;Eta;TPCsigmaEle",100,0.,10.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
   //histos->UserHistogram("Track","TPCnSigma_Eta","TPC number of sigmas Electrons vs Eta;Eta;TPCsigmaEle",800,-0.8.,0.8.,800,-40,40,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
   //histos->UserHistogram("Track","TPCdEdx_Eta","TPC dEdx vs Eta;Eta;TPCsigmaEle",800,-0.8,0.8,800,0,200,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  //histos->UserHistogram("Track","TPCnSigma_Phi","TPC number of sigmas Electrons vs Phi;Phi;TPCsigmaEle",640,0.,6.4,800,-40,40,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","TPCdEdx_Phi","TPC dEdx vs Phi;Phi;TPCsigmaEle",640,0.,6.4,800,0,200,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal);

  

  //add histograms to Pair classes

   //  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
   //                     500,0,4,AliDielectronVarManager::kM);

  // histos->UserHistogram("Pair","InvMass_low","Inv.Mass;Inv. Mass [GeV];#pairs",
  //                     500,0,0.5,AliDielectronVarManager::kM);
  /*
  histos->UserHistogram("Pair","kDeltaEta","kDeltaEta;kDeltaEta;#pairs",
                        160,0.,1.6,AliDielectronVarManager::kDeltaEta);
  histos->UserHistogram("Pair","kDeltaPhi","kDeltaPhi;kDeltaPhi;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kDeltaPhi);
  histos->UserHistogram("Pair","PhiV",";PhiV;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMass_bin2","Inv.Mass;Inv. Mass [GeV];#pairs",
	"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.75 , 4.0",AliDielectronVarManager::kM);
  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        500, 0. , 4., 100 , 0., 5. ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPt );
  */
  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        400, 0. , 4., 1000 , 0., 10. ,
			AliDielectronVarManager::kM , AliDielectronVarManager::kPt );
  histos->UserHistogram("Pair",
                        "InvMass_Eta","InvMass_Eta;InvMass;Eta",
                        400, 0. , 4., 200 , -1., 1. ,
			AliDielectronVarManager::kM , AliDielectronVarManager::kEta );
  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                        400, 0. , 4., 320 , 0., 3.2 ,
			AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );
  histos->UserHistogram("Pair",
                        "InvMass_dPhi","InvMass_PhivPair;InvMass;PhivPair",
                        400, 0. , 4., 320 , 0., 3.2 ,
			AliDielectronVarManager::kM , AliDielectronVarManager::kDeltaPhi );
  histos->UserHistogram("Pair",
                        "InvMass_OpenAng","OpAngle_InvMass;Opening angle;Invariant Mass",
                         400 , 0. , 4. , 320, 0. , 3.2,
			AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);


  /*
  histos->UserHistogram("Pair",
                        "Y_PtPair","InvMass_PhivPair;InvMass;PhivPair",
                        120, -1.2. , 1.2, 100 , 0., 5. ,
			AliDielectronVarManager::kY , AliDielectronVarManager::kPt );
  histos->UserHistogram("Pair",
                        "OpAngle","Opening angle;Opening angle;#pairs",
                        320, 0. , 3.2, 
                         AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair",
                        "OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",
                        320, 0. , 3.2, 500 , 0. , 4. ,
			AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair",
                        "Phi","Phi;counts;Phi",
                        320, 0. , 6.4, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair",
                        "Y","Y;counts;Y",
                        120, -1.2 , 1.2, AliDielectronVarManager::kY);
  histos->UserHistogram("Pair",
			"Pt","Pt;counts;Pt",
			500, 0 , 10, AliDielectronVarManager::kPt);
  */ 
 die->SetHistogramManager(histos);
  
}




void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  //  cf->AddVariable(AliDielectronVarManager::kM,500,0.,4.);
  cf->AddVariable(AliDielectronVarManager::kM,"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.75 , 4.0");//data
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,0.1);
  //  cf->AddVariable(AliDielectronVarManager::kM,500,0.,0.1);
  cf->AddVariable(AliDielectronVarManager::kY,20,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kR,500,0.,100.);
  cf->AddVariable(AliDielectronVarManager::kPhi,32, 0., 3.2);
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kPt,"0.,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,1.75,2.0,3.,5.");
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,62,0.,6.2);
  cf->AddVariable(AliDielectronVarManager::kPhivPair,64, 0., 6.4);
  cf->AddVariable(AliDielectronVarManager::kPairType,4,-0.5,3.5);


  cf->AddVariable(AliDielectronVarManager::kEta,40,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kCharge,40,-2.,2.,kTRUE);
  //leg 
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,36,0.,360.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,100,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,20,-3.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,"40,50,55,60,65,68,70,72,75,80,90,100,110,120",kTRUE);
  // cf->AddVariable(AliDielectronVarManager::kTPCsignal,200,0.,200.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kNFclsTPCr,200,0.,200.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kNFclsTPCfCross,200,0.,2.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,200,0.,10.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kImpactParXY,400,-.5,.5,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kImpactParZ,400,-2.,2.,kTRUE);
  //cf->AddVariable(AliDielectronVarManager::kNclsITS,10,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNFclsTPCrFrac,10,0.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);


  if(cutDefinition==0){
    cf->SetStepForMCtruth();
  }
  cf->SetStepForAfterAllCuts();
  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);


  //primary
  AliDielectronSignalMC* DielePrimary=new
	AliDielectronSignalMC("Diele Primary","low mass primary dielectron pairs");
  DielePrimary->SetLegPDGs(11,-11);
  DielePrimary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  DielePrimary->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //DielePrimary->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(DielePrimary);

  // secondary  
  AliDielectronSignalMC* secsec = new 
    AliDielectronSignalMC("Secondarypairs","secondary electron pairs");      // pairs made from conversion (may be also from 2 different conversions)
  secsec->SetLegPDGs(11,-11);
  secsec->SetCheckBothChargesLegs(kTRUE,kTRUE);
  secsec->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  //  secsec->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(secsec);
  
  // secondary-primary
  AliDielectronSignalMC* DielePriSec=new
	AliDielectronSignalMC("Diele Prim-Sec","low mass prim-sec dielectron pairs");
  DielePriSec->SetLegPDGs(11,-11);
  DielePriSec->SetCheckBothChargesLegs(kTRUE,kTRUE);
  DielePriSec->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  //DielePriSec->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(DielePriSec);

  /*
  AliDielectronSignalMC* DielePrimary=new
	AliDielectronSignalMC("Diele Primary","low mass primary dielectron pairs");
  DielePrimary->SetLegPDGs(11,-11);
  DielePrimary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  DielePrimary->SetLegSources(AliDielectronSignalMC::kFinalState,
	  AliDielectronSignalMC::kFinalState);
  DielePrimary->SetMothersRelation(AliDielectronSignalMC::kSame);
  DielePrimary->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(DielePrimary);
  */
    
  AliDielectronSignalMC* pi0 = new AliDielectronSignalMC("pi0dalitz","pi0dalitz");
  pi0->SetLegPDGs(11,-11);
  pi0->SetMotherPDGs(111,111);
  pi0->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0->SetLegSources(AliDielectronSignalMC::kFinalState,AliDielectronSignalMC::kFinalState);
  pi0->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0);

   AliDielectronSignalMC* etaSig = new AliDielectronSignalMC("Eta", "etaSignal"); ///eta dalitz pairs 
  etaSig->SetLegPDGs(11,-11);
  etaSig->SetMotherPDGs(221,221);
  etaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaSig);

  AliDielectronSignalMC* etaprimeSig = new AliDielectronSignalMC("Etaprime", "etaprimeSignal"); ///etaprime pairs 
  etaprimeSig->SetLegPDGs(11,-11);
  etaprimeSig->SetMotherPDGs(331,331);
  etaprimeSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  etaprimeSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  etaprimeSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  etaprimeSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaprimeSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  etaprimeSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(etaprimeSig);

  AliDielectronSignalMC* rhoSig = new AliDielectronSignalMC("Rho", "rhoSignal"); ///rho pairs 
  rhoSig->SetLegPDGs(11,-11);
  rhoSig->SetMotherPDGs(113,113);
  rhoSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  rhoSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  rhoSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  rhoSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  rhoSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  rhoSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(rhoSig);

  AliDielectronSignalMC* omegaSig = new AliDielectronSignalMC("Omega", "omegaSignal"); ///omega pairs 
  omegaSig->SetLegPDGs(11,-11);
  omegaSig->SetMotherPDGs(223,223);
  omegaSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  omegaSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  omegaSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  omegaSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  omegaSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  omegaSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(omegaSig);
  
  AliDielectronSignalMC* phiSig = new AliDielectronSignalMC("Phi", "phiSignal"); ///phi pairs 
  phiSig->SetLegPDGs(11,-11);
  phiSig->SetMotherPDGs(333,333);
  phiSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  phiSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  phiSig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  phiSig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  phiSig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  phiSig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(phiSig);

  AliDielectronSignalMC* dieleOpenCharm = new AliDielectronSignalMC("dieleOpenCharm", "dieleOpenCharm");
  dieleOpenCharm->SetLegPDGs(11,-11);
  dieleOpenCharm->SetMotherPDGs(402,402);
  dieleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState,AliDielectronSignalMC::kFinalState);
  dieleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  dieleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  dieleOpenCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(dieleOpenCharm);
  
  AliDielectronSignalMC* diEleCharm = new AliDielectronSignalMC("diEleCharm","di-electrons from charm");  // dielectrons originating from charm hadrons (not neccessary from same mother)
  diEleCharm->SetLegPDGs(11,-11);
  diEleCharm->SetMotherPDGs(403,403);
  diEleCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleCharm);
 
  
}





