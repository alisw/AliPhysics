
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);


TString names= ("pt400;pt400DCA;V0gamma;PhiV;strPID");

	Bool_t kRot = 0;
	Bool_t kMix = 1;

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Config_lowmass(Int_t cutDefinition=1)
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
 SetupCuts(die,cutDefinition);



  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //

  InitHistograms(die,cutDefinition);
  InitCF(die,cutDefinition);

  
  return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //



	//options
        die->SetPreFilterAllSigns();

	//pairing with TLorentzVector
	die->SetUseKF(kFALSE);

	//tracks only
//	if(cutDefinition == 0)die->SetNoPairing(); 

	//track cuts
        die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
	//pid cuts
	die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));


           if(cutDefinition >= 0){
                AliDielectronVarCuts *RapidityCut=new AliDielectronVarCuts("RapidityCut","RapidityCut");
                RapidityCut->AddCut(AliDielectronVarManager::kY, -0.9 , 0.9);
                die->GetPairFilter().AddCuts(RapidityCut);
           }


	   if(cutDefinition==2){//try the V0 finder
	     AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
	     noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
	     die->GetTrackFilter().AddCuts(noconv);
	   }

         if(cutDefinition == 3){
         AliDielectronVarCuts *phiVcut=new AliDielectronVarCuts("phiVcut","Phi_V cut <2.5 rad");
	 phiVcut->AddCut(AliDielectronVarManager::kPhivPair, 2.5, TMath::ACos(-1.0), kTRUE);
	 die->GetPairPreFilter().AddCuts(phiVcut);
           }


}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){

  AliDielectronPID *pid = new AliDielectronPID();


  if(cutDefinition >= 0){//400
        //TPC
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.4,100.,kFALSE);
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.,0.4,100.,kTRUE);
        //TOF
        pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.4, 3., kFALSE);
        }


  if(cutDefinition == 4){//str PID
        //TPC
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,0.,3.,0.4,100.,kFALSE);
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.,0.4,100.,kTRUE);
        //TOF
        pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.4, 3., kFALSE);
        }


 return pid;

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{

  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *fesdTrackCuts1 = new AliESDtrackCuts;

  //global

  //  if(cutDefinition == 0)fesdTrackCuts1->SetPtRange( 0.2 , 100. );
  fesdTrackCuts1->SetPtRange( 0.4 , 100. );
  fesdTrackCuts1->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts1->SetDCAToVertex2D(kTRUE);
  if(cutDefinition == 1){
  fesdTrackCuts1->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts1->SetMaxDCAToVertexXY(1.);
  }else{
  fesdTrackCuts1->SetMaxDCAToVertexZ(2.5);
  fesdTrackCuts1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  }
  fesdTrackCuts1->SetEtaRange( -0.8 , 0.8 );

//  fesdTrackCuts1->SetMaxNsigmaToVertex( 2. );
//  fesdTrackCuts1->SetMaxCovDiagonalElements( 2 , 2, 0.5, 0.5 , 2 );

  //ITS
  fesdTrackCuts1->SetRequireITSRefit(kTRUE);
  fesdTrackCuts1->SetMinNClustersITS(3); 
  //fesdTrackCuts1->SetMaxChi2PerClusterITS(4);
  fesdTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  //kAny kBoth kFirst kSecond

  //TPC
  fesdTrackCuts1->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts1->SetMinNClustersTPC(80);
  fesdTrackCuts1->SetMinNCrossedRowsTPC(100);
  fesdTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
  fesdTrackCuts1->SetMaxChi2PerClusterTPC(4);
  fesdTrackCuts1->SetMaxFractionSharedTPCClusters(0.5); 

  return fesdTrackCuts1;

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
  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta","Eta; Eta ;#tracks",100,-1.,1.,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi; Phi ;#tracks",640,0.,6.4,AliDielectronVarManager::kPhi);


  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",100,-1.,1.,320,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_pt","Eta vs Pt;Eta;Pt",100,-1.,1.,500,0.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","ImpParXY_Pt","ImpParXY_Pt; Pt; ÍmpParXY",500,0.,10.,500,-5.,5.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ_Pt","ImpParZ_Pt; Pt; ÍmpParZ",500,0.,10.,500,-5.,5.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParZ);


  //track checks (ITS)
  histos->UserHistogram("Track","ITSchi2Cl_Mom","ITS Chi2 vs Momentum;Mom;ITS chi2",500,0.,5.,50,0.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsITS_Mom",";Mom;kNclsITS",500,0.,5.,7,0,7,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsITS);

  //track checks (TPC)
  histos->UserHistogram("Track","TPCsignalNfrac_Mom",";fraction TPCSignalN/TPCncls vs Momentum;Mom;TPCSignalN/TPCncls",500,0.,5.,60,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignalNfrac);
  histos->UserHistogram("Track","TPCchi2Cl_Mom","TPC Chi2 vs Momentum;Mom;TPC Chi2",500,0.,10.,100,0,5,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","TPCclsDiff_Mom","kTPCclsDiff vs Momentum;Mom;kTPCclsDiff",500,0.,10.,100,-10,10,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","FindableTPCcls_Mom","kNFclsTPC vs Momentum;Mom;kNFclsTPC",500,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNFclsTPC);
  histos->UserHistogram("Track","TPCcls_Mom","kNclsTPC vs Momentum;Mom;kNclsTPC",500,0.,10.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","kNclsSFracTPC_Mom","kNclsSFracTPC vs Momentum;Mom;kTPCclsSFrac",500,0.,10.,60,0.,0.12,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","kNFclsTPCrFrac_Mom","kNFclsTPCrFrac vs Momentum;Mom;kNFclsTPCrFrac",500,0.,10.,60,0.,1.2.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNFclsTPCrFrac);

  //track checks (TOF)
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta",500,0.,5.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFPIDBit_Mom","kTOFPIDBit vs Momentum;Mom;TOFPIDbit",500,0.,5.,2,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFPIDBit);

  //track checks (TRD)
  histos->UserHistogram("Track","NclsTRD_Mom","kNclsTRD vs Momentum;Mom;NclsTRD",500,0.,5.,20,0.,20,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsTRD);

  histos->UserHistogram("Track","TRDntracklets_Mom","TRDntracklets vs Momentum;Mom;TRDnTracklets",500,0.,5.,20,0.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDntracklets);

  histos->UserHistogram("Track","TRDpidProb_Electrons_Mom","kTRDpidProb_Electrons vs Momentum;Mom;TRDpidProb_Electrons",500,0.,5.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle);

   histos->UserHistogram("Track","kTRDpidProb_Pions_Mom","kTRDpidProb_Pions vs Momentum;Mom;kTRDpidProb_Pions",500,0.,5.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio);

  //pid checks
  histos->UserHistogram("Track","ITSnSigma_MomPio","ITS number of sigmas Pion vs Momentum;Mom;ITSsigmaPion",500,0.,5.,1000,-20,20,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaPio);
  histos->UserHistogram("Track","TPCnSigma_MomPio","TPC number of sigmas Pions vs Momentum;Mom;TPCsigmaPion",500,0.,5.,1000,-20,20,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);

  histos->UserHistogram("Track","ITSnSigma_MomEle","ITS number of sigmas Electrons vs Momentum;Mom;ITSsigmaEle",500,0.,5.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle",1000,0.,10.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TOFnSigma_MomEle","TOF number of sigmas Electrons vs Momentum;Mom;TOFsigmaEle",500,0.,5.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);

  histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal", 500,0,5,800,0,200,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal", 500,0,10.,800,0,200,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);

  //
  //add histograms to Pair classes
  //

  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0,4,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","PhiV",";PhiV;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kPhivPair);


   histos->UserHistogram("Pair","PhiV_Pt",";Pt;PhiV",
			 100,0.,10.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);

   histos->UserHistogram("Pair","InvMass_bin2","Inv.Mass;Inv. Mass [GeV];#pairs",
	"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.75 , 4.0",AliDielectronVarManager::kM);


   histos->UserHistogram("Pair","InvMass_bin3","Inv.Mass;Inv. Mass [GeV];#pairs",
	"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.65 ,3.8, 4.0",AliDielectronVarManager::kM);

   histos->UserHistogram("Pair","InvMass_bin4","Inv.Mass;Inv. Mass [GeV];#pairs",
	"0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
 	 0.5 , 0.6, 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
	 0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
	2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
	3.3 , 3.5, 3.65 ,3.8, 4.0",AliDielectronVarManager::kM);

  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        500, 0. , 4., 100 , 0., 10. ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPt );

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
                        -1.2, 1.2 , 2.4, AliDielectronVarManager::kY);

  die->SetHistogramManager(histos);
}




void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());


   //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.5, 3., 3.5, 4.0, 4.5, 5.");
  //  cf->AddVariable(AliDielectronVarManager::kP,"0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.5, 3., 3.5, 4.0, 4.5, 5.");
  cf->AddVariable(AliDielectronVarManager::kM,100,0.,4.);
  cf->AddVariable(AliDielectronVarManager::kY,24,-1.2,1.2);
  //  cf->AddVariable(AliDielectronVarManager::kPhi,36, 0., 360.);

  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,62,0.,6.2);
  cf->AddVariable(AliDielectronVarManager::kPhivPair,64, 0., 6.4);
  cf->AddVariable(AliDielectronVarManager::kPairType,4,-0.5,3.5);
  //leg 
  cf->AddVariable(AliDielectronVarManager::kPt,"0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.5, 3., 3.5, 4.0, 4.5, 5.",kTRUE);
    cf->AddVariable(AliDielectronVarManager::kP,"0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.5, 3., 3.5, 4.0, 4.5, 5.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,20,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,36,0.,360.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kY,50,-1.,1.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,40,-2.,2.,kTRUE);

  //  cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,50,-10.,10.,kTRUE);
  //  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"0, 70, 75, 80, 85, 90, 100, 120, 160",kTRUE);

  die->SetCFManagerPair(cf);
  
}
