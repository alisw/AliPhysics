
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);


//start systematics
TString names= ("Base;SPDboth;TPCshift;TOFtight;PhiVhard");


	Bool_t kRot = 0;
	Bool_t kMix = 1;

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Config_lowmasspp(Int_t cutDefinition=1)
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


  // eta correction
//  SetEtaCorrection();

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

	//track cuts
        die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
	//pid cuts
	die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));

  
	if(cutDefinition >= 0){
             AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
             noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
             die->GetTrackFilter().AddCuts(noconv);
	     
	     AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
             PhiV->AddCut(AliDielectronVarManager::kM, 0. , 0.15);
             PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2. , 3.2);
             die->GetPairPreFilter().AddCuts(PhiV);

	     
	}
		
		
	if(cutDefinition == 4){
	
	     AliDielectronVarCuts *PhiVhard = new AliDielectronVarCuts("PhiVhard","PhiVhard");//mass and PhiVhard together (hard PhiVhard)
             PhiVhard->AddCut(AliDielectronVarManager::kM, 0. , 0.05);
             PhiVhard->AddCut(AliDielectronVarManager::kPhivPair, 1.5 , 3.2);
             die->GetPairPreFilter().AddCuts(PhiVhard);
	}
		



}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){

  AliDielectronPID *pid = new AliDielectronPID();

 if(cutDefinition == 0 || cutDefinition == 1 || cutDefinition == 4){//TPC+TOF pT200
   //stronger electron inclusion
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.2,100.,kFALSE);
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.,0.2,100.,kTRUE);
        //TOF
        pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.4, 5., kFALSE);
 }//complete PID


 if(cutDefinition == 2){//TPC+TOF pT200 shift in sigma_e (TPC)
   //stronger electron inclusion
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-2.,2.5,0.2,100.,kFALSE);
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.,0.2,100.,kTRUE);
        //TOF
        pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3.,3.,0.4, 5., kFALSE);
 }//complete PID
 

 if(cutDefinition == 3){//TPC+TOF pT200 tighter sigma_e (TOF)
   //stronger electron inclusion
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-1.5,3.,0.2,100.,kFALSE);
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.,0.2,100.,kTRUE);
        //TOF
        pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-2.,2.,0.4, 5., kFALSE);
 }//complete PID



 return pid;

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;


  //global
  
  fesdTrackCuts->SetPtRange( 0.2 , 100. );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );

  //ITS
  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersITS(3);
  if(cutDefinition == 1)fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth); 
  else fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);


  //TPC
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersTPC(80);
  fesdTrackCuts->SetMinNCrossedRowsTPC(100);
  fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
  fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  return fesdTrackCuts;


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
  histos->UserHistogram("Track","kNclsSFracTPC_Mom","kNclsSFracTPC vs Momentum;Mom;kTPCclsSFrac",500,0.,10.,1000,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","kNFclsTPCrFrac_Mom","kNFclsTPCrFrac vs Momentum;Mom;kNFclsTPCrFrac",500,0.,10.,60,0.,1.2.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kNFclsTPCrFrac);

  //track checks (TOF)
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta",500,0.,5.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFPIDBit_Mom","kTOFPIDBit vs Momentum;Mom;TOFPIDbit",500,0.,5.,2,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFPIDBit);

  //track checks (PID)
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

  histos->UserHistogram("Pair","InvMass_low","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0,0.5,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","InvMass10","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0.,5.,AliDielectronVarManager::kM);

  histos->UserHistogram("Pair","kDeltaEta","kDeltaEta;kDeltaEta;#pairs",
                        160,0.,1.6,AliDielectronVarManager::kDeltaEta);

  histos->UserHistogram("Pair","kDeltaEta_low","kDeltaEta;kDeltaEta;#pairs",
                        500,0.,0.5,AliDielectronVarManager::kDeltaEta);

  histos->UserHistogram("Pair","kDeltaPhi","kDeltaPhi;kDeltaPhi;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kDeltaPhi);

  histos->UserHistogram("Pair",
                        "kDeltaEta_kDeltaPhi","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        160, 0. , 1.6, 320 , 0., 6.4 ,
                         AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );

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


  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        500, 0. , 4., 100 , 0., 5. ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPt );

  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                        500, 0. , 4., 320 , 0., 3.2 ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );

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
  cf->AddVariable(AliDielectronVarManager::kY,20,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kPhi,32, 0., 3.2);
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,5.);


  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,62,0.,6.2);
  cf->AddVariable(AliDielectronVarManager::kPhivPair,64, 0., 6.4);
  cf->AddVariable(AliDielectronVarManager::kPairType,4,-0.5,3.5);
  
  
  
  //leg 
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,36,0.,360.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,100,-1.,1.,kTRUE);

  
/*  
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,200,0.,200.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kImpactParXY,40,-2.,2.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kImpactParZ,40,-2.,2.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsITS,10,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNFclsTPCrFrac,10,0.,1.,kTRUE);
*/


  die->SetCFManagerPair(cf);
  
}





void SetEtaCorrection()
{
  if (AliDielectronPID::GetEtaCorrFunction()) return;

  TString list=gSystem->Getenv("LIST");
  
  TFile f("$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root");
  if (!f.IsOpen()) return;
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(list)){
      printf("Using Eta Correction Function: %s\n",kName.Data());
      AliDielectronPID::SetEtaCorrFunction((TF1*)f.Get(kName.Data()));
    }
  }
}



