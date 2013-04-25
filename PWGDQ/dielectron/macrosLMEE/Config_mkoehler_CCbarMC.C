
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);

void AddMCsignal(AliDielectron *die);

void SetEtaCorrection();

//TString names= ("Electrons;QA;EMCal_AnyLeg");
TString names= ("Electrons;QA");

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

AliDielectron* Config_mkoehler_CCbarMC(Int_t cutDefinition=1,Bool_t isAOD)
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

  // cut setup
  SetupTrackCuts(die,cutDefinition);



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
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //




        AddMCsignal(die);

	//track cuts
        if(cutDefinition>0){
        die->GetTrackFilter().AddCuts(SetupESDtrackCuts(cutDefinition));
	}

	//pid cuts
//	die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));


       
	if(cutDefinition >= 0){

	  
	    AliDielectronVarCuts *mcpid =new AliDielectronVarCuts("mcpid","mcpid");
	    mcpid->SetCutType(AliDielectronVarCuts::kAny);
	    mcpid->AddCut(AliDielectronVarManager::kPdgCode, 11);
	    mcpid->AddCut(AliDielectronVarManager::kPdgCode, -11);
	    die->GetTrackFilter().AddCuts(mcpid);

	   AliDielectronVarCuts *MCnoConv =new AliDielectronVarCuts("MCnoConv","MCnoConv");
       	   MCnoConv->AddCut(AliDielectronVarManager::kPdgCodeMother, 22,kTRUE);
	   die->GetTrackFilter().AddCuts(MCnoConv);	     
	     
	}//Base



	if(cutDefinition == 2){//EMCal cut
	  AliDielectronVarCuts *mycut = new AliDielectronVarCuts("CutEMCAL","cut for EMCal");
  	  mycut->AddCut(AliDielectronVarManager::kEMCALnSigmaEle,-3.5,10.);
	  mycut->AddCut(AliDielectronVarManager::kEMCALE,3.5,100.);
	  mycut->AddCut(AliDielectronVarManager::kEMCALEoverP,0.7,1.3);

  
	  AliDielectronPairLegCuts *varpair=new AliDielectronPairLegCuts();
	  varpair->GetLeg1Filter().AddCuts(mycut);
	  varpair->GetLeg2Filter().AddCuts(mycut);
	  varpair->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
//	  die->GetPairFilter().AddCuts(varpair);
	}

}

//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){

  AliDielectronPID *pid = new AliDielectronPID();

  

 if(cutDefinition >= 1){
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.5,4.5,0.2,100.,kFALSE);
 }//complete PID



 
 return pid;

}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts;

  //global
  fesdTrackCuts->SetPtRange( 0.8 , 100. );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  //ITS
  fesdTrackCuts->SetRequireITSRefit(kTRUE);
  fesdTrackCuts->SetMinNClustersITS(3); 
  fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);//kFirst
  
  //TPC
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
//  fesdTrackCuts->SetMinNClustersTPC(80);
//  fesdTrackCuts->SetMinNCrossedRowsTPC(120);
//  fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
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
  histos->UserHistogram("Track","ImpParXY","ImpParXY; ÍmpParXY ;#tracks",500,-5.,5.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","ImpParZ","ImpParZ; ÍmpParZ ;#tracks",500,-5.,5.,AliDielectronVarManager::kImpactParZ);

  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",100,-1.,1.,320,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_pt","Eta vs Pt;Eta;Pt",100,-1.,1.,500,0.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kPt);



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

  histos->UserHistogram("Pair","Y",";Rapidity;#pairs",
                        120,-1.2,1.2,AliDielectronVarManager::kY);
			
  histos->UserHistogram("Pair","Eta",";PseudoRapidity;#pairs",
                        120,-1.2,1.2,AliDielectronVarManager::kEta);
			
  histos->UserHistogram("Pair","Phi",";PseudoRapidity;#pairs",
                        320,0,3.2,AliDielectronVarManager::kPhi);
			
  histos->UserHistogram("Pair","Eta_Phi",";PseudoRapidity;Phi",
                        120,-1.2,1.2,320,0,3.2,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
			
			
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0,4,AliDielectronVarManager::kM);


  histos->UserHistogram("Pair","InvMassBin2000","Inv.Mass;Inv. Mass [GeV];#pairs",
                        2000,0.,10.,AliDielectronVarManager::kM);

   histos->UserHistogram("Pair","InvMass_bin2","Inv.Mass;Inv. Mass [GeV];#pairs",
        "0. , 0.025, 0.05 , 0.075 ,0.1 ,0.15 , 0.2 , 0.25 , 0.3 , 
        0.4 ,  0.5 , 0.6, 0.65 , 0.688 , 0.725, 0.75, 0.775, 0.8 , 0.85 ,
         0.95,  0.975 , 1.0 , 1.025 , 1.05, 1.125 , 1.25 , 1.5 , 1.75 , 2.0 , 
        2.25, 2.5 , 2.75 , 2.85, 2.95,3.05, 3.1 , 3.15 , 
        3.3 , 3.5, 3.75 , 4.0",AliDielectronVarManager::kM);

  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                        500, 0. , 4., 320 , 0., 3.2 ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );
			 
  histos->UserHistogram("Pair","PhiV",";PhiV;#pairs",
                        320,0.,6.4,AliDielectronVarManager::kPhivPair);


   histos->UserHistogram("Pair","PhiV_Pt",";Pt;PhiV",
			 100,0.,10.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhivPair);
			 
			 
  histos->UserHistogram("Pair",
                        "InvMass_Pt","InvMass_Pt;InvMass;Pt",
                        500, 0. , 4., 100 , 0., 5. ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPt );

  histos->UserHistogram("Pair",
                        "OpAngle","Opening angle;Opening angle;#pairs",
                        320, 0. , 3.2, 
                         AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair",
                        "OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",
                        320, 0. , 3.2, 500 , 0. , 4. ,
                         AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);



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



  die->SetHistogramManager(histos);
}






void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

   
    //pair variables
   cf->AddVariable(AliDielectronVarManager::kM, 500,0.,10.);
	
   cf->AddVariable(AliDielectronVarManager::kY,24,-1.2,1.2);
   cf->AddVariable(AliDielectronVarManager::kPhi,64, -3.2, 3.2);
   cf->AddVariable(AliDielectronVarManager::kPt,100,0.,5.);
   
   cf->AddVariable(AliDielectronVarManager::kPairType,4,-0.5,3.5);
   cf->AddVariable(AliDielectronVarManager::kNumberOfDaughters,5,0,5);
   cf->AddVariable(AliDielectronVarManager::kPhivPair,640,0.,3.2);

  //leg 
  cf->AddVariable(AliDielectronVarManager::kPt,100,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,36,0.,360.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,100,-1.,1.,kTRUE);


  cf->AddVariable(AliDielectronVarManager::kImpactParXY,40,-2.,2.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kImpactParZ,40,-4.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsITS,10,0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNFclsTPCrFrac,10,0.,1.,kTRUE);


  //TPC pid
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,80,-4.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,8,1.,4.5,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,8,0.,4.,kTRUE);

  if (cutDefinition == 0){
    cf->SetStepForMCtruth();
  }
  
//     cf->SetStepsForEachCut();
//     cf->SetStepsForCutsIncreasing();


  die->SetCFManagerPair(cf);
  
}


void AddMCsignal(AliDielectron *die){

  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("diEleOpenCharm","di-electrons from open charm");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenCharm);


  AliDielectronSignalMC* diEleOpenBeauty = new AliDielectronSignalMC("diEleOpenBeauty","di-electrons from open beauty");  // dielectrons originating f$
  diEleOpenBeauty->SetLegPDGs(11,-11);
  diEleOpenBeauty->SetMotherPDGs(503,503);
  diEleOpenBeauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenBeauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBeauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBeauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenBeauty);
  //signal pairs

   AliDielectronSignalMC* pi0 = new AliDielectronSignalMC("pi0","pi0");
   pi0->SetLegPDGs(11,-11);
   pi0->SetMothersRelation(AliDielectronSignalMC::kSame);
   pi0->SetMotherPDGs(111,111);//pi0
   pi0->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   pi0->SetFillPureMCStep(kTRUE);
   pi0->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(pi0);

   AliDielectronSignalMC* eta = new AliDielectronSignalMC("eta","eta");
   eta->SetLegPDGs(11,-11);
   eta->SetMothersRelation(AliDielectronSignalMC::kSame);
   eta->SetMotherPDGs(221,221);//eta
   eta->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   //   eta->SetFillPureMCStep(kTRUE);
   eta->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(eta);

   AliDielectronSignalMC* etaprime = new AliDielectronSignalMC("etaprime","etaprime");
   etaprime->SetMothersRelation(AliDielectronSignalMC::kSame);
   etaprime->SetMotherPDGs(331,331);//etaprime
   etaprime->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   //   etaprime->SetFillPureMCStep(kTRUE);
   etaprime->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(etaprime);

   AliDielectronSignalMC* omega = new AliDielectronSignalMC("omega","omega");
   omega->SetLegPDGs(11,-11);
   omega->SetMothersRelation(AliDielectronSignalMC::kSame);
   omega->SetMotherPDGs(223,223);//omega
   omega->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   omega->SetFillPureMCStep(kTRUE);
   omega->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(omega);

   AliDielectronSignalMC* phi = new AliDielectronSignalMC("phi","phi");
   phi->SetLegPDGs(11,-11);
   phi->SetMothersRelation(AliDielectronSignalMC::kSame);
   phi->SetMotherPDGs(333,333);//phi
   phi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   phi->SetFillPureMCStep(kTRUE);
   phi->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(phi);


   AliDielectronSignalMC* jpsi = new AliDielectronSignalMC("jpsi","jpsi");
   jpsi->SetLegPDGs(11,-11);
   jpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
   jpsi->SetMotherPDGs(443,443);//jpsi
   jpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//   jpsi->SetFillPureMCStep(kTRUE);
   jpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   die->AddSignalMC(jpsi);
   
   
 AliDielectronSignalMC* diEleContinuum = new AliDielectronSignalMC("diEleContinuum","di-electron continuum");     // all di-electrons originating in the collision
  diEleContinuum->SetLegPDGs(11,-11);
  diEleContinuum->SetMotherPDGs(0,0,22,22);
  diEleContinuum->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleContinuum->SetMothersRelation(AliDielectronSignalMC::kSame);
  diEleContinuum->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleContinuum->SetFillPureMCStep(kTRUE);
//  die->AddSignalMC(diEleContinuum);

  //
  // contamination from misidentification and conversions
  //

  //conversion
  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");      // pairs made from conversion
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  //  conversionElePairs->SetMothersRelation(AliDielectronSignalMC::kSame);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
//  die->AddSignalMC(conversionElePairs);

  AliDielectronSignalMC* conversionLeg = new AliDielectronSignalMC("conversionLeg","pairs with a leg from conversion");      // pairs made from conversion
  conversionLeg->SetLegPDGs(11,-11);
  conversionLeg->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionLeg->SetLegSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kSecondary);
  conversionLeg->SetMotherPDGs(22,0);
//  die->AddSignalMC(conversionLeg);



  // background from secondary electrons
  AliDielectronSignalMC* secondaryElectrons = new AliDielectronSignalMC("secondaryElectrons","Secondary electrons");   // all di-electrons from secondary electrons (interaction with detector)
  secondaryElectrons->SetLegPDGs(11,-11);
  secondaryElectrons->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  secondaryElectrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(secondaryElectrons);

  AliDielectronSignalMC* primarySecElePairs = new AliDielectronSignalMC("primarySecElePairs","Primary+Secondary electron pairs");  // primary-secondary pairs
  primarySecElePairs->SetLegPDGs(11,-11);
  primarySecElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  primarySecElePairs->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
//  die->AddSignalMC(primarySecElePairs);

  // misidentification
  AliDielectronSignalMC* allEleMisIdPairs = new AliDielectronSignalMC("allEleMisIdPairs","all electron+misid. pairs");  // one true electron + a mis-id electron (all sources included)
  allEleMisIdPairs->SetLegPDGs(11,11,kFALSE,kTRUE);
  allEleMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(allEleMisIdPairs);

  AliDielectronSignalMC* allMisIdMisIdPairs = new AliDielectronSignalMC("allMisIdMisIdPairs","all misid.+misid. pairs");  // mis-id + mis-id
  allMisIdMisIdPairs->SetLegPDGs(11,11,kTRUE,kTRUE);
  allMisIdMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(allMisIdMisIdPairs);

  AliDielectronSignalMC* elePionPairs = new AliDielectronSignalMC("elePionPairs","electron+pion pairs");    // true electron + mis-id pion
  elePionPairs->SetLegPDGs(11,211);
  elePionPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(elePionPairs);

  AliDielectronSignalMC* eleKaonPairs = new AliDielectronSignalMC("eleKaonPairs","electron+kaon pairs");   // true electron + mis-id kaon
  eleKaonPairs->SetLegPDGs(11,321);
  eleKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(eleKaonPairs);

  AliDielectronSignalMC* eleProtonPairs = new AliDielectronSignalMC("eleProtonPairs","Electron+proton pairs");  // true electron + mis-id proton
  eleProtonPairs->SetLegPDGs(11,2212);
  eleProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(eleProtonPairs);

  AliDielectronSignalMC* piPiPairs = new AliDielectronSignalMC("piPiPairs","pion+pion pairs");    // mis-id pion + mis-id pion
  piPiPairs->SetLegPDGs(211,211);
  piPiPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(piPiPairs);

  AliDielectronSignalMC* piKaonPairs = new AliDielectronSignalMC("piKaonPairs","pion+kaon pairs");  // mis-id pion + mis-id kaon
  piKaonPairs->SetLegPDGs(211,321);
  piKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(piKaonPairs);

  AliDielectronSignalMC* piProtonPairs = new AliDielectronSignalMC("piProtonPairs","pion+proton pairs");  // mis-id pion + mis-id proton
  piProtonPairs->SetLegPDGs(211,2212);
  piProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(piProtonPairs);

  AliDielectronSignalMC* kaonKaonPairs = new AliDielectronSignalMC("kaonKaonPairs","kaon+kaon pairs");  // mis-id kaon + mis-id kaon
  kaonKaonPairs->SetLegPDGs(321,321);
  kaonKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(kaonKaonPairs);

  AliDielectronSignalMC* kaonProtonPairs = new AliDielectronSignalMC("kaonProtonPairs","kaon+proton pairs");   // mis-id kaon + mis-id proton
  kaonProtonPairs->SetLegPDGs(321,2212);
  kaonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(kaonProtonPairs);

  AliDielectronSignalMC* protonProtonPairs = new AliDielectronSignalMC("protonProtonPairs","proton+proton pairs");  // mis-id proton + mis-id proton
  protonProtonPairs->SetLegPDGs(2212,2212);
  protonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(protonProtonPairs);

  AliDielectronSignalMC* muonAllPairs = new AliDielectronSignalMC("muonAllPairs","muon+everything pairs");        // mis-id muon + something else (electron, pion, kaon, proton)
  muonAllPairs->SetLegPDGs(13,13,kFALSE,kTRUE);
  muonAllPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  die->AddSignalMC(muonAllPairs);



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




