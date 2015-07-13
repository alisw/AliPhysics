void InitHistograms(AliDielectron *die, Int_t cutDefinition);
//void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupCuts(AliDielectron *die, Int_t cutDefinition);

AliESDtrackCuts *SetupTrackCuts(Int_t cutDefinition);
AliAnalysisCuts *SetupPreFilterTrackCuts(Int_t cutDefinition);
AliDielectronPID *SetPIDcuts(Int_t cutDefinition);
AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition);
AliDielectronEventCuts* GetEventCuts();


TString names ("MB_cut23_pf;MB_cut23");

//TString names ("USfiltering_100mrad_mass30MeV_cut23;USfiltering_25mrad_mass60MeV_cut23;USfiltering_50mrad_mass60MeV_cut23");
//TString names ("USfiltering_100mrad_mass60MeV_cut23;USfiltering_50mrad_mass100MeV_cut23;USfiltering_100mrad_mass100MeV_cut23");
//TString names ("multi_0_30_pf_100mrad_mass60MeV_cut23;multi_30_60_pf_100mrad_mass60MeV_cut23;multi_60_1000_pf_100mrad_mass60MeV_cut23");


     Bool_t kRot = 0;
     Bool_t kMix = 1;

     ULong64_t triggerMask = AliTrigger::kINT7;
     Bool_t randomizeDau = kTRUE;
     
     
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Bool_t MCenabled=kFALSE;
const Int_t nPF = 1;

AliDielectron* Config_lowmasspPb(Int_t cutDefinition=0)
{
  //
  // Setup the instance of AliDielectron
  //
//   MCenabled=hasMC;
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
        mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
	mix->SetDepth(30);
	//mix->SetMixUncomplete(kTRUE);
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
//  InitCF(die,cutDefinition);

  // eta correction
  // SetEtaCorrection();

  die->SetNoPairing(kFALSE);
  
  return die;

}

//______________________________________________________________________________________
void SetupCuts(AliDielectron *die, Int_t cutDefinition)
{
  // Setup the track cuts

	//pairing with TLorentzVector
	die->SetUseKF(kFALSE);

//      AliDielectronVarCuts *evTrackletCut = new AliDielectronVarCuts("evTrackletCut","evTrackletCut");
//      if(cutDefinition == 0) evTrackletCut->AddCut(AliDielectronVarManager::kNaccItsPureEsd05Corr,0.,30.-0.2e-5,kFALSE);
//      if(cutDefinition == 1) evTrackletCut->AddCut(AliDielectronVarManager::kNaccItsPureEsd05Corr,30.-.1e-5,60.-0.2e-5,kFALSE);
//      if(cutDefinition == 2) evTrackletCut->AddCut(AliDielectronVarManager::kNaccItsPureEsd05Corr,60.-0.1e-5,10000.,kFALSE);
//      evTrackletCut->AddCut(AliDielectronVarManager::kNaccItsPureEsd05Corr,60.-0.1e-5,10000.,kFALSE);
//      die->GetEventFilter().AddCuts(evTrackletCut);

      AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
      noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
      
      if(cutDefinition < nPF){

	die->GetTrackFilter().AddCuts(SetupPreFilterTrackCuts(cutDefinition));
	//die->GetTrackFilter().AddCuts(SetPreFilterPIDcuts(cutDefinition));
	
	//pairPrefilter
	AliAnalysisCuts* pairPreCuts=0x0;
	
	AliDielectronVarCuts* pairCutsInvM  = new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
        AliDielectronVarCuts* pairCutsOpAng = new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
	
	pairCutsInvM ->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
        pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.050);
	
        AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
        pairCutsCG->AddCut(pairCutsInvM);
        pairCutsCG->AddCut(pairCutsOpAng);
        pairPreCuts = pairCutsCG;
	die->GetPairPreFilter().AddCuts(pairPreCuts);

	//FinalTrackCuts after prefiltering
	die->GetPairPreFilterLegs().AddCuts(SetPIDcuts(cutDefinition));
	die->GetPairPreFilterLegs().AddCuts(SetupTrackCuts(cutDefinition));
	
	die->GetPairPreFilterLegs().AddCuts(noconv);
      
//        AliDielectronVarCuts *PhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
//        PhiV->AddCut(AliDielectronVarManager::kM,        0.0 , 0.1,kTRUE);
//        PhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.0 , 3.2,kTRUE);
//        PhiV->SetCutType(1);
//        die->GetPairFilter().AddCuts(PhiV);
        die->SetPreFilterUnlikeOnly(kTRUE);
	//if(cutDefinition == 0) die->SetPreFilterUnlikeOnly(kTRUE);
	//if(cutDefinition == 1) die->SetPreFilterAllSigns(kTRUE);
      }
      else{
        die->GetTrackFilter().AddCuts(SetupTrackCuts(cutDefinition));
        die->GetTrackFilter().AddCuts(SetPIDcuts(cutDefinition));
        die->GetTrackFilter().AddCuts(noconv);
      }
}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliDielectronPID *SetPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID *pid = new AliDielectronPID();

    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  
 return pid;
}

//______________________________________________________________________________________
AliAnalysisCuts *SetupTrackCuts(Int_t cutDefinition){


  
  // same as ESD cuts without ITS chi2/Ncls
  
  AliDielectronCutGroup* trackCuts=0x0;
  
  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0); // means at least 2 with PID
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,  100.0, 160.0); // crossed rows
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,    120.0, 160.0); // N cls TPC
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.9, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
  
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
  trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  
  
  AliDielectronCutGroup *cgTrackCuts = new AliDielectronCutGroup("cgTrackCuts","cgTrackCuts",AliDielectronCutGroup::kCompAND);
  cgTrackCuts->AddCut(trackCutsDiel);
  cgTrackCuts->AddCut(trackCutsAOD);
  trackCuts = cgTrackCuts;
  
  return trackCuts;
}

AliDielectronPID *SetPreFilterPIDcuts(Int_t cutDefinition){
  
  AliDielectronPID *pid = new AliDielectronPID();

//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
if(cutDefinition == 0)    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
if(cutDefinition == 1)    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
   
 return pid;
}

AliAnalysisCuts *SetupPreFilterTrackCuts(Int_t cutDefinition){
 
  AliDielectronCutGroup* trackCuts=0x0;
  
  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0); // means at least 2 with PID
  
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(1<<2); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
  //trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  
  AliDielectronCutGroup *cgTrackCuts = new AliDielectronCutGroup("cgTrackCuts","cgTrackCuts",AliDielectronCutGroup::kCompAND);
  cgTrackCuts->AddCut(trackCutsDiel);
  cgTrackCuts->AddCut(trackCutsAOD);
  trackCuts = cgTrackCuts;
  

}

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
  //histos->SetReservedWords("Track;Pair");
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack");
  //histos->SetReservedWords("Track");  

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
  if(cutDefinition < nPF)    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
  }
  
  //ME and track rot
  if (die->GetMixingHandler()) {
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
  }
  if (die->GetTrackRotator()) {
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
  }
if(cutDefinition < nPF){  
  //PreFilter Classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Create Classes for Rejected Tracks/Pairs:
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
    // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
    histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
  }
}
  /*
  
  

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
  //if(kRot)histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));//Rot
*/
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",1,0,1,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","ZVertex","ZVertex;ZVertex/cm",480,-12.,12.,AliDielectronVarManager::kZvPrim);
/* 
  histos->UserHistogram("Event","NaccItsPureEsd05Corr","NaccItsPureEsd05Corr;NaccItsPureEsd05Corr/cm",151,-0.5,150.5,AliDielectronVarManager::kNaccItsPureEsd05Corr);
  histos->UserHistogram("Event","NaccItsPureEsd05","NaccItsPureEsd05;NaccItsPureEsd05/cm",151,-0.5,150.5,AliDielectronVarManager::kNaccItsPureEsd05);
  histos->UserHistogram("Event","kNaccTrckltsEsd05Corr","kNaccTrckltsEsd05Corr;kNaccTrckltsEsd05Corr/cm",151,-0.5,150.5,AliDielectronVarManager::kNaccTrckltsEsd05Corr);
  histos->UserHistogram("Event","kNaccTrckltsEsd05","kNaccTrckltsEsd05;kNaccTrckltsEsd05/cm",151,-0.5,150.5,AliDielectronVarManager::kNaccTrckltsEsd05);
 
  histos->UserHistogram("Event","NTrk_Nacc05","NTrk vs. Nacc05; NTrk; Nacc05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kNTrk,AliDielectronVarManager::kNacc05);
  
  histos->UserHistogram("Event","MultV0A_Nacc","MultV0A vs. Nacc; MultV0A; Nacc",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0A,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","MultV0C_Nacc","MultV0C vs. Nacc; MultV0C; Nacc",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0C,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","MultV0_Nacc", "MultV0  vs. Nacc; MultV0;  Nacc",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0, AliDielectronVarManager::kNacc);
 
  histos->UserHistogram("Event","MultV0A_Nacc05","MultV0A vs. Nacc05; MultV0A; Nacc05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0A,AliDielectronVarManager::kNacc05);
  histos->UserHistogram("Event","MultV0C_Nacc05","MultV0C vs. Nacc05; MultV0C; Nacc05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0C,AliDielectronVarManager::kNacc05);
  histos->UserHistogram("Event","MultV0_Nacc05", "MultV0  vs. Nacc05; MultV0;  Nacc05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kMultV0, AliDielectronVarManager::kNacc05);
    
  histos->UserHistogram("Event","Nacc05_NaccTrckltsEsd05","Nacc05 vs. NaccTrckltsEsd05; Nacc05; NaccTrckltsEsd05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kNacc05,AliDielectronVarManager::kNaccTrckltsEsd05);
  histos->UserHistogram("Event","Nacc05_NaccItsTpcEsd05", "Nacc05 vs. NaccItsTpcEsd05; Nacc05; NaccItsTpcEsd05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kNacc05,AliDielectronVarManager::kNaccItsTpcEsd05);
  histos->UserHistogram("Event","Nacc05_NaccItsPureEsd05","Nacc05 vs. NaccItsPureEsd05; Nacc05; NaccItsPureEsd05",1000,-0.5,999.5,1000,-0.5,999.5,AliDielectronVarManager::kNacc05,AliDielectronVarManager::kNaccItsPureEsd05);
  histos->UserHistogram("Event","NaccItsPureEsd05","NaccItsPureEsd05; NaccItsPureEsd05",1000,-0.5,999.5,AliDielectronVarManager::kNaccItsPureEsd05);
*/
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",1000,0.,10.,AliDielectronVarManager::kPt);
//  histos->UserHistogram("Track","P","P;P [GeV];#tracks",500,0.,10.,AliDielectronVarManager::kP);
  histos->UserHistogram("Track","Eta_phi","Eta vs Phi;Eta;Phi",90,-0.9,0.9,160,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Pt_phi","Pt vs Phi;Pt;Phi [GeV];#tracks",500,0.,5.,320,0.,6.4,AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
 // histos->UserHistogram("Track","ImpParXY","ImpParXY; ImpParXY ;#tracks",500,-5.,5.,AliDielectronVarManager::kImpactParXY);
 // histos->UserHistogram("Track","ImpParZ","ImpParZ; ImpParZ ;#tracks",500,-5.,5.,AliDielectronVarManager::kImpactParZ);

  histos->UserHistogram("Track","NclsITS","NclsITS ;kNclsITS ;#tracks",6,2,8,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","CrossedRowsOverFindable","CrRowsOverFindable; CrRows/FindableCls ;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCchi2perCls","TPCchi2perCls; TPCchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","CrossedRows","CrossedRows; CrossedRows ;#tracks",200,0,200,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","NClusterTPC","NClusterTPC; NClusterTPC ;#tracks",200,0,200,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","NClusterITS","NClusterITS; NClusterITS ;#tracks",10,0,10,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2perCls","ITSchi2perCls; ITSchi2perCls ;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
 // histos->UserHistogram("Track","TPCsignalNfrac",";fraction TPCSignalN/TPCncls ;TPCSignalN/TPCncls;#tracks",60,0.,1.2,AliDielectronVarManager::kTPCsignalNfrac);
 // histos->UserHistogram("Track","TPCsignalN","TPCsignalN;TPCsignalN;#tracks",200,0.,200,AliDielectronVarManager::kTPCsignalN);
/*
  histos->UserHistogram("Track","TPCnSigma_MomEle","TPC number of sigmas Electrons vs Momentum;Mom;TPCsigmaEle",1000,0.,10.,800,-40,40,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","ITSdEdx_P","dEdx;P [GeV];ITS signal (arb units) vs Momentum;Mom;ITSsignal", 400,0.,8.,450,0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","TPCdEdx_P","dEdx;P [GeV];TPC signal (arb units) vs Momentum;Mom;TPCsignal", 400,0.,8.,450,0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TOFbeta_Mom","kTOFbeta vs Momentum;Mom;TOFbeta"                           , 400,0.,8.,120,0.,  1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
*/
//  histos->UserHistogram("Track","IsPureITS",";IsPureITS ;#tracks",3,-1.5,1.5,AliDielectronVarManager::kITSpureSA);
  //
  //add histograms to Pair classes
  //

//  histos->UserHistogram("Pair","InvMass_pPt","Inv.Mass vs. pair p_{T};Inv. Mass (GeV/c^{2});pair p_{T} (GeV/c)",
//                        800,0.,8.,800,0.,8.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  //histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
  //                      650,0.,6.5, 200,0.,10., 100,0.,TMath::Pi(),
  //                      AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);

  
  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        800,0.,8., 400,0.,8., 32,0.,TMath::Pi(),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
			
  //histos->UserHistogram("Pair","InvMass_pPt_FineBinning","Inv.Mass vs. pair p_{T};Inv. Mass (GeV/c^{2});pair p_{T} (GeV/c)",
  //                      4000,0.,0.4,50,0.,5.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);	
/*			
  histos->UserHistogram("Pair","Eta_phi_pair","Eta vs Phi (pair);Eta;Phi",100,-1.,1.,320,0.,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Pair",
                        "kDeltaEta_kDeltaPhi","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        160, 0. , 1.6, 160 , 0., 3.2 ,
                         AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );
*/
//   histos->UserHistogram("Pair",
//                         "OpAngle","opAngle;counts;opAngle",
//                         3000, -1.5. , 1.5, AliDielectronVarManager::kOpeningAngle);
			
  histos->UserHistogram("Pair",
                        "InvMass_PhivPair","InvMass_PhivPair;InvMass;PhivPair",
                         50, 0. , 0.5, 160 , 0., 3.2 ,
                         AliDielectronVarManager::kM , AliDielectronVarManager::kPhivPair );
/*
  histos->UserHistogram("Pair",
                        "OpAngle_PhivPair","OpAngle_PhivPair;Opening Angle;PhivPair",
                        320, 0. , 3.2, 320 , 0., 3.2,
                         AliDielectronVarManager::kOpeningAngle , AliDielectronVarManager::kPhivPair );
*/
  histos->UserHistogram("Pair",
			"OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",
			500, 0. , 2.5, 350 , 0. , TMath::Pi() ,
			AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  /*
  histos->UserHistogram("Pair",
                        "DeltaEta_DeltaPhi_fine","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        600, 0. , 0.3, 640 , 0., .32 ,
                        AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );
  histos->UserHistogram("Pair",
                        "OpAngle_InvMass_fine","InvMass_openingAngle;Invariant Mass;opening angle",
                        1000, 0. , 0.05, 1000 , 0. , 0.05 ,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);*/
//  histos->UserHistogram("Pair",
//                        "Y","Y;counts;Y",
//                        120, -1.2 , 1.2, AliDielectronVarManager::kY);

   // 2D-PID
  //histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle","2D PID - TPC vs ITS;P [GeV];TPC number of sigmas Electrons;ITS number of sigmas Electrons",100,0.,5.,320,-12.,20.,300,-10.,20.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle,kTRUE);
  //histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle","2D PID - TPC vs TOF;P [GeV];TPC number of sigmas Electrons;TOF number of sigmas Electrons",100,0.,5.,320,-12.,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle,kTRUE);

if(cutDefinition < nPF){ 
  histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
//  histos->UserHistogram("Pre","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
//  histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
//  histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
    
//  histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
//                        400,0.,8.,450,0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal);
//  histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
//                        400,0.,8.,450,0.,150. ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
//  histos->UserHistogram("Pre","IsPureITS",";IsPureITS ;#tracks",3,-1.5,1.5,AliDielectronVarManager::kITSpureSA);
			
  histos->UserHistogram("RejPair","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  /*
  histos->UserHistogram("RejPair",
                        "DeltaEta_DeltaPhi","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        160, 0. , 1.6, 320 , 0., 3.2 ,
                         AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );
  histos->UserHistogram("RejPair",
			"OpAngle_InvMass","OpAngle_InvMass;Opening angle;Invariant Mass",
			320, 0. , 3.2, 500 , 0. , 5. ,
			AliDielectronVarManager::kOpeningAngle,AliDielectronVarManager::kM);*/
  histos->UserHistogram("RejPair",
                        "DeltaEta_DeltaPhi_fine","kDeltaEta_kDeltaPhi;kDeltaEta;kDeltaPhi",
                        300, 0. , 0.3, 640 , 0., .32 ,
                        AliDielectronVarManager::kDeltaEta , AliDielectronVarManager::kDeltaPhi );
  histos->UserHistogram("RejPair",
                        "OpAngle_InvMass_fine","InvMass_openingAngle;Invariant Mass;opening angle",
                        500, 0. , 0.05, 500 , 0. , 0.05 ,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("RejTrack","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
//  histos->UserHistogram("RejTrack","IsPureITS",";IsPureITS ;#tracks",3,-1.5,1.5,AliDielectronVarManager::kITSpureSA);
  histos->UserHistogram("Track_Legs","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt); 
}
  die->SetHistogramManager(histos);

}

AliDielectronEventCuts* GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1); 
  
  return eventCuts;
}


