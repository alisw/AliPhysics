//#include "PWGDQ/dielectron/macrosLMEE/LMEECutLib.C"

void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void EnableMC();

TString names=("noPairing;TPCTOFCentwidenoRej;TPCTOFSemiCent1noRej;TPCTOFSemiCent2noRej;TPCTOFPerinoRej;TPCTOFCentRP;TPCTOFSemiCent1RP;TPCTOFSemiCent2RP;TPCTOFPeriRP;TPCTOFCentMag;TPCTOFSemiCent1Mag;TPCTOFSemiCent2Mag;TPCTOFPeriMag;TPCTOFCentnoTOF;NoPIDNoPairing;TPCTOFCentnoRej");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t MCenabled=kFALSE;


AliDielectron* ConfigAsakoLMEEPbPb2011AOD(Int_t cutDefinition, Bool_t hasMC=kFALSE, Bool_t ESDanalysis=kFALSE)
{

  Int_t selectedPID=-1;
  Int_t selectedCentrality=-1;
  Int_t selectedPairCut=-1;
  //  Int_t selectedPairMCut=-1;
  Bool_t rejectionStep=kFALSE;
  Bool_t PairCut=kFALSE; //add by asako
  LMEECutLibAsako*  LMCL = new LMEECutLibAsako();

  //
  // Setup the instance of AliDielectron
  //

  MCenabled=hasMC;
  // create the actual framework object

  TString name=Form("%02d",cutDefinition);
  if ((cutDefinition)<arrNames->GetEntriesFast()){
	name=arrNames->At((cutDefinition))->GetName();
  }

  //thisCut only relevant for MC:
  AliDielectron *die =
	new AliDielectron(Form
					  ("%s",name.Data()),
					  Form("Track cuts: %s",name.Data()));


  //Setup AnalysisSelection:
  if (cutDefinition==0) {
	//not yet implemented
  }
  else if (cutDefinition==1) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
	rejectionStep = kFALSE;
	PairCut = kFALSE;

  }
  else if (cutDefinition==2) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral1;
	rejectionStep = kFALSE;
	PairCut=kFALSE;
  }
  else if (cutDefinition==3) {
    selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral2;
    rejectionStep = kFALSE;
	PairCut=kFALSE;
  }

  else if (cutDefinition==4) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Peripheral;
	rejectionStep = kFALSE;
	PairCut=kFALSE;
  }

  ///////////////////////////////////////////////////////////
  else if (cutDefinition==5) {
    selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
	selectedPairCut = LMEECutLibAsako::kPbPb2011RP;
	//selectedPairMCut = LMEECutLibAsako::kPbPb2011MassAll;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  else if (cutDefinition==6) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral1;
	selectedPairCut = LMEECutLibAsako::kPbPb2011RP;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  else if (cutDefinition==7) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral2;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	selectedPairCut = LMEECutLibAsako::kPbPb2011RP;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  else if (cutDefinition==8) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Peripheral;
	selectedPairCut = LMEECutLibAsako::kPbPb2011RP;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }

  /////////////////////////////////////////////////////////
  else if (cutDefinition==9) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
	selectedPairCut = LMEECutLibAsako::kPbPb2011Mag;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassAll;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  
  else if (cutDefinition==10) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral1;
	selectedPairCut = LMEECutLibAsako::kPbPb2011Mag;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  else if (cutDefinition==11) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011SemiCentral2;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	selectedPairCut = LMEECutLibAsako::kPbPb2011Mag;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }
  else if (cutDefinition==12) {
	selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOFwide;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Peripheral;
	selectedPairCut = LMEECutLibAsako::kPbPb2011Mag;
	// selectedPairMCut = LMEECutLibAsako::kPbPb2011MassMiddle;
	rejectionStep = kFALSE;
	PairCut=kTRUE;
  }

  

  ///////////////////////////////////////////////////////
  else if (cutDefinition==13) {
	//selectedPID = LMEECutLib::kPbPb2011TPCandTOFwide;
	selectedPID = LMEECutLibAsako::kPbPb2011TPC;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
	rejectionStep = kFALSE;
	PairCut=kFALSE;
  }
  else if (cutDefinition==14) {
	selectedPID = LMEECutLibAsako::kPbPb2011NoPID;
	selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
	rejectionStep = kFALSE;
	PairCut=kFALSE;
  }

  else if (cutDefinition==15) {
    selectedPID = LMEECutLibAsako::kPbPb2011TPCandTOF;
    selectedCentrality = LMEECutLibAsako::kPbPb2011Central;
    rejectionStep = kFALSE;
    PairCut = kFALSE;
  }


  else Semi{
	cout << " =============================== " << endl;
	cout << " ==== INVALID CONFIGURATION ==== " << endl;
	cout << " =============================== " << endl;
  }


  //Now configure task

  //Apply correct Pre-Filter Scheme, if necessary
  die->SetPreFilterAllSigns();

  //switch off KF PArticle:
  die->SetUseKF(kFALSE);

  if (selectedPID == LMEECutLibAsako::kPbPb2011NoPID) {
	die->SetNoPairing();
  }

  //  if (ESDanalysis) {
  //die->GetTrackFilter().AddCuts( LMCL->GetESDTrackCutsAna(selectedPID) );
  // }

  die->GetEventFilter().AddCuts(LMCL->GetCentralityCuts(selectedCentrality));
  die->GetTrackFilter().AddCuts( LMCL->GetTrackCutsAna(selectedPID) );
  die->GetTrackFilter().AddCuts( LMCL->GetPIDCutsAna(selectedPID) );

  if(PairCut){
	if (rejectionStep) {
	  die->GetPairPreFilterLegs().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
	  die->GetPairPreFilter().AddCuts( LMCL->GetPairPreFilterCuts(selectedPairCut));
	  die->GetPairFilter().AddCuts( LMCL->GetPairCuts(selectedPairCut));
	}
	else {   
	  die->GetPairFilter().AddCuts( LMCL->GetPairCutsInOut(selectedPairCut));
	  //  die->GetPairFilter().AddCuts( LMCL->GetPairCuts4(selectedPairMCut));
	}
  }
  



  AliDielectronTrackRotator *rot= 0x0;
  /*AliDielectronTrackRotator *rot= LMCL->GetTrackRotator(selectedPID);
  die->SetTrackRotator(rot);
  */
  AliDielectronMixingHandler *mix=LMCL->GetMixingHandler(selectedPID);
  die->SetMixingHandler(mix);

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);

  // the last definition uses no cuts and only the QA histograms should be filled!
  //  InitCF(die,cutDefinition);

  return die;
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
  histos->SetReservedWords("Track;Pair;Pre;RejTrack;RejPair");

  //Event class
  //  if (cutDefinition==nDie-1) 
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

  //PreFilter Classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
	histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
  }


  //Create Classes for Rejected Tracks/Pairs:
  for (Int_t i=0; i<2; ++i){
	histos->AddClass(Form("RejTrack_%s",AliDielectron::TrackClassName(i)));
  }
  for (Int_t i=0; i<3; ++i){
	histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
  }

  /*
  //track rotation

  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  */
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
						1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]","0,10,20,40,80,100,101",
						AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]",100, 0, 100,
						AliDielectronVarManager::kCentrality);

  histos->UserHistogram("Event","v0ACrpH2_old","VZERO-AC;v0ACrpH2",
						100,-2.0,2.0,
						AliDielectronVarManager::kv0ACrpH2);

  histos->UserHistogram("Event","v0ArpH2_old","VZERO-A;v0ArpH2",
						100,-2.0,2.0,
						AliDielectronVarManager::kv0ArpH2);
  histos->UserHistogram("Event","v0CrpH2_old","VZERO-C;v0CrpH2",
						100,-2.0,2.0,
						AliDielectronVarManager::kv0CrpH2);

  histos->UserHistogram("Event","Corr_v0ACrpH2_old","CORR VZERO-AC RP;#Psi_{2}^{V0A} (rad.);#Psi_{2}^{V0C} (rad.)",
						100,-2.0,2.0.,100,-2.0,2.0,
						AliDielectronVarManager::kv0ArpH2,AliDielectronVarManager::kv0CrpH2);

  histos->UserHistogram("Event","v0ACrpH2","VZERO-AC;v0ACrpH2",
                        100,0.,3.15,
                        AliDielectronVarManager::kv0ACrpH2);

  histos->UserHistogram("Event","v0ArpH2","VZERO-A;v0ArpH2",
                        100,0.,3.15,
                        AliDielectronVarManager::kv0ArpH2);
  histos->UserHistogram("Event","v0CrpH2","VZERO-C;v0CrpH2",
                        100,0.,3.15,
                        AliDielectronVarManager::kv0CrpH2);

  histos->UserHistogram("Event","Corr_v0ACrpH2","CORR VZERO-AC RP;#Psi_{2}^{V0A} (rad.);#Psi_{2}^{V0C} (rad.)",
                        100,0.,3.15.,100,0.,3.15,
                        AliDielectronVarManager::kv0ArpH2,AliDielectronVarManager::kv0CrpH2);
  
	//add histograms to Track classes
  /*  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px","Px;Px [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py","Py;Py [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz","Pz;Pz [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPz);

  histos->UserHistogram("Track","NclsSFracTPC","NclsSFracTPC; NclsSFracTPC;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff","TPCclsDiff; TPCclsDiff;#tracks",200,0,10.,AliDielectronVarManager::kTPCclsDiff);
  */
  histos->UserHistogram("Track","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
    400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
    400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
  /*
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
    400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
    400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
    400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);

  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons;#tracks",
    400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,kTRUE);
  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions;#tracks",
    400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,kTRUE);

  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
    400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
    400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);

  histos->UserHistogram("Track","TOFbeta","TOF beta;P [GeV];TOF beta;#tracks",
    400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);

  
  histos->UserHistogram("Track","Eta","Eta; Eta;#tracks",
    200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi; Phi;#tracks",
    200,0.,3.15,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
    200,-2,2,200,0,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dXY_dZ","dXY dZ Map; dXY; dZ;#tracks",
    200,-2,2,200,-2,2.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);


  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParZ);

    histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",100,0.,1.,AliDielectronVarManager::kNFclsTPCfCross);
	  histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",159,0.,159.,AliDielectronVarManager::kNFclsTPCr);

  	
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);
		  histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;ITS number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsITS);

		    histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC chi2;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
			  histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS chi2;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);

			    histos->UserHistogram("Track","TPCnCls_kNFclsTPCr","nTPC vs nTPCr;nTPC vs nTPCr;#tracks",159,0.,159.,159,0.,159.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);

				  histos->UserHistogram("Track","kNFclsTPCr_pT","nTPCr vs pt;nTPCr vs pt;#tracks",159,0.,159.,200,0.,20.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kPt);
				    
  */
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMassAll","Inv.Mass;Inv. Mass [GeV];#pairs",
						500,0.0,5.00,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","InvMassLow","Inv.Mass;Inv. Mass [GeV];#pairs",
						300,0.0,0.03,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","InvMassMiddle","Inv.Mass;Inv. Mass [GeV];#pairs",
						180,0.12,0.3,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","InvMassHigh","Inv.Mass;Inv. Mass [GeV];#pairs",
						200,0.3,0.5,AliDielectronVarManager::kM);

    
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
						100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
						100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","PhivPair","PhiV ;#",
						100,0.,3.15,AliDielectronVarManager::kPhivPair );
    
  histos->UserHistogram("Pair","Pt","Pt;Pt [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Px","Px;Px [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Pair","Py","Py;Py [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Pair","Pz","Pz;Pz [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Pair","Phi","Phi;Phi[rad];#counts",100,-3.15,3.15,AliDielectronVarManager::kPhi );
    
    
  histos->UserHistogram("Pair","DeltaPhiv0ArpH2","Phi;Phi[rad];#counts",
						100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0ArpH2);
  histos->UserHistogram("Pair","DeltaPhiv0CrpH2","Phi;Phi[rad];#counts",
						100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0CrpH2);
  histos->UserHistogram("Pair","DeltaPhiv0ACrpH2","Phi;Phi[rad];#counts",
						100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0ACrpH2);

    
  histos->UserHistogram("Pair","PairPlanev0ACrpH2Angle","Phi;Phi[rad];#counts",
						100,0,1.6,AliDielectronVarManager::kPairPlanev0rpH2Angle);
  histos->UserHistogram("Pair","PairPlaneMagAngle","Phi;Phi[rad];#counts",
						100,0,1.6,AliDielectronVarManager::kPairPlaneMagAngle);
  histos->UserHistogram("Pair","PairPlaneAngle","Phi;Phi[rad];#counts",
						100,0,1.6,AliDielectronVarManager::kPairPlaneAngle);


  histos->UserHistogram("Pair","RotationPair_x","Px_Rot;Px_Rot[GeV];#tracks",
						300,0,30,AliDielectronVarManager::kRotPairx);
  histos->UserHistogram("Pair","RotationPair_y","Px_Rot;Py_Rot[GeV];#tracks",
						300,0,30,AliDielectronVarManager::kRotPairy);
  histos->UserHistogram("Pair","RotationPair_z","Px_Rot;Pz_Rot[GeV];#tracks",
						300,0,30,AliDielectronVarManager::kRotPairz);

  
    
  //2D Histo Plot
  histos->UserHistogram("Pair","InvMassALLPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
						1000,0.0,5.0,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMassALL2PairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
						500,0.0,0.5,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMassLowPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
						300,0.0,0.03,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMassMiddlePairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
						180,0.12,0.3,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMassHighPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
						200,0.3,0.5,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
    
  
  /*
  histos->UserHistogram("Pair","InvMassALLPhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
						1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMassALL2PhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
						500,0.0,0.50,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMassLowPhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
						300,0.0,0.03,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMassMiddlePhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
						180,0.12,0.3,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMassHighPhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
						200,0.3,0.50,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);

  histos->UserHistogram("Pair","InvMassAllOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
						1000,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMassAll2OpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
						500,0.0,0.5,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMassLowOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
						300,0.0,0.03,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMassMiddleOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
						180,0.12,0.3,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMassHighOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
						200,0.3,0.5,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  */


  histos->UserHistogram("Pair","InvMassAllPairplaneMagAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						1000,0.0,5.0,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagAngle);
  histos->UserHistogram("Pair","InvMassAll2PairplaneMagAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagAngle);
  histos->UserHistogram("Pair","InvMassLoweeplaneMagAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						300,0.0,0.03,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagAngle);
  histos->UserHistogram("Pair","InvMassMiddleeelaneMagAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						180,0.12,0.3,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagAngle);
  histos->UserHistogram("Pair","InvMassHighPairplaneMagAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						200,0.3,0.5,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagAngle);

  histos->UserHistogram("Pair","InvMassAllPairplaneRPAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						1000,0.0,5.0,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlanev0rpH2Angle);
  histos->UserHistogram("Pair","InvMassAll2PairplaneRPAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						500,0.0,0.5,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlanev0rpH2Angle);
  histos->UserHistogram("Pair","InvMassLoweeplaneRPAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						300,0.0,0.03,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlanev0rpH2Angle);
  histos->UserHistogram("Pair","InvMassMiddleeelaneRPAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						180,0.12,0.3,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlanev0rpH2Angle);
  histos->UserHistogram("Pair","InvMassHighPairplaneRPAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
						200,0.3,0.5,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlanev0rpH2Angle);
    
  histos->UserHistogram("Pair","InvMassAllPairplaneAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        1000,0.0,5.0,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle);
  histos->UserHistogram("Pair","InvMassAll2PairplaneAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle);
  histos->UserHistogram("Pair","InvMassLoweeplaneAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        300,0.0,0.03,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle);
  histos->UserHistogram("Pair","InvMassMiddleeelaneAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        180,0.12,0.3,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle);
  histos->UserHistogram("Pair","InvMassHighPairplaneAngle","ee plane and Mag Vector Angle vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                        200,0.3,0.5,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle);

     

  //add histograms to Track classes
  /*  histos->UserHistogram("Pre","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);

  histos->UserHistogram("Pre","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
						400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Pre","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
						400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);


  histos->UserHistogram("Pre","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
						400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Pre","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
						400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Pre","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
						400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);

  histos->UserHistogram("Pre","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons;#tracks",
						400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,kTRUE);
  histos->UserHistogram("Pre","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions;#tracks",
						400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,kTRUE);

  histos->UserHistogram("Pre","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
						400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Pre","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
						400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);

  histos->UserHistogram("Pre","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
						200,-2,2,200,0,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Pre","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Pre","ZVertex ","ZVertex ;ZVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kZv);
  histos->UserHistogram("Pre","XVertex ","XVertex ;XVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kXv);
  histos->UserHistogram("Pre","YVertex ","YVertex ;YVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kYv);

  histos->UserHistogram("Pre","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);
  */    
  //add histograms to Pair classes For Rejected Pairs:
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,200,0,20);
  // cf->AddVariable(AliDielectronVarManager::kPt,200,0,20);// added
  cf->AddVariable(AliDielectronVarManager::kM,201,-0.01,4.01); //20Mev Steps
  //cf->AddVariable(AliDielectronVarManager::kM, 1000,0,5.00); //5Mev Steps 
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //cf->AddVariable(AliDielectronVarManager::kCentrality,100, 0, 100);


   cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,10.0,30.0,40.0,60.,80.,100.");

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,200,0.,20.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,1000,0.0.,1000.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,500,0.0.,500.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);

  //only in this case write MC truth info
  if (MCenabled) {
	cf->SetStepForMCtruth();
	cf->SetStepsForMCtruthOnly();
	cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
	cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}

//--------------------------------------
void EnableMC() {
  MCenabled=kTRUE;
}



