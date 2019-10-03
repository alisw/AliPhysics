void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void EnableMC();
void SetSignals(AliDielectron *die);

//////////////////////////////////////////////////////
///// lists of cuts 
////  |cutDefinition | CENT | PID | Pair |
///// | 0 | - | - | - |
///// | 1 | 0-10 | kPbPb2011TPCandTOFwide | - | 
///// | 2 | 10-30 | kPbPb2011TPCandTOFwide | - | 
///// | 3 | 30-50 | kPbPb2011TPCandTOFwide | - | 
///// | 4 | 0-10 | kPbPb2011TPCandTOFwide | OP & SetPreFilterAllSigns | 
///// | 5 | 10-30 | kPbPb2011TPCandTOFwide | OP & SetPreFilterAllSigns | 
///// | 6 | 30-50 | kPbPb2011TPCandTOFwide | OP & SetPreFilterAllSigns | 
///// | 7 | 0-10 | kPbPb2011TPCandTOFwide | Phiv & SetPreFilterAllSigns | 
///// | 8 | 10-30 | kPbPb2011TPCandTOFwide | Phiv & SetPreFilterAllSigns | 
///// | 9 | 30-50 | kPbPb2011TPCandTOFwide | Phiv & SetPreFilterAllSigns | 
///// | 10 | 0-10 | kPbPb2011TPCandTOFwide | OP & No SetPreFilterAllSigns | 
///// | 11 | 10-30 | kPbPb2011TPCandTOFwide | OP & No SetPreFilterAllSigns | 
///// | 12 | 30-50 | kPbPb2011TPCandTOFwide | OP & No SetPreFilterAllSigns | 
///// | 13 | 0-10 | kPbPb2011TPCandTOFwide | Phiv & No SetPreFilterAllSigns | 
///// | 14 | 10-30 | kPbPb2011TPCandTOFwide | Phiv & No SetPreFilterAllSigns | 
///// | 15 | 30-50 | kPbPb2011TPCandTOFwide | Phiv & No SetPreFilterAllSigns | 
//////////////////////////////////////////////////////////


TString names=("noPairing;TPCTOFCentnoPair;TPCTOFSemiCent1noPair;TPCTOFSemiCent2noPair;TPCTOFCentOpPair;TPCTOFSemiCent1OpPair;TPCTOFSemiCent2OpPair;TPCTOFCentPvPair;TPCTOFSemiCent1PvPair;TPCTOFSemiCent2PvPair;TPCTOFCentOpPair2;TPCTOFSemiCent1OpPair2;TPCTOFSemiCent2OpPair2;TPCTOFCentPvPair2;TPCTOFSemiCent1PvPair2;TPCTOFSemiCent2PvPair2");
                
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t MCenabled=kFALSE;


AliDielectron* ConfigTakuLMEEPbPb2011(Int_t cutDefinition, Bool_t withMC=kFALSE, Bool_t CFenable=kFALSE)
{

  Int_t selectedPID=-1;
  Int_t selectedCentrality=-1;
  Int_t selectedPairCut=-1;
  Bool_t rejectionStep=kFALSE;
  Bool_t HasPairCut = kFALSE;
  LMEECutLibTaku*  LMCL = new LMEECutLibTaku();

  //
  // Setup the instance of AliDielectron
  //

  MCenabled=withMC;
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
  
  
  if (MCenabled){
    die->SetHasMC(kTRUE);
  }
  
  //Setup AnalysisSelection:
  if (cutDefinition==0) {
    //not yet implemented
  }
  //////////////////////////////////////////////
  else if (cutDefinition==1) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011Central;
    rejectionStep = kFALSE;
    HasPairCut = kFALSE;
  }
  else if (cutDefinition==2) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral1;
    rejectionStep = kFALSE;
    HasPairCut = kFALSE;
  }
  else if (cutDefinition==3) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral2;
    rejectionStep = kFALSE;
    HasPairCut = kFALSE;
  }
  //////////////////////////////////////////////
  else if (cutDefinition==4) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011Central;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==5) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral1;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==6) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral2;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  //////////////////////////////////////////////
  else if (cutDefinition==7) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011Central;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==8) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral1;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==9) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral2;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kTRUE;
    HasPairCut =  kTRUE;
  }
  //////////////////////////////////////////////
  else if (cutDefinition==10) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011Central;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kFALSE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==11) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral1;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kFALSE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==12) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral2;
    selectedPairCut = LMEECutLibTaku::kPbPb2011OP;
    rejectionStep = kFALSE;
    HasPairCut =  kTRUE;
  }
  //////////////////////////////////////////////
  else if (cutDefinition==13) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011Central;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kFALSE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==14) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral1;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kFALSE;
    HasPairCut =  kTRUE;
  }
  else if (cutDefinition==15) {
    selectedPID = LMEECutLibTaku::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibTaku::kPbPb2011SemiCentral2;
    selectedPairCut = LMEECutLibTaku::kPbPb2011Phiv;
    rejectionStep = kFALSE;    
    HasPairCut =  kTRUE;
  }
  //////////////////////////////////////////////



  //Legacy cuts, check consistence w/ 1 & 4, then remove
  /*
  else if (cutDefinition==7) {
    selectedPID = LMEECutLib::kPbPb2011TPCandTOF;
    selectedCentrality = LMEECutLib::kPbPb2011Central;
    rejectionStep = kFALSE;
  }
  else if (cutDefinition==8) {
    selectedPID = LMEECutLib::kPbPb2011TPCandTOF;
    selectedCentrality = LMEECutLib::kPbPb2011Central;
    rejectionStep = kTRUE;
  }
  */
  else Semi{
    cout << " =============================== " << endl;
    cout << " ==== INVALID CONFIGURATION ==== " << endl;
    cout << " =============================== " << endl;
  }


  //Now configure task

  /*
  if (rejectionStep) {
    die->GetTrackFilter().AddCuts(LMCL->GetPIDCutsPre(selectedPID) );
    die->GetPairPreFilterLegs().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
    die->GetPairPreFilter().AddCuts(LMCL->GetPairCuts(selectedPID) );
  }
  else { //No Prefilter, no Pairfilter
  }
  */


  //Apply correct Pre-Filter Scheme, if necessary
  die->SetPreFilterAllSigns(rejectionStep);

  die->GetTrackFilter().AddCuts( LMCL->GetPIDCutsAna(selectedPID) );
  die->GetEventFilter().AddCuts(LMCL->GetCentralityCuts(selectedCentrality));

  if(HasPairCut==kTRUE){
    if(rejectionStep==kTRUE){
      die->GetPairPreFilterLegs().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
      die->GetPairPreFilter().AddCuts( LMCL->GetPairPreFilterCuts(selectedPairCut));
      die->GetPairFilter().AddCuts( LMCL->GetPairCuts(selectedPairCut));
    }else{
      die->GetPairFilter().AddCuts( LMCL->GetPairCuts(selectedPairCut));
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
  if (CFenable) {
    SetSignals(die);
    InitCF(die,cutDefinition);
  }
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
  histos->SetReservedWords("Track;Pair;Pre;RejTrack;RejPair:MyPair:MyPairV0");

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
  
  
  //My Pair classes
  // to fill also mixed event histograms loop until 10

  const char* fgkPairClassNames[7] = {
    "unlike",
    "like_pp",
    "like_ee",
    "mixunlike_pe",
    "mixunlike_ep",
    "mixlike_pp",
    "mixlike_ee"
  };
  for (Int_t i=0; i<7; ++i){
    histos->AddClass(Form("MyPair_%s",fgkPairClassNames[i]));
  }
  for (Int_t i=0; i<7; ++i){
    histos->AddClass(Form("MyPairV0_%s",fgkPairClassNames[i]));
  }


  /*
  //track rotation

  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  */
  //add histograms to event class
  histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
			1,0.,1.,AliDielectronVarManager::kNevents);
  //histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]","0,10,20,40,80,100,101",
  //AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]",100, 0, 100,
			AliDielectronVarManager::kCentrality);
  

  //add histograms to Track classes, also fills RejTrack
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","NclsSFracTPC","NclsSFracTPC; NclsSFracTPC;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff","TPCclsDiff; TPCclsDiff;#tracks",200,0,10.,AliDielectronVarManager::kTPCclsDiff);

  histos->UserHistogram("Track","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
	  400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
	  400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

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

  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
			200,-1,1,200,0,6.3,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);

  histos->UserHistogram("Track","TPCnCls_kNFclsTPCr","nTPC vs nTPCr;nTPC vs nTPCr;#tracks",159,0.,159.,159,0.,159.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);

  histos->UserHistogram("Track","kNFclsTPCr_pT","nTPCr vs pt;nTPCr vs pt;#tracks",159,0.,159.,200,0.,20.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kPt);

  //add histograms to Pair classes, also fills RejPair
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
			1000,0.0,5.00,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
			100,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
			100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  //2D Histo Plot
  histos->UserHistogram("Pair","InvMassPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
			1000,0.0,5.0,200,0.,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMassPhivPair","PhivPair vs Inv. Mass;Inv. Mass [GeV]; Phiv",
			1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMassOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV] ; OpeningAngle[rad]",
			1000,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  //add histograms to Track classes
  histos->UserHistogram("Pre","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);

  histos->UserHistogram("Pre","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
	  400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Pre","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
	  400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Pre","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
			200,-1,1,200,0,6.3,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Pre","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Pre","ZVertex ","ZVertex ;ZVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kZv);
  histos->UserHistogram("Pre","XVertex ","XVertex ;XVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kXv);
  histos->UserHistogram("Pre","YVertex ","YVertex ;YVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kYv);

  histos->UserHistogram("Pre","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
			400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Pre","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
			400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Pre","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
			400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);
  histos->UserHistogram("Pre","TPCnSigmaPro_P","TPC number of sigmas Proton;P [GeV];TPC number of sigmas Pions;#tracks",
			400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,kTRUE);


//  histos->UserHistogram("Pre","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);


  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////// my own histograms 
  histos->UserHistogram("MyPair","InvMyMassPairPt","Inv.MyMass vs PairPt;Inv. MyMass [GeV]; pT [GeV]",
			1000,0.0,5.0,200,0.,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

  histos->UserHistogram("MyPair","InvMyMassOpeningAngle","Opening Angle vs Inv.MyMass;Inv. MyMass [GeV] ; OpeningAngle[rad]",
			1000,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("MyPair","InvMyMassPsiPair","PsiPair vs Inv.MyMass;Inv. MyMass [GeV]; Phiv",
			1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPsiPair);

  histos->UserHistogram("MyPair","InvMyMassPhivPair","PhivPair vs Inv.MyMass;Inv. MyMass [GeV]; Phiv",
			1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);

  histos->UserHistogram("MyPair","InvMyMassR","R vs Inv.MyMass;Inv. MyMass [GeV]; R[cm]",
			1000,0.0,5.0,200,0.,30,AliDielectronVarManager::kM,AliDielectronVarManager::kR);

  histos->UserHistogram("MyPair","InvMyMassdPhi","dPhi of pairs vs Inv.MyMass;Inv. MyMass [GeV]; dPhi[rad]",
			1000,0.0,5.0,200,0.,3.1415,AliDielectronVarManager::kM,AliDielectronVarManager::kCosPointingAngle);

  histos->UserHistogram("MyPair","PhivR","R vs Phiv; R[cm]; Phiv[rad]",
			200,0.,30, 200,0,4, AliDielectronVarManager::kR,AliDielectronVarManager::kPsiPair);
  

  histos->UserHistogram("MyPairV0","InvV0MyMassPairPt","Inv.MyMass vs PairPt;Inv. MyMass [GeV]; pT [GeV]",
			1000,0.0,5.0,200,0.,20.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

  histos->UserHistogram("MyPairV0","InvV0MyMassOpeningAngle","Opening Angle vs Inv.MyMass;Inv. MyMass [GeV]; OA[rad]",
			1000,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("MyPairV0","InvV0MyMassPsiPair","PsiPair vs Inv.MyMass;Inv. MyMass [GeV]; Psi[rad]",
			1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPsiPair);

  histos->UserHistogram("MyPairV0","InvV0MyMassPhivPair","PhivPair vs Inv.MyMass;Inv. MyMass [GeV]; Phiv [rad]",
			1000,0.0,5.0,200,0.,4,AliDielectronVarManager::kM,AliDielectronVarManager::kPhivPair);

  histos->UserHistogram("MyPairV0","InvV0MyMassR","R vs Inv.MyMass;Inv. MyMass [GeV]; R[cm]",
			1000,0.0,5.0,200,0.,30,AliDielectronVarManager::kM,AliDielectronVarManager::kR);
  histos->UserHistogram("MyPairV0","PsiPairR","R vs PsiPair; R[cm]; PsiPair[rad]",
			200,0.,30, 200,0,4, AliDielectronVarManager::kR,AliDielectronVarManager::kPsiPair);
  histos->UserHistogram("MyPairV0","InvV0MyMassdPhi","dPhi of pairs vs Inv.MyMass;Inv. MyMass [GeV]; dPhi[rad]",
			1000,0.0,5.0,200,0.,3.1415,AliDielectronVarManager::kM,AliDielectronVarManager::kCosPointingAngle);
  
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////






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
  cf->AddVariable(AliDielectronVarManager::kPt,200,0,20);
  cf->AddVariable(AliDielectronVarManager::kM, 1000,0,5.00); //5Mev Steps
  cf->AddVariable(AliDielectronVarManager::kY,100,-2.,2.);
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);

  cf->AddVariable(AliDielectronVarManager::kCentrality,100, 0, 100);
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,320,0.,3.2);
  cf->AddVariable(AliDielectronVarManager::kPsiPair,320,0.,3.2);
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,200,0.,20.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPt,200,0.,20.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,1000,0.0.,1000.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,500,0.0.,500.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,100,-10.0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSnSigmaEle,100,-10.0.,10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kY,100,-2.,2.,kTRUE);
  //only in this case write MC truth info
  if (MCenabled) {
    cf->SetStepForMCtruth();
    //cf->SetStepsForMCtruthOnly();
    cf->SetStepForAfterAllCuts();
    cf->AddVariable(AliDielectronVarManager::kHaveSameMother,5,-2,2);
    //cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    //cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);



  /*
 AliDielectronSignalMC* lowMassDiele=new
    AliDielectronSignalMC("lowMassDiele","low mass dielectron pairs");
  lowMassDiele->SetLegPDGs(11,-11);
  lowMassDiele->SetCheckBothChargesLegs(kTRUE,kTRUE);
  lowMassDiele->SetLegSources(AliDielectronSignalMC::kPrimary,
      AliDielectronSignalMC::kPrimary);
  lowMassDiele->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(lowMassDiele);

  AliDielectronSignalMC* secondary=new
    AliDielectronSignalMC("secondary","secondary electrons pairs");
  secondary->SetLegPDGs(11,-11);
  secondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  secondary->SetLegSources(AliDielectronSignalMC::kSecondary,
      AliDielectronSignalMC::kSecondary);
  die->AddSignalMC(secondary);


  AliDielectronSignalMC* finalState=new
    AliDielectronSignalMC("finalState","finalState electrons pairs");
  finalState->SetLegPDGs(11,-11);
  finalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
  finalState->SetLegSources(AliDielectronSignalMC::kFinalState,
      AliDielectronSignalMC::kFinalState);
  die->AddSignalMC(finalState);
  */

}

//--------------------------------------
void EnableMC() {
  MCenabled=kTRUE;
}


//--------------------------------------
void SetSignals(AliDielectron *die)
{


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

  // physical backgrounds (electrons from other sources)
  AliDielectronSignalMC* diEleContinuum = new AliDielectronSignalMC("diEleContinuum","di-electron continuum");     // all di-electrons originating in the collision
  diEleContinuum->SetLegPDGs(11,-11);
  diEleContinuum->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleContinuum->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(diEleContinuum);
  
  AliDielectronSignalMC* diEleCharm = new AliDielectronSignalMC("diEleCharm","di-electrons from charm");  // dielectrons originating from charm hadrons (not neccessary from same mother)
  diEleCharm->SetLegPDGs(11,-11);
  diEleCharm->SetMotherPDGs(403,403);
  diEleCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleCharm);

  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("diEleOpenCharm","di-electrons from open charm");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharm);


  // background from secondary electrons
  AliDielectronSignalMC* secondaryElectrons = new AliDielectronSignalMC("secondaryElectrons","Secondary electrons");   // all di-electrons from secondary electrons (interaction with detector)
  secondaryElectrons->SetLegPDGs(11,-11);
  secondaryElectrons->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  secondaryElectrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(secondaryElectrons);

  /*
  AliDielectronSignalMC* primarySecElePairs = new AliDielectronSignalMC("primarySecElePairs","Primary+Secondary electron pairs");  // primary-secondary pairs
  primarySecElePairs->SetLegPDGs(11,-11);
  primarySecElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  primarySecElePairs->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  die->AddSignalMC(primarySecElePairs);
  */

  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");      // pairs made from conversion (may be also from 2 different conversions)
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  die->AddSignalMC(conversionElePairs);

  // misidentification
  /*
  AliDielectronSignalMC* allEleMisIdPairs = new AliDielectronSignalMC("allEleMisIdPairs","all electron+misid. pairs");  // one true electron + a mis-id electron (all sources included)
  allEleMisIdPairs->SetLegPDGs(11,11,kFALSE,kTRUE);
  allEleMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(allEleMisIdPairs);

  AliDielectronSignalMC* allMisIdMisIdPairs = new AliDielectronSignalMC("allMisIdMisIdPairs","all misid.+misid. pairs");  // mis-id + mis-id
  allMisIdMisIdPairs->SetLegPDGs(11,11,kTRUE,kTRUE);
  allMisIdMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(allMisIdMisIdPairs);

  AliDielectronSignalMC* elePionPairs = new AliDielectronSignalMC("elePionPairs","electron+pion pairs");    // true electron + mis-id pion
  elePionPairs->SetLegPDGs(11,211);
  elePionPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(elePionPairs);

  AliDielectronSignalMC* eleKaonPairs = new AliDielectronSignalMC("eleKaonPairs","electron+kaon pairs");   // true electron + mis-id kaon
  eleKaonPairs->SetLegPDGs(11,321);
  eleKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(eleKaonPairs);

  AliDielectronSignalMC* eleProtonPairs = new AliDielectronSignalMC("eleProtonPairs","Electron+proton pairs");  // true electron + mis-id proton
  eleProtonPairs->SetLegPDGs(11,2212);
  eleProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(eleProtonPairs);

  AliDielectronSignalMC* piPiPairs = new AliDielectronSignalMC("piPiPairs","pion+pion pairs");    // mis-id pion + mis-id pion
  piPiPairs->SetLegPDGs(211,211);
  piPiPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piPiPairs);

  AliDielectronSignalMC* piKaonPairs = new AliDielectronSignalMC("piKaonPairs","pion+kaon pairs");  // mis-id pion + mis-id kaon
  piKaonPairs->SetLegPDGs(211,321);
  piKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piKaonPairs);

  AliDielectronSignalMC* piProtonPairs = new AliDielectronSignalMC("piProtonPairs","pion+proton pairs");  // mis-id pion + mis-id proton
  piProtonPairs->SetLegPDGs(211,2212);
  piProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piProtonPairs);

  AliDielectronSignalMC* kaonKaonPairs = new AliDielectronSignalMC("kaonKaonPairs","kaon+kaon pairs");  // mis-id kaon + mis-id kaon
  kaonKaonPairs->SetLegPDGs(321,321);
  kaonKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(kaonKaonPairs);

  AliDielectronSignalMC* kaonProtonPairs = new AliDielectronSignalMC("kaonProtonPairs","kaon+proton pairs");   // mis-id kaon + mis-id proton
  kaonProtonPairs->SetLegPDGs(321,2212);
  kaonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(kaonProtonPairs);

  AliDielectronSignalMC* protonProtonPairs = new AliDielectronSignalMC("protonProtonPairs","proton+proton pairs");  // mis-id proton + mis-id proton
  protonProtonPairs->SetLegPDGs(2212,2212);
  protonProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(protonProtonPairs);

  AliDielectronSignalMC* muonAllPairs = new AliDielectronSignalMC("muonAllPairs","muon+everything pairs");        // mis-id muon + something else (electron, pion, kaon, proton)
  muonAllPairs->SetLegPDGs(13,13,kFALSE,kTRUE);
  muonAllPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(muonAllPairs);
  */


  AliDielectronSignalMC* pi0Sig = new AliDielectronSignalMC("pi0", "pi0Signal"); ///pi0 dalitz pairs 
  pi0Sig->SetLegPDGs(11,-11);
  pi0Sig->SetMotherPDGs(111,111);
  pi0Sig->SetMothersRelation(AliDielectronSignalMC::kSame);
  pi0Sig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pi0Sig->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  pi0Sig->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pi0Sig->SetCheckBothChargesMothers(kTRUE,kTRUE);
  pi0Sig->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(pi0Sig);


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

  /*
  AliDielectronSignalMC* convSig = new AliDielectronSignalMC("conv", "convSignal");
  convSig->SetLegPDGs(11,-11);
  convSig->SetMothersRelation(AliDielectronSignalMC::kSame);
  convSig->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); //in this case, ee needs to be true=IsPhysicalPrimary(label)
  convSig->SetCheckBothChargesLegs(kFALSE,kFALSE);
  die->AddSignalMC(convSig);
  */
}

