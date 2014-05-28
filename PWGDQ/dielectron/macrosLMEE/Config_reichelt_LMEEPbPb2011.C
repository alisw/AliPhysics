//#include "AliDielectron.h"
//#include "PWGDQ/dielectron/macrosLMEE/LMEECutLib.C"

void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
//adjust for 1,4,8
TString names=("ITSTPCTOFif_trkSPDfirst_1_kSemi;ITSTPCTOFif_trkSPDfirst5cls_4_kSemi;ITS2gevTPCTOFif_trkSPDfirst_5_tight_kSemi");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();
Bool_t MCenabled=kFALSE;//needed for LMEEcutlib

AliDielectron* Config_reichelt_LMEEPbPb2011(Int_t cutDefinition, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  Int_t selectedPID=-1;
  Int_t selectedCentrality=-1;
  Bool_t rejectionStep=kFALSE;
  //Bool_t doPairing=kTRUE;
  LMEECutLib*  LMcutlib = new LMEECutLib();
  
  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  
  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));
  //die->SetHasMC(hasMC);
  MCenabled=hasMC;
  
  cout << "cutDefinition = " << cutDefinition << endl;
  // Setup Analysis Selection
  if (cutDefinition==0) {           // Config for Technical Preliminaries for QM2014 (no prefilter used!)
    selectedPID = LMEECutLib::kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1;
    selectedCentrality = LMEECutLib::kPbPb2011SemiCentral;
  }
  else if (cutDefinition==1) {      // Config for systematic checks 1
    selectedPID = LMEECutLib::kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4;
    selectedCentrality = LMEECutLib::kPbPb2011SemiCentral;
  }
  else if (cutDefinition==2) {      // Config for systematic checks 2
    selectedPID = LMEECutLib::kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight;
    selectedCentrality = LMEECutLib::kPbPb2011SemiCentral;
  }
  else {
    cout << " =============================== " << endl;
    cout << " ==== INVALID CONFIGURATION ==== " << endl;
    cout << " cutDefinition = " << cutDefinition << endl;
    cout << " =============================== " << endl;
  }
  // __________________________________________________
  // POSSIBLE FURTHER SETTINGS
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // deactivate pairing to check track cuts or run with loose pid cuts:
  //die->SetNoPairing();
  // apply correct Pre-Filter Scheme, if necessary
  //die->SetPreFilterAllSigns();
  //die->SetPreFilterUnlikeOnly();
  // --------------------------------------------------
  
  //
  // Now configure task
  //
  // add centrality selection to event cuts
  die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedCentrality) );
  // switch off KF Particle
  die->SetUseKF(kFALSE);
  
  // --------------------------------------------------
  // with Rejection Step (Prefilter)
  // --------------------------------------------------
  if (rejectionStep) 
  {
    if (isESD) {
      die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna(selectedPID) );
      die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetESDTrackCutsAna(selectedPID) ); // this is redundant!?
    }
    // set initial track filter.
    // the function 'GetPIDCutsPre()' must also call 'GetTrackCutsPre()'!
	  die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsPre(selectedPID) );
    // set Prefilter. "remove all tracks from the Single track arrays that pass the cuts in this filter" (comment in AliDielectron.cxx)
    // cuts = REJECTION!!!
    die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) ); 
    // "apply leg cuts after the pre filter" (comment in AliDielectron.cxx)
    
    
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
	  die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
    
    
    // set Pairfilter.
    // cuts = SELECTION!!!
	  die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID,kFALSE) );
    
	}
  // --------------------------------------------------
  // without Rejection Step
  // --------------------------------------------------
	else 
  {
	  if (isESD) {
      die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna(selectedPID) );
	  }
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
	  die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
	  die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID,kFALSE) );
	}
  // --------------------------------------------------
  
  
  AliDielectronTrackRotator *rot= 0x0;
  //To save time and as it is not 100% test, rotation switched off
  /*AliDielectronTrackRotator *rot= LMcutlib->GetTrackRotator(selectedPID);
   die->SetTrackRotator(rot);
   */
  AliDielectronMixingHandler *mix=LMcutlib->GetMixingHandler(selectedPID);
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
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack");
  
  //Event class
  histos->AddClass("Event"); // all classes will be stored in 'THashList fHistoList'
  
  //Track classes
  //to fill also track info from 2nd event loop until 3
  // in AliDielectron.cxx: fgkTrackClassNames[4] = {"ev1+","ev1-","ev2+","ev2-"};
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  // fgkPairClassNames[11] = {
  //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
  //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
  //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
  //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
  // };
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    // Legs of final Pairs. Both charges together. No duplicate entries.
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
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
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
    // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
    histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
  }
  
  /*
   //track rotation
   
   histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
   histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
   */
	//add histograms to event class
	histos->UserHistogram("Event","nEvents","Number of processed events after cuts;#events",1,0.,1.,AliDielectronVarManager::kNevents);
	histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","centrality","N events vs centrality;centrality [%];#events",100,0,100,AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","nESDTracks","ESD tracks;ESD tracks;#events",1000,0,10000,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","Nacc","Number of accepted tracks;Nacc;#events",1200,0,1200,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","epTPC","TPC event plane angle (uncorr);EP angle TPC (uc);#events",320,-3.2,3.2,AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","V0AC event plane angle;EP angle V0AC;#events",320,-3.2,3.2,AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","event plane angle V0AC vs TPC;EP angle TPC (uc);EP angle V0AC",
                        160,-1.6,1.6,160,-1.6,1.6,AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
  
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px","Px;Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py","Py;Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz","Pz;Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn","TPC inner P vs P; P [GeV]; TPC inner P [GeV]",
                        160,0.,8.,160,0.,8.,AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  
  // ITS
  Bool_t makeLogx = kFALSE;
  histos->UserHistogram("Track","ITS_dEdx_P","ITS dEdx;P [GeV];ITS signal (arb units)",
                        160,0.,8.,700,0.,700.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal,makeLogx);
  histos->UserHistogram("Track","ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV];ITS number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  histos->UserHistogram("Track","ITSnSigmaPio_P","ITS number of sigmas Pions;P [GeV];ITS number of sigmas Pions",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio,makeLogx);
histos->UserHistogram("Track","ITSnSigmaKao_P","ITS number of sigmas Kaons;P [GeV];ITS number of sigmas Kaons",
                      160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao,makeLogx);
//histos->UserHistogram("Track","ITSnSigmaPro_P","ITS number of sigmas Protons;P [GeV];ITS number of sigmas Protons",
//                      160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro,makeLogx);
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P","TPC dEdx;P [GeV];TPC signal (arb units)",
                        160,0.,8.,120,0.,120.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,makeLogx);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;TPC signal (arb units)",
                        80,0.,4.,80,-4.,4.,50,50.,100.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal,makeLogx);
  histos->UserHistogram("Track","TPC_dEdx_P_run","TPC dEdx;P [GeV];TPC signal (arb units);run number",
                        AliDielectronHelper::MakeLinBinning(80,0.,4.), 
                        AliDielectronHelper::MakeLinBinning(50,50.,100.), 
                        AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600"), 
                        AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
  //                        AliDielectronHelper::MakeArbitraryBinning("170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 167988, 167987, 167800"), 
  
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,makeLogx);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,makeLogx);
  histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,makeLogx);
  // TRD
  // TRD variables need lot of computing time. since a performance update by Julian, they will not be computed if not needed! (status 7.3.14)
  //  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons",
  //                        160,0.,8.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,makeLogx);
  //  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions",
  //                        160,0.,8.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,makeLogx);
  // TOF
  histos->UserHistogram("Track","TOFnSigmaEle_P","TOF number of sigmas Electrons;P [GeV];TOF number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
  histos->UserHistogram("Track","TOFnSigmaPio_P","TOF number of sigmas Pions;P [GeV];TOF number of sigmas Pions",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio,makeLogx);
  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,makeLogx);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,makeLogx);
  histos->UserHistogram("Track","TOFbeta","TOF beta;P [GeV];TOF beta",
                        160,0.,8.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,makeLogx);
  // 2D-PID
  histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle","2D PID - TPC vs ITS;P [GeV];TPC number of sigmas Electrons;ITS number of sigmas Electrons",
                        50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                        AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle","2D PID - TPC vs TOF;P [GeV];TPC number of sigmas Electrons;TOF number of sigmas Electrons",
                        50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                        AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
  
  // Eta and Phi
  histos->UserHistogram("Track","Eta","Eta;Eta;#tracks",
                        200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi;Phi;#tracks",
                        320,0.,6.4,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map;Eta;Phi",
                        100,-1,1,320,0,6.4,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
                        
  histos->UserHistogram("Track","TPC_dEdx_Eta","TPC dEdx;Eta;TPC signal (arb units)",
                        100,-1,1,120,0.,120.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPC_dEdx_Eta_P","TPC dEdx;Eta;TPC signal (arb units); TPC inner P [GeV]",
                        100,-1,1,60,0.,120.,80,0.,4.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons",
                        100,-1,1,100,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta_P","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons; TPC inner P [GeV]",
                        100,-1,1,80,-4.,4.,80,0.,4.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","TPCnSigmaKao_Eta","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons",
                        100,-1,1,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  histos->UserHistogram("Track","TPCnSigmaKao_Eta_P","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons; TPC inner P [GeV]",
                        50,-1,1,100,-10.,10.,80,0.,4.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao,AliDielectronVarManager::kPIn);
  
  // DCA
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",
                        200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",
                        400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","dXY dZ Map;dXY;dZ",
                        100,-1.,1.,300,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  
  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC","Fraction of shared clusters assigned in the TPC;TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff","TPC cluster difference;TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN","Number of PID Clusters TPC;TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  histos->UserHistogram("Track","TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters;TPC number clusters;TPC crossed rows",
                        160,-0.5,159.5,160,-0.5,159.5,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCcrossedRows_Pt","TPC crossed rows vs Pt;Pt [GeV];TPC crossed rows",
                        160,0.,8.,160,-0.5,159.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt","Number of Crossed Rows TPC over Findable vs Pt;Pt [GeV];TPC crossed rows over findable",
                        160,0.,8.,120,0.,1.2,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta","Number of Crossed Rows TPC over Findable vs Eta;Eta;TPC crossed rows over findable",
                        100,-1,1,120,0.,1.2,AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi","Number of Crossed Rows TPC over Findable vs Phi;Phi;TPC crossed rows over findable",
                        320,0.,6.4,120,0.,1.2,AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
  
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        500,0.,5.,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","PairPt","PairPt;Pair Pt [GeV];#pairs",
                        160,0.,8.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        200,-2.,2.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle;#pairs",
                        320,0.,3.2,AliDielectronVarManager::kOpeningAngle);
  
  //2D and 3D histograms
  histos->UserHistogram("Pair","InvMass_PairPt","PairPt vs InvMass;Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                        500,0.,5., 160,0.,8.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Eta_Phi_Pair","Phi vs Eta (pair);Eta;Phi;#pairs",
                        100,-1.,1., 320,0.,6.4,
                        AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair","InvMass:PairPt:PhivPair;Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        150,0.,1.5, 40,0.,4., 64,0.,3.2,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle","InvMass:PairPt:OpeningAngle;Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                        150,0.,1.5, 40,0.,4., 64,0.,3.2,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle","InvMass:PhivPair:OpeningAngle;Inv. Mass [GeV];PhiV;Opening Angle",
                        50,0.,0.5, 160,0.,3.2, 40,0.,0.8,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
  
  //opening angle and PhiV
  histos->UserHistogram("Pair","InvMass_OpeningAngle","Opening Angle vs InvMass;Inv. Mass [GeV];Opening Angle;#pairs",
                        500,0.,5., 160,0.,3.2,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_PhivPair","PhiV vs InvMass;Inv. Mass [GeV];PhiV;#pairs",
                        500,0.,5., 160,0.,3.2,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","PairPt_OpeningAngle","Opening Angle vs PairPt;Pair Pt [GeV];Opening Angle;#pairs",
                        160,0.,8., 160,0.,3.2,
                        AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","PairPt_PhivPair","PhiV vs PairPt;Pair Pt [GeV];PhiV;#pairs",
                        160,0.,8., 160,0.,3.2,
                        AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  histos->UserHistogram("Pair","OpeningAngle_PhivPair","PhiV vs OpeningAngle;Opening Angle;PhiV;#pairs",
                        160,0.,3.2, 160,0.,3.2,
                        AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
  
  //centrality
  histos->UserHistogram("Pair","InvMass_Centrality","Centrality vs InvMass;Inv. Mass [GeV];Centrality;#pairs",
                        500,0.,5., 102,-1,101,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Pair","PairPt_Centrality","Centrality vs PairPt;Pair Pt [GeV];Centrality;#pairs",
                        160,0.,8., 102,-1,101,
                        AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);
  
  
/*  
  //add histograms to Track classes
  histos->UserHistogram("Pre","Pt","Pt;Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Pre","Px","Px;Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Pre","Py","Py;Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Pre","Pz","Pz;Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  
  histos->UserHistogram("Pre","ITS_dEdx_P","ITS dEdx;P [GeV];ITS signal (arb units)",
                        160,0.,8.,700,0.,700.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal,makeLogx);
  histos->UserHistogram("Pre","TPC_dEdx_P","TPC dEdx;P [GeV];TPC signal (arb units)",
                        160,0.,8.,120,0.,120.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,makeLogx);
  
  histos->UserHistogram("Pre","ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV];ITS number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  //histos->UserHistogram("Pre","ITSnSigmaPio_P","ITS number of sigmas Pions;P [GeV];ITS number of sigmas Pions",
  //                      160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio,makeLogx);
  //histos->UserHistogram("Pre","ITSnSigmaKao_P","ITS number of sigmas Kaons;P [GeV];ITS number of sigmas Kaons",
  //                      160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao,makeLogx);
  //histos->UserHistogram("Pre","ITSnSigmaPro_P","ITS number of sigmas Protons;P [GeV];ITS number of sigmas Protons",
  //                      160,0.,8.,200,-10.,10.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro,makeLogx);
  
  histos->UserHistogram("Pre","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons",
                        160,0.,8.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  histos->UserHistogram("Pre","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,makeLogx);
  histos->UserHistogram("Pre","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,makeLogx);
  histos->UserHistogram("Pre","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons",
                        160,0.,8.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,makeLogx);
*/  
  
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,30.,50.,80.,100.");
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,160,0.,8.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,350,0.,700.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,60,0.,120.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
  
  //only in this case write MC truth info
  if (MCenabled) { // more elegant: die->GetHasMC() 
    cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }
  
  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}
