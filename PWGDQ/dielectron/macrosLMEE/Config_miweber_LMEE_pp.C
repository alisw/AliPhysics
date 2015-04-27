//#include "AliDielectron.h"
//#include "PWGDQ/dielectron/macrosLMEE/LMEECutLib.C"

void      InitHistograms(AliDielectron *die, Int_t cutDefinition);
void      InitCF(AliDielectron* die, Int_t cutDefinition);
TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD *GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx};

TString names=("pp2010All_pidITSTPCTOFif_trkSPDfirst_1;pp2010All_pidITSTPCTOFif_trkSPDfirsttight_1;pp2010Low_pidITSTPCTOFif_trkSPDfirst_1;pp2010Mid_pidITSTPCTOFif_trkSPDfirst_1;pp2010High_pidITSTPCTOFif_trkSPDfirst_1");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();
Bool_t MCenabled=kFALSE;//needed for LMEEcutlib
Bool_t isQAtask=kFALSE;

AliDielectron* Config_miweber_pp(Int_t cutDefinition, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  isQAtask=kFALSE;
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
  if (cutDefinition==0) {
    selectedPID = LMEECutLib::kpp2010All_pidITSTPCTOFif_trkSPDfirst_1;
    selectedCentrality = LMEECutLib::kpp2010All;
    die->SetPreFilterAllSigns();   
    rejectionStep=kTRUE;
    //LMcutlib->SetEtaCorrection(die, selectedPID, selectedCentrality, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  }
  else if (cutDefinition==1) {
    selectedPID = LMEECutLib::kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1;
    selectedCentrality = LMEECutLib::kpp2010All;
    die->SetPreFilterAllSigns();   
    rejectionStep=kTRUE;
    //LMcutlib->SetEtaCorrection(die, selectedPID, selectedCentrality, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  }
  else if (cutDefinition==2) {
    selectedPID = LMEECutLib::kpp2010All_pidITSTPCTOFif_trkSPDfirst_1;
    selectedCentrality = LMEECutLib::kpp2010Low;
    die->SetPreFilterAllSigns();   
    rejectionStep=kTRUE;
    //LMcutlib->SetEtaCorrection(die, selectedPID, selectedCentrality, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  }
  else if (cutDefinition==3) {
    selectedPID = LMEECutLib::kpp2010All_pidITSTPCTOFif_trkSPDfirst_1;
    selectedCentrality = LMEECutLib::kpp2010Mid;
    die->SetPreFilterAllSigns();   
    rejectionStep=kTRUE;
    //LMcutlib->SetEtaCorrection(die, selectedPID, selectedCentrality, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
  }
  else if (cutDefinition==4) {
    selectedPID = LMEECutLib::kpp2010All_pidITSTPCTOFif_trkSPDfirst_1;
    selectedCentrality = LMEECutLib::kpp2010High;
    die->SetPreFilterAllSigns();   
    rejectionStep=kTRUE;
    //LMcutlib->SetEtaCorrection(die, selectedPID, selectedCentrality, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kEta);
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
  //die->SetPreFilterAllSigns();   rejectionStep=kTRUE;
  //die->SetPreFilterUnlikeOnly(); rejectionStep=kTRUE;
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
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID, kFALSE) );
    
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
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID, kFALSE) );
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
  // InitCF(die,cutDefinition);
  
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
  histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","RunNumber","",AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600"),AliDielectronVarManager::kRunNumber);
  histos->UserHistogram("Event","Centrality","","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","centrality","",100,0,100,AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","nESDTracks","",8000,0,80000,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","Nacc","",8000,0,8000,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","RefMultTPConly","",8000,0,8000,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","epTPC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","",120,-TMath::PiOver2(),TMath::PiOver2(),120,-TMath::PiOver2(),TMath::PiOver2(),AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
  
  
  //add histograms to Track classes
  // axis labels are set to the values in 'AliDielectronVarManager.cxx' if the histogram title is empty or starts with ';'! [see AliDielectronHistos::AdaptNameTitle(...)]
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",
                        GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  
  // ITS
  histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                        GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
  histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
//histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
//                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                        GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                        GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  // 3D
  histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                        GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                        GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                        AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  
  if (isQAtask) {
    histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                          GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPC_dEdx_P_RunNumber",";p_{in} (GeV/c);TPC signal (arb units);run",
                          GetVector(kP2D), GetVector(kTPCdEdx), GetVector(kRuns), 
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal);
    
    //    histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
    //                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
    //                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    //    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
    //                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
    //                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_Nacc",";N_{TPC ref};n#sigma_{ele}^{TPC};N_{acc}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  }
  
  // TRD
  // TRD variables need lot of computing time. since a performance update by Julian, they will not be computed if not needed! (status 7.3.14)
  //  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;p_{in} (GeV/c);TRD prob Electrons",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle);
  //  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;p_{in} (GeV/c);TRD prob Pions",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio);
  
  // TOF
  histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
                        GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  if (isQAtask) {
    histos->UserHistogram("Track","TOFnSigmaPio_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P",";p_{in} (GeV/c);TOF number of sigmas Kaons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P",";p_{in} (GeV/c);TOF number of sigmas Protons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  }    
  // 2D-PID
  if (isQAtask) {
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                          50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                          50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);
  }
  
  // Eta and Phi
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
                        
  // DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  
  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  if (isQAtask) {
    histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                          160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
                          GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
                          GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                          100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                          120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
  }
  
  if (!isQAtask) 
  {
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","",240,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
    
    //2D and 3D histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                          100,-1.,1., 120,0.,TMath::TwoPi(),
                          AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                          GetVector(kMee), GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                          GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
    
    //opening angle and PhiV
    histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                          GetVector(kMee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                          GetVector(kMee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                          GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                          GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                          GetVector(kOpAng), GetVector(kPhiV), 
                          AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
    
    //centrality
    histos->UserHistogram("Pair","InvMass_Centrality",";Inv. Mass [GeV];Centrality;#pairs",
                          GetVector(kMee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
    histos->UserHistogram("Pair","PairPt_Centrality",";Pair Pt [GeV];Centrality;#pairs",
                          GetVector(kPtee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);
    
    
    //add histograms to Track classes
    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Pre","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
    
    histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                          GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                          GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    
    histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
    
    histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
  
  }// end: (!isQAtask)
  
  die->SetHistogramManager(histos);
}



TVectorD *GetVector(Int_t var) 
{
  switch (var) 
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
      
    case kSigmaEle:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
      else          return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else          return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(120,  0.,120.);
      else          return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);
      
    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 
                                                                   0.90, 0.94, 0.98, 1.02, 1.06, 
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00 
                                                                   ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50 
                                                                   ");
    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 
                                                                   5.00, 6.00, 7.00, 8.00 
                                                                   ");
    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 
                                                                   2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
                                                                   4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00 
                                                                   ");
      //2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 
      //2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 
      //3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 
      
  case kRuns:   return AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600");
    
    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  } 
  //if ( var.EqualTo("p_2D"      , kIgnoreCase) ) return AliDielectronHelper::MakeLinBinning(160,0.,8.);
}



TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  
  //  return vec;
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
