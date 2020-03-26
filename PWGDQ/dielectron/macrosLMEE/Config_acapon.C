// ROOT6 modifications
#ifdef __CLING__
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <PWGDQ/dielectron/macrosLMEE/LMEECutLib_acapon.C>
#endif
void InitHistograms(AliDielectron *die, Bool_t doPairing, Bool_t trackVarPlots, Int_t whichDetPlots, Bool_t v0plots, Bool_t plots3D, TString cutDefinition);
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee = 0, kPtee, kPt, kRuns, kPhiV, kOpAng, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx, kCent, kPhi2D};

AliDielectron* Config_acapon(TString cutDefinition,
                             Bool_t hasMC,
                             Bool_t SDDstatus,
                             Bool_t doPairing,
                             Bool_t applyPairCuts,
                             Bool_t doEventMixing,
                             Bool_t trackVarPlots,
                             Int_t whichDetPlots,
                             Bool_t v0plots,
                             Bool_t setITScorr,
                             Bool_t setTPCcorr,
                             Bool_t setTOFcorr,
                             Bool_t plots3D)
{
  // Setup the instance of AliDielectron
  LMEECutLib*  LMcutlib = new LMEECutLib(SDDstatus);

  // Init AliDielectron
  AliDielectron* die = new AliDielectron(Form("%s",cutDefinition.Data()), Form("AliDielectron with cuts: %s",cutDefinition.Data()));

  if(setTPCcorr && !hasMC){
    LMcutlib->SetEtaCorrectionTPC(die, AliDielectronVarManager::kP,
                                       AliDielectronVarManager::kEta,
                                       AliDielectronVarManager::kRefMultTPConly);
  }
  if(setITScorr){
    LMcutlib->SetEtaCorrectionITS(die, AliDielectronVarManager::kP,
                                       AliDielectronVarManager::kEta,
                                       AliDielectronVarManager::kRefMultTPConly, hasMC);

  }
  if(setTOFcorr){
    LMcutlib->SetEtaCorrectionTOF(die, AliDielectronVarManager::kP,
                                       AliDielectronVarManager::kEta,
                                       AliDielectronVarManager::kRefMultTPConly, hasMC);

  }

  // Deactivate pairing to check track cuts or run with loose pid cuts:
  if(!doPairing){
    die->SetNoPairing();
  }
  if(hasMC){
    die->SetHasMC(hasMC);
  }

  // Event mixing handler. Will be set after cut sets are set up due to flag
  // described below
  AliDielectronMixingHandler* mix = 0x0;
  // One "standard" setting used for mixing unless doing specific mixing tests
  // Flag will be switched if one of those cut sets are chosen
  Bool_t nonStandardMixing = kFALSE;

  die->SetPreFilterUnlikeOnly(kTRUE);

  std::cout << "cutDefinition = " << cutDefinition << std::endl;
  // ######### QA CUTS ##############
  // Simple cuts to trim outliers
  if(cutDefinition == "kAll"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kAllSpecies, LMEECutLib::kAllSpecies));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts( LMcutlib->GetPairCuts(LMEECutLib::kAllSpecies) );
    }
  } // Used for basic QA of data sets (simple cuts and PID)
  else if(cutDefinition == "kElectrons"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kElectrons, LMEECutLib::kElectrons));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kElectrons));
    }
  }

  // ##### CUTS USED TO CREATE TTREES ##########
  else if(cutDefinition == "kTTreeCuts"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kTTreeCuts, LMEECutLib::kTTreeCuts));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kTTreeCuts));
    }
  }

  // ######### CUTS TO OBTAIN CORRECTION MAPS FOR DATA AND MC  ###############
  else if(cutDefinition == "kV0_TPCcorr"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kV0_trackCuts, LMEECutLib::kV0_TPCcorr));
  }
  else if(cutDefinition == "kV0_ITScorr"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kV0_trackCuts, LMEECutLib::kV0_ITScorr));
  }
  else if(cutDefinition == "kV0_TOFcorr"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kV0_trackCuts, LMEECutLib::kV0_TOFcorr));
  }
  else if(cutDefinition == "kMCpdgSel"){
    die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCuts(LMEECutLib::kMCsel, LMEECutLib::kPdgSel) );
  }
  // Standard analysis withSDD analysis cuts and MVA ePID
  else if(cutDefinition == "kCutSet1"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // Standard noSDD analysis cuts using MVA ePID(req. ITS PID)
  else if(cutDefinition == "kCutSet3"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kScheidCuts, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // Standard analysis cuts with hadron rejection ePID
  else if(cutDefinition == "kScheidCuts"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kScheidCuts, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // Standard analysis cuts with 7TeV diElec ePID scheme
  else if(cutDefinition == "k7TeVpaper"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::k7TeVtrack, LMEECutLib::k7TeVPID));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // ######## Traditional Cut Set #################
  // Standard PID cut set taken from a Run 1 analysis
  else if(cutDefinition == "kTheoPID"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kTheoPID));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // Standard run1 track+PID cuts. Use V0 finder as well as conversion cuts
  else if(cutDefinition == "kTheoPIDv0finder"){
    // Applies very loose track cuts and no PID
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kCutSet1));
    die->GetTrackFilter().AddCuts(LMcutlib->GetV0finder()); // Dummy argument
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kTOFreq"){ // Copy of TheoPID cut setting however TOF always required
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kTOFreq));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  // Two cut settings to check PID efficiency using V0 electrons
  // (does not work for MC, checked 2019.05.08)
  else if(cutDefinition == "kV0_TTreeCutPID"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kV0_trackCuts, LMEECutLib::kTTreeCuts));
  }
  else if(cutDefinition == "kV0_MVAePID"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kV0_trackCuts, LMEECutLib::kMVA1));
  }
  // ######## Different R factor bin mixing schemes #################
  else if(cutDefinition == "kMixScheme1"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
    nonStandardMixing = kTRUE;
    mix = LMcutlib->GetMixingHandler(LMEECutLib::kMixScheme1);
  }
  else if(cutDefinition == "kMixScheme2"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
    nonStandardMixing = kTRUE;
    mix = LMcutlib->GetMixingHandler(LMEECutLib::kMixScheme2);
  }
  else if(cutDefinition == "kMixScheme3"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
    nonStandardMixing = kTRUE;
    mix = LMcutlib->GetMixingHandler(LMEECutLib::kMixScheme3);
  }
  else if(cutDefinition == "kMixScheme4"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
    nonStandardMixing = kTRUE;
    mix = LMcutlib->GetMixingHandler(LMEECutLib::kMixScheme4);
  }
  else if(cutDefinition == "kMixScheme5"){
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kCutSet1, LMEECutLib::kMVA1));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
    nonStandardMixing = kTRUE;
    mix = LMcutlib->GetMixingHandler(LMEECutLib::kMixScheme5);
  }
  // Produces plots using MCtruth information to select dielectron pairs
  else if(cutDefinition == "kDCAdists"){
    // Applies very loose track cuts and standard PID
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kTTreeCuts, LMEECutLib::kTheoPID));
    LMcutlib->SetSignalsMC(die);
  }
  // Cut sets to to vary fITSshared cut
  else if(cutDefinition == "kITSshared1"){ // One shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared1, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kITSshared2"){ // Two shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared2, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kITSshared3"){ // Three shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared3, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kITSshared4"){ // Four shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared4, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kITSshared5"){ // Five shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared5, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else if(cutDefinition == "kITSshared6"){ // Six shared hits in ITS
    die->GetTrackFilter().AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kITSshared6, LMEECutLib::kHadRej));
    if(applyPairCuts){
      die->GetPairFilter().AddCuts(LMcutlib->GetPairCuts(LMEECutLib::kCutSet1));
    }
  }
  else{
    cout << " =============================== " << endl;
    cout << " ==== INVALID CONFIGURATION ==== " << endl;
    cout << " cutDefinition = " << cutDefinition << endl;
    cout << " =============================== " << endl;
  }

  // KF (whatever that means) is depreceated and will return incorrect results
  // The default setting is on though because......yep.....
  die->SetUseKF(kFALSE);

  if(doEventMixing){
    if(!nonStandardMixing){
      mix = LMcutlib->GetMixingHandler(LMEECutLib::kCutSet1);
    }
    die->SetMixingHandler(mix);
  }

  InitHistograms(die, doPairing, trackVarPlots, whichDetPlots, v0plots, plots3D, cutDefinition);

  return die;
}

//______________________________________________________________________________________

void InitHistograms(AliDielectron *die, Bool_t doPairing, Bool_t trackVarPlots, Int_t whichDetPlots, Bool_t v0plots, Bool_t plots3D, TString cutDefinition){

    // Setup histogram Manager
    AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(),die->GetTitle());

    // Initialise histogram classes
    histos->SetReservedWords("Track;Pair");//;Track_Legs");//RejPair;RejTrack");

    // Event class
    histos->AddClass("Event");

    // Track classes
    //0,1: +- ev1, 2,3: +- ev2
    for(Int_t i = 0; i < 2; ++i){
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
    if(doPairing){
      for(Int_t i = 0; i < 3; ++i){
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      }
      // Mixed event pairs if mixing handler present
      if(die->GetMixingHandler()){
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
      }
    }

    // Add MC signal histograms
    if(die->GetMCSignals()){
      for(Int_t i = 0; i < die->GetMCSignals()->GetEntriesFast(); ++i){
        histos->AddClass(Form("Track_%s_%s", AliDielectron::PairClassName(1), die->GetMCSignals()->At(i)->GetName()));
        if(doPairing){
          histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
        }
      }
    }

    TH1::AddDirectory(kFALSE);
    // Add histograms to event class
    histos->UserHistogram("Event", "nEvents",        "", 1,   0.,  1.,  AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event", "Centrality",     "", 100, 0,   100, AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event", "nESDTracks",     "", 500, 0,   500, AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event", "zVertexPrimary", "", 122, -11, 11,  AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;N of Vertex Contributors;N of events",
                          200, -0.5, 199.5, AliDielectronVarManager::kNVtxContrib);
    //------ Num. tracks -----/
    histos->UserHistogram("Event", "Accepted tracks", "", 50,  0,    50,   AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event", "Ntracks",         "", 100, -0.5, 99.5, AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NtracksVsVtxZ","", 150, -15, 15, 50, -0.5, 49.5, AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event", "RefMultTPConly", "", 300, 0, 300, AliDielectronVarManager::kRefMultTPConly);
    //------ Pile-up check plots
    histos->UserHistogram("Event", "SPD clusters vs tracklets", "",
                          75, 0, 150, 30, 0, 60, AliDielectronVarManager::kNaccTrcklts10, AliDielectronVarManager::kITSLayerFirstCls);
    histos->UserHistogram("Event","NTPCcls",   "",
                          500, 0, 1000, AliDielectronVarManager::kNTPCclsEvent);
    histos->UserHistogram("Event","NTPCtrkswITSout","",
                          500, 0, 1000, AliDielectronVarManager::kNTPCtrkswITSout);
    // Using new centrality estimator (run2 V0M)
    histos->UserHistogram("Event","NTPCclsEventRun2",   "kNTPCclsEvent;Centrality/%;kNTPCclsEvent",
                          202, -1., 100., 500, 0, 1000, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);
    histos->UserHistogram("Event","NTPCtrkswITSoutRun2","kNTPCtrkswITSout;Centrality/%;kNTPCtrkswITSout",
                          202, -1., 100., 500, 0, 1000, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCtrkswITSout);
    // Using new centrality estimator (run1 V0M)
    histos->UserHistogram("Event","NTPCclsEventRun1",   "kNTPCclsEvent;Centrality/%;kNTPCclsEvent",
                          202, -1., 100., 500, 0, 1000, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kNTPCclsEvent);
    histos->UserHistogram("Event","NTPCtrkswITSoutRun1","kNTPCtrkswITSout;Centrality/%;kNTPCtrkswITSout",
                          202, -1., 100., 500, 0, 1000, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kNTPCtrkswITSout);

    //--------- V0 plots ------------------------//
    histos->UserHistogram("Event","MultV0","Multiplicity V0;V0M amplitude",                       4000, -0.5, 3999.5, AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","EqMultV0","Equalized Multiplicity V0;Equalized V0M amplitude", 4000, -0.5, 3999.5, AliDielectronVarManager::kEqMultV0);
    histos->UserHistogram("Event","ChMultV0","Charged Multiplicity V0;Charged V0M amplitude",     1000, -0.5, 999.5,  AliDielectronVarManager::kVZEROchMult);
    histos->UserHistogram("Event","CentralityV0Mrun2","Centrality V0;V0M percentile",   102, -1, 101, AliDielectronVarManager::kCentralityNew); // V0M in run2
    histos->UserHistogram("Event","CentralityV0Mrun1","Centrality V0;V0M percentile",   102, -1, 101, AliDielectronVarManager::kCentrality); // V0M in run1
    histos->UserHistogram("Event","CentralityV0A","Centrality V0;V0A percentile",   102, -1, 101, AliDielectronVarManager::kCentralityV0A);
    histos->UserHistogram("Event","CentralityV0C","Centrality V0;V0C percentile",   102, -1, 101, AliDielectronVarManager::kCentralityV0C);
    histos->UserHistogram("Event","CentralityZNA","Centrality V0;V0ZNA percentile", 102, -1, 101, AliDielectronVarManager::kCentralityZNA);
    histos->UserHistogram("Event","CentralitySPD","Centrality V0;V0SPD percentile", 102, -1, 101, AliDielectronVarManager::kCentralitySPD);
    histos->UserHistogram("Event","CentralityCL0","Centrality V0;CL0 percentile",   102, -1, 101, AliDielectronVarManager::kCentralityCL0);
    histos->UserHistogram("Event","CentralityCL1","Centrality V0;CL1 percentile",   102, -1, 101, AliDielectronVarManager::kCentralityCL1);

    // 2D V0 plots
    histos->UserHistogram("Event","V0AvsV0C","V0A;V0C",
                          102, -1, 101, 102, -1, 101, AliDielectronVarManager::kCentralityV0C, AliDielectronVarManager::kCentralityV0A);
    histos->UserHistogram("Event","V0MvsV0C","V0M;V0C",
                          102, -1, 101, 102, -1, 101, AliDielectronVarManager::kCentralityV0C, AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","V0MvsV0A","V0M;V0A",
                          102, -1, 101, 102, -1, 101, AliDielectronVarManager::kCentralityV0A, AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","RefMulTPConlytVsMult","#Charged Tracks Multiplicity (%);Ref. Mult TPC only",
                          100, 0, 100, 600, 0, 600, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kRefMultTPConly);
    histos->UserHistogram("Event","RefMulOvRefMultTPConlytVsMult","#Charged Tracks Multiplicity (%);Ref. Mult Over Ref Mult TPC only",
                          100, 0, 100, 600, 0, 600, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kRefMultOvRefMultTPConly);
    histos->UserHistogram("Event","NumTrackletsVsMult05","#Charged Tracks Multiplicity (%), #eta < |0.5|;Num. SPD tracklets ",
                          100, 0, 100, 600, 0, 600, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNaccTrckltsEsd05);
    histos->UserHistogram("Event","NumTrackletsVsMult10","#Charged Tracks Multiplicity (%), #eta < |1|;Num. SPD tracklets ",
                          100, 0, 100, 600, 0, 600, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNaccTrckltsEsd10);

    //------Pile up check------
    histos->UserHistogram("Event","SPDclustsVsSPDtracklets","SPD: clusters vs tracklets #eta < |1|; Tracklets; Cluster", 600, 0, 600, 1000, 0, 1000,
                          AliDielectronVarManager::kNaccTrcklts10, AliDielectronVarManager::kITSLayerFirstCls);

    if(trackVarPlots){
      // Add histograms to Track classes
      histos->UserHistogram("Track","Pt",";Pt (GeV);#tracks",200,0,10.,AliDielectronVarManager::kPt);
      histos->UserHistogram("Track","Px",";Px (GeV);#tracks",200,0,10.,AliDielectronVarManager::kPx);
      histos->UserHistogram("Track","Py",";Py (GeV);#tracks",200,0,10.,AliDielectronVarManager::kPy);
      histos->UserHistogram("Track","Pz",";Pz (GeV);#tracks",200,0,10.,AliDielectronVarManager::kPz);
      histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",200,0,10,AliDielectronVarManager::kPIn);
      // Eta and Phi
      histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
      histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
      histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
      // DCA
      histos->UserHistogram("Track","dXY","",400,-2.,2.,AliDielectronVarManager::kImpactParXY);
      histos->UserHistogram("Track","dZ" ,"",600,-4.,4.,AliDielectronVarManager::kImpactParZ);
      histos->UserHistogram("Track","dXYsig","",100,0,20,AliDielectronVarManager::kImpactParXYsigma);
      histos->UserHistogram("Track","dZsig" ,"",100,0,20,AliDielectronVarManager::kImpactParZsigma);
      histos->UserHistogram("Track","SPD clusters vs. tracklets",";tracklets;SPD clusters",
                            150,0,150,6,0,6,AliDielectronVarManager::kNTrk ,AliDielectronVarManager::kITSLayerFirstCls);
      /* histos->UserHistogram("Track","DCA_{xy} vs p_T","",300,0,5,100,0,0.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY); */
      /* histos->UserHistogram("Track","DCA_{Z} vs p_T","",300,0,5,100,-1,1,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY); */

      // Track cut variables for trackQA
      // ITS
      histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsITS);
      histos->UserHistogram("Track","ITSnClsClusterMap",";ITS cluster map;#tracks",100, 0.0, 100.0,AliDielectronVarManager::kITSclusterMap);
      histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",110,0.,11.,AliDielectronVarManager::kITSchi2Cl);
      histos->UserHistogram("Track","nITSshared","#shared ITS clusters", 7, 0, 7, AliDielectronVarManager::kNclsSITS);
      histos->UserHistogram("Track","fracITSshared","frac. shared ITS clusters", 120, 0,  1.2, AliDielectronVarManager::kNclsSITS);

      // TPC
      histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNclsTPC);
      histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
      histos->UserHistogram("Track","nTPCshared","#shared TPC clusters", 170, 0, 170, AliDielectronVarManager::kNclsSTPC);
      histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
      histos->UserHistogram("Track","TPCnCrossed",";TPC findable clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNFclsTPCr);
      histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",240,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
      histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                            160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
      histos->UserHistogram("Track","TPCcrossedRows_P",";Pt (GeV);TPC crossed rows",
                            GetVector(kPt), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kP,AliDielectronVarManager::kNFclsTPCr);
      histos->UserHistogram("Track","TPCcrossedRowsOverFindable_P",";P (GeV);TPC crossed rows over findable",
                            GetVector(kPt), BinsToVector(120,0.,1.2), AliDielectronVarManager::kP,AliDielectronVarManager::kNFclsTPCfCross);
      histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                            100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
      histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                            120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);

      // 1D PID plots for runByRun QA
      histos->UserHistogram("Track","nSigITSeRaw",";raw n#sigma^{ITS}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEleRaw);
      histos->UserHistogram("Track","nSigITSe",";n#sigma^{ITS}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","nSigTPCeRaw",";raw n#sigma^{TPC}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEleRaw);
      histos->UserHistogram("Track","nSigTPCe",";n#sigma^{TPC}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","nSigTOFeRaw",";raw n#sigma^{TOF}_{e};#tracks", GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEleRaw);
      histos->UserHistogram("Track","nSigTOFe",";n#sigma^{TOF}_{e};#tracks",        GetVector(kSigmaEle),   AliDielectronVarManager::kTOFnSigmaEle);
      histos->UserHistogram("Track","nSigTPCpi",";n#sigma^{TPC}_{#pi};#tracks",     GetVector(kSigmaOther), AliDielectronVarManager::kTPCnSigmaPio);

    }

    // ITS
    if((whichDetPlots & 1) == 1){
      histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                            GetVector(kPt), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
      histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","ITSnSigmaEleRaw_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEleRaw);
      histos->UserHistogram("Track","ITSnSigmaEle_Eta",";Eta;n#sigma_{ele}^{ITS}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","ITSnSigmaEle_Phi",";Phi;n#sigma_{ele}^{ITS}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","ITSnSigmaEle_Cent", ";Centrality;n#sigma_{ele}^{ITS}",
                            GetVector(kCent), GetVector(kSigmaEle), AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","ITSnSigmaEle_RunNumber",";run;n#sigma_{ele}^{ITS}",
                            GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kITSnSigmaEle);
    }
    // TPC
    if((whichDetPlots & 2) == 2){
      histos->UserHistogram("Track","TPC_dEdx_P",";p (GeV/c);TPC signal (arb units)",
                            GetVector(kPt), GetVector(kTPCdEdx), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCsignal);
      histos->UserHistogram("Track","TPCnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{TPC}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEleRaw_P",";p (GeV/c);n#sigma_{ele}^{TPC}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEleRaw);
      histos->UserHistogram("Track","TPCnSigmaEle_Pt",";p_{T} (GeV/c);n#sigma_{ele}^{TPC}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEle_Phi",";Phi;n#sigma_{ele}^{TPC}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                            GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
      histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                            GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
      histos->UserHistogram("Track","TPCnSigmaEle_Cent", ";Centrality;n#sigma_{ele}^{ITS}",
                            GetVector(kCent), GetVector(kSigmaEle), AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kTPCnSigmaEle);
    }
    // TOF
    if((whichDetPlots & 3) == 3){
      histos->UserHistogram("Track","TOFbeta_P",";p (GeV/c);TOF beta",
                            GetVector(kPt), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta);
      histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{elec}^{TOF}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle);
      histos->UserHistogram("Track","TOFnSigmaEleRaw_P",";p_{in} (GeV/c);n#sigma_{elec}^{TOF}",
                            GetVector(kPt), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEleRaw);
      histos->UserHistogram("Track","TOFnSigmaEle_Eta",";Eta;n#sigma_{elec}^{TOF}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);
      histos->UserHistogram("Track","TOFnSigmaEle_Phi",";Phi;n#sigma_{elec}^{TOF}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);
      histos->UserHistogram("Track","TOFnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TOF}",
                            GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTOFnSigmaEle);
      histos->UserHistogram("Track","TOFnSigmaEle_Cent", ";Centrality;n#sigma_{ele}^{ITS}",
                            GetVector(kCent), GetVector(kSigmaEle), AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kTOFnSigmaEle);
    }
    // 3D plots
    if(plots3D){
      histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                            50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                            AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
      histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                            50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                            AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);

      histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                            GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kPt),
                            AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
      histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                            GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                            AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
      histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                            GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                            AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
      histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
                            GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
                            AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
      histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
                            AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
      histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_Nacc",";N_{TPC ref};n#sigma_{ele}^{TPC};N_{acc}",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                            AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    }

    // Histograms for POST PID calibration
    // Define the min/max limits for each of the four variables:
    // P, numTrack, {DET}nSigma{ParticleType}, eta
    if(v0plots){
      const Int_t dimensions    = 4;
      Int_t bins[dimensions]    = {100, 40, 40, 16};
      Double_t xmin[dimensions] = {0., 0., -4, -0.8};
      Double_t xmax[dimensions] = {10., 800., 4, 0.8};
      // Define the histograms to be plotted using refMultTPConly
      UInt_t value_refMultTPC_ITSnSigmaEle[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kITSnSigmaEle, AliDielectronVarManager::kEta};
      UInt_t value_refMultTPC_TPCnSigmaEle[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kEta};
      UInt_t value_refMultTPC_TOFnSigmaEle[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTOFnSigmaEle, AliDielectronVarManager::kEta};
      histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value_refMultTPC_ITSnSigmaEle);
      histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value_refMultTPC_TPCnSigmaEle);
      histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value_refMultTPC_TOFnSigmaEle);
    }

    if(doPairing){
      // Add histograms to Pair classes (all 3D histograms)
      histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass (GeV);Pair Pt (GeV);PhiV",
                            GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
      /* histos->UserHistogram("Pair","InvMass_PairPt_Rapdity",";Inv. Mass (GeV);Pair Pt (GeV);Y_{ee}", */
      /*                       GetVector(kMee), GetVector(kPtee), BinsToVector(200, -2, 2), */
      /*                       AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY); */

      // Multiplicity
      histos->UserHistogram("Pair", "InvMass_PairPt_CentralityV0M", ";Inv. Mass (GeV);Pair Pt (GeV);CentralityV0M",
                            GetVector(kMee), GetVector(kPtee), GetVector(kCent),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityNew);
      histos->UserHistogram("Pair", "InvMass_PairPt_CentralityV0A", ";Inv. Mass (GeV);Pair Pt (GeV);CentralityV0A",
                            GetVector(kMee), GetVector(kPtee), GetVector(kCent),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityV0A);
      histos->UserHistogram("Pair", "InvMass_PairPt_CentralityV0C", ";Inv. Mass (GeV);Pair Pt (GeV);CentralityV0C",
                            GetVector(kMee), GetVector(kPtee), GetVector(kCent),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityV0C);
    }// End doPairing histograms

    // V0 feature histograms
    if(cutDefinition == "kV0_allAcc"){
      histos->UserHistogram("Pair", "CosPointingAngle", "", BinsToVector(100, 0, 1),    AliDielectronVarManager::kCosPointingAngle);
      histos->UserHistogram("Pair", "Chi2NDF",          "", BinsToVector(100, 0, 100),  AliDielectronVarManager::kChi2NDF);
      histos->UserHistogram("Pair", "LegDist",          "", BinsToVector(400, 0, 0.1),    AliDielectronVarManager::kLegDist);
      histos->UserHistogram("Pair", "R",                "", BinsToVector(1000, 0, 200), AliDielectronVarManager::kR);
      histos->UserHistogram("Pair", "PsiTrack",         "", BinsToVector(100, 0, TMath::Pi()), AliDielectronVarManager::kPsiPair);
      histos->UserHistogram("Pair", "kM",               "", BinsToVector(2000, 0, 20),  AliDielectronVarManager::kM);

      histos->UserHistogram("Pair", "ArmAlpha_armPt", "", BinsToVector(400, -2.5, 2.5), BinsToVector(500, 0, 3),
                            AliDielectronVarManager::kArmAlpha, AliDielectronVarManager::kArmPt);
    }
    die->SetHistogramManager(histos);
}

TVectorD* GetVector(Int_t var){

  switch(var){

    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(100, 0, 2*TMath::Pi());
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);

    case kSigmaEle:
      return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      return AliDielectronHelper::MakeLinBinning( 100,-10.,20.);
    case kTPCdEdx:
      return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);

    // Mass bins.
    case kMee:
      return AliDielectronHelper::MakeArbitraryBinning("0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,"
                                                        "0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,"
                                                        "0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,"
                                                        "0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,"
                                                        "0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,"
                                                        "0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,"
                                                        "0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,"
                                                        "0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,"
                                                        "0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,"
                                                        "0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,"
                                                        "1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,"
                                                        "1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,"
                                                        "2.10,2.20,2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,"
                                                        "3.01,3.02,3.03,3.04,3.05,3.06,3.07,3.08,3.09,3.10,"
                                                        "3.11,3.12,3.30,3.50,4.00,4.50,5.00");
    // Pair pt bins
    case kPtee:
      return AliDielectronHelper::MakeArbitraryBinning("0.000,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,"
                                                       "0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.550,"
                                                       "0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.050,1.100,"
                                                       "1.150,1.200,1.250,1.300,1.350,1.400,1.450,1.500,1.550,1.600,1.650,"
                                                       "1.700,1.750,1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.200,"
                                                       "2.250,2.300,2.350,2.400,2.450,2.500,2.600,2.700,2.800,2.900,3.000,"
                                                       "3.100,3.200,3.300,3.400,3.500,3.600,3.700,3.800,3.900,4.000,4.100,"
                                                       "4.200,4.300,4.400,4.500,5.000,5.500,6.000,6.500,7.000,8.000,10.00");
    // Single leg pt bins
    case kPt:
      return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,"
                                                       "0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,"
                                                       "1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,"
                                                       "2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,"
                                                       "3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,"
                                                       "4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 5.00, 5.50, 6.00, 6.50,"
                                                       "7.00, 8.00, 10.0");
    case kCent:
      return AliDielectronHelper::MakeArbitraryBinning("0, 0.5, 5.0, 10, 20, 40, 60, 80, 100");

    case kRuns:
      return AliDielectronHelper::MakeArbitraryBinning("265300, 265309, 265332, 265334, 265336, 265338, 265339, 265342,"
                                                       "265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385,"
                                                       "265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425,"
                                                       "265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525,"
                                                       "265530");

    default: std::cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << std::endl;
      break;
  }
  return 0x0;
}

TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max){

  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
}
