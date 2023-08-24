//===============================================================================================================
// addtask for J/psi - hadron correlation analysis in pp 13TeV (last updated: 29/10/2019)
//
// objetives of this add task:
//  - data:
//    - event histograms
//    - electron and associated track histograms
//    - pair histograms (same+mixed event)
//    - correlation histograms (same+mixed event)
//  - MC:
//    - event histograms
//    - electron track histograms
//    - pair histograms
//    -> Jpsi signal shape
//    -> Jpsi efficiency
//
//===============================================================================================================

// MC signals, in order of appearance in AddTask_laltenka_dst_correlations.C
enum MCFilters {
  kJpsiInclusive=0,
  kJpsiNonPrompt,
  kJpsiPrompt,
  kJpsiRadiative,
  kJpsiNonRadiative,
  kJpsiNonPromptRadiative,
  kJpsiNonPromptNonRadiative,
  kJpsiDecayElectron,
  kJpsiNonPromptDecayElectron,
  kJpsiPromptDecayElectron,
  kJpsiRadiativeDecayElectron,
  kJpsiNonRadiativeDecayElectron,
  kJpsiDecayPhoton
};

// trigger modes
enum triggerModes {
  kMB=0,
  kV0HM,
  kSPDHM,
  kEMC7,
  kEG1,
  kEG2,
  kDG1,
  kDG2,
  kEG1DG1,
  kEG2DG2,
  kNTriggerModes=10
};

// systematics modes
enum systematicsModes {
  kStandard=0,
  kElectronPIDstrict,
  kElectronPIDopen,
  kHadronSPDstrict,
  kSystematicsFull,
  kNSystematicsModes=5
};

// running modes
enum runningModes {
  kFull=0,
  kCorr,
  kQA,
  kNRunningModes=3
};

// settings flags
Bool_t gDoEfficiencyCorrection  = kTRUE;
Bool_t gJpsiEfficiencyUsed      = kFALSE;
Bool_t gHadronEfficiencyUsed    = kFALSE;
Bool_t gUseEG1DG1PidMC          = kFALSE;
Bool_t gUseEG2DG2PidMC          = kFALSE;

// bins for mixing variables
const Int_t gNVtxZLimits              = 21;
Double_t    gVtxZLimits[gNVtxZLimits] = {-10.,-9.,-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

// forward declarations
TString SetRunNumbers(TString prod="");
void    SetJpsiEfficiencyMap(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB, Int_t systematicsMode=kStandard);
void    SetHadronEfficiencyMap(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB, Int_t systematicsMode=kStandard);
void    SetEventCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB);
void    SetClusterCuts(AliReducedAnalysisJpsi2eeCorrelations* task);
void    SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB, Int_t systematicsMode=kStandard);
void    SetAssociatedTrackCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t systematicsMode=kStandard);
void    SetTrackPrefilterCuts(AliReducedAnalysisJpsi2eeCorrelations* task);
void    SetPairCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB, Int_t systematicsMode=kStandard);
void    SetPairPrefilterCuts(AliReducedAnalysisJpsi2eeCorrelations* task);
void    SetClusterTrackMatcher(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB);
void    SetMCSignalCutsLeg(AliReducedAnalysisJpsi2eeCorrelations* task);
void    SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2eeCorrelations* task);
void    SetupMixingHandlers(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice=kMB, Int_t systematicsMode=kStandard);
void    SetupHistogramManager(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t runningMode, TString runNumbers="");
void    DefineHistograms(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t runningMode, TString runNumbers="");

//_______________________________________________________________________________________________________________
AliAnalysisTask* AddTask_laltenka_jpsi2ee_correlations(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prodString="LHC16l") {

  //
  // arguments:
  //  - isAliroot = kTRUE for ESD/AOD analysis, = KFALSE for reduced trees
  //  - runMode = 1 for AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents, = 2 for AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree
  //
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): (isAliRoot, runMode) :: (%d,%d)\n", isAliRoot, runMode);

  // settings according to production
  //-----------------------------------------------------------------------------------------------------------
  TObjArray*  prodStringArr = prodString.Tokenize("_");
  TString     prod          = prodStringArr->At(0)->GetName();
  
  // MC choice
  Bool_t setMC = kFALSE;
  if (prod.Contains("LHC17f6") || prod.Contains("LHC17f9") || prod.Contains("LHC17d17") || prod.Contains("LHC17f5") ||
      prod.Contains("LHC17d3") || prod.Contains("LHC17e5") || prod.Contains("LHC18f1") || prod.Contains("LHC18d8") ||
      prod.Contains("LHC17d16") || prod.Contains("LHC17d18"))
    setMC = kTRUE; // 2016 gen. purp. MC
  if (prod.Contains("LHC17h2a") || prod.Contains("LHC17h2b") || prod.Contains("LHC17h2d") || prod.Contains("LHC17h2e") ||
      prod.Contains("LHC17h2f") || prod.Contains("LHC17h2g") || prod.Contains("LHC17h2h2") || prod.Contains("LHC17h2i2") ||
      prod.Contains("LHC17h2j") || prod.Contains("LHC17h2k"))
    setMC = kTRUE; // 2016 J/psi inj. MC
  if (prod.Contains("LHC18d3") || prod.Contains("LHC17h1") || prod.Contains("LHC18c12") || prod.Contains("LHC17k4") ||
      prod.Contains("LHC17h11") || prod.Contains("LHC18c13") || prod.Contains("LHC18a8") || prod.Contains("LHC17l5") ||
      prod.Contains("LHC18a9") || prod.Contains("LHC18a1"))
    setMC = kTRUE; // 2017 gen. purp. MC
  if (prod.Contains("LHC18b1a"))
    setMC = kTRUE; // 2017 J/psi inj. MC
  if (prod.Contains("LHC18g4") || prod.Contains("LHC18g5") || prod.Contains("LHC18g6") || prod.Contains("LHC18h2") ||
      prod.Contains("LHC18h4") || prod.Contains("LHC18j1") || prod.Contains("LHC18j4") || prod.Contains("LHC18k1") ||
      prod.Contains("LHC18k2") || prod.Contains("LHC18k3"))
    setMC = kTRUE; // 2018 gen. purp. MC
  if (prod.Contains("LHC19c6"))
    setMC = kTRUE; // 2018 J/psi inj. MC
  
  // trigger choice
  Int_t   triggerChoice = kMB;
  TString triggerString = "";
  if (prodStringArr->GetEntries()>1) {
    triggerString = prodStringArr->At(1)->GetName();
    if (!triggerString.CompareTo("MB"))     triggerChoice = kMB;          // AliReducedVarManager::kINT7
    if (!triggerString.CompareTo("V0HM"))   triggerChoice = kV0HM;        // AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD, V0  (trigger class)
    if (!triggerString.CompareTo("SPDHM"))  triggerChoice = kSPDHM;       // AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD, SPD (trigger class)
    if (!triggerString.CompareTo("EMC7"))   triggerChoice = kEMC7;        // AliReducedVarManager::kEMC7
    if (!triggerString.CompareTo("EG1"))    triggerChoice = kEG1;         // AliReducedVarManager::kEMCEGA, EG1
    if (!triggerString.CompareTo("EG2"))    triggerChoice = kEG2;         // AliReducedVarManager::kEMCEGA, EG2
    if (!triggerString.CompareTo("DG1"))    triggerChoice = kDG1;         // AliReducedVarManager::kEMCEGA, DG1
    if (!triggerString.CompareTo("DG2"))    triggerChoice = kDG2;         // AliReducedVarManager::kEMCEGA, DG2
    if (!triggerString.CompareTo("EG1DG1")) triggerChoice = kEG1DG1;      // AliReducedVarManager::kEMCEGA, EG1 | DG1
    if (!triggerString.CompareTo("EG2DG2")) triggerChoice = kEG2DG2;      // AliReducedVarManager::kEMCEGA, EG2 | DG2
  }
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): trigger choice = %d (%s)\n", triggerChoice, triggerString.Data());
  
  // set trigger to MB for MC
  if (setMC && triggerChoice!=kMB) {
    printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): trigger choice (%s) not MB but must be MB for MC, setting trigger to MB... \n", triggerString.Data());
    if (triggerChoice==kEG1 || triggerChoice==kDG1 || triggerChoice==kEG1DG1) gUseEG1DG1PidMC = kTRUE;
    if (triggerChoice==kEG2 || triggerChoice==kDG2 || triggerChoice==kEG2DG2) gUseEG2DG2PidMC = kTRUE;
    triggerChoice = kMB;
  }
  
  // systematics choice
  Int_t   systematicsMode   = kStandard;
  TString systematicsString = "";
  if (prodStringArr->GetEntries()>2) {
    systematicsString = prodStringArr->At(2)->GetName();
    if (!systematicsString.CompareTo("standard"))           systematicsString = "";
    if (!systematicsString.CompareTo("systematics"))        systematicsMode   = kSystematicsFull;
    if (!systematicsString.CompareTo("electronPIDstrict"))  systematicsMode   = kElectronPIDstrict;
    if (!systematicsString.CompareTo("electronPIDopen"))    systematicsMode   = kElectronPIDopen;
    if (!systematicsString.CompareTo("hadronSPDstrict"))    systematicsMode   = kHadronSPDstrict;
  }
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): systematics mode = %d (%s)\n", systematicsMode, systematicsString.Data());
  
  // mode choice
  Int_t   runningMode   = kFull;
  TString runningString = "";
  if (prodStringArr->GetEntries()>3) {
    runningString = prodStringArr->At(3)->GetName();
    if (!runningString.CompareTo("corr")) runningMode = kCorr;
    if (!runningString.CompareTo("QA"))   runningMode = kQA;
  }
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): running mode = %d (%s)\n", runningMode, runningString.Data());

  // initialize analysis task
  //-----------------------------------------------------------------------------------------------------------
  AliReducedAnalysisJpsi2eeCorrelations* jpsi2eeCorrAnalysis = new AliReducedAnalysisJpsi2eeCorrelations("Jpsi2eeCorrelationAnalysis","Jpsi->ee - hadron correlation analysis");
  //-----------------------------------------------------------------------------------------------------------
  jpsi2eeCorrAnalysis->SetLoopOverTracks(kTRUE);
  jpsi2eeCorrAnalysis->SetRunPairing(kTRUE);
  jpsi2eeCorrAnalysis->SetRunPrefilter(kTRUE);
  jpsi2eeCorrAnalysis->SetFillCaloClusterHistograms(kTRUE);
  //-----------------------------------------------------------------------------------------------------------
  if (setMC) {  // only J/psi signal extraction in MC
    jpsi2eeCorrAnalysis->SetStoreJpsiCandidates(kFALSE);
    jpsi2eeCorrAnalysis->SetRunEventMixing(kFALSE);
    jpsi2eeCorrAnalysis->SetRunLikeSignPairing(kFALSE);
    jpsi2eeCorrAnalysis->SetRunCorrelation(kFALSE);
    jpsi2eeCorrAnalysis->SetRunCorrelationMixing(kFALSE);
    jpsi2eeCorrAnalysis->SetUseLikeSignPairs(kFALSE);
  } else {      // full correlation analysis
    jpsi2eeCorrAnalysis->SetStoreJpsiCandidates(kTRUE);
    jpsi2eeCorrAnalysis->SetRunEventMixing(kTRUE);
    jpsi2eeCorrAnalysis->SetRunLikeSignPairing(kTRUE);
    jpsi2eeCorrAnalysis->SetRunCorrelation(kTRUE);
    jpsi2eeCorrAnalysis->SetRunCorrelationMixing(kTRUE);
    jpsi2eeCorrAnalysis->SetUseLikeSignPairs(kTRUE);
  }
  //-----------------------------------------------------------------------------------------------------------
  jpsi2eeCorrAnalysis->SetRunOverMC(setMC);

  // set run numbers for runwise histograms
  //-----------------------------------------------------------------------------------------------------------
  TString runNumbers = SetRunNumbers(prod);
  
  // set hadron efficiency maps
  //-----------------------------------------------------------------------------------------------------------
  if (!setMC && gDoEfficiencyCorrection) {
    SetJpsiEfficiencyMap(   jpsi2eeCorrAnalysis, triggerChoice, systematicsMode);
    SetHadronEfficiencyMap( jpsi2eeCorrAnalysis, triggerChoice, systematicsMode);
  }
    
  // set analysis and prefilter cuts
  //-----------------------------------------------------------------------------------------------------------
  SetEventCuts(jpsi2eeCorrAnalysis, triggerChoice);
  SetClusterCuts(jpsi2eeCorrAnalysis);
  SetJpsiElectronTrackCuts(jpsi2eeCorrAnalysis, triggerChoice, systematicsMode);
  SetTrackPrefilterCuts(jpsi2eeCorrAnalysis);
  SetPairCuts(jpsi2eeCorrAnalysis, triggerChoice, systematicsMode);
  SetPairPrefilterCuts(jpsi2eeCorrAnalysis);
  SetClusterTrackMatcher(jpsi2eeCorrAnalysis, triggerChoice);
  SetAssociatedTrackCuts(jpsi2eeCorrAnalysis, systematicsMode);

  // set MC signal selection
  //-----------------------------------------------------------------------------------------------------------
  Bool_t isMC = jpsi2eeCorrAnalysis->GetRunOverMC();
  if (isMC) {
    SetMCSignalCutsLeg(jpsi2eeCorrAnalysis);
    SetMCSignalCutsJPsi(jpsi2eeCorrAnalysis);
  }

  // set MC pT weights
  //-----------------------------------------------------------------------------------------------------------
  if (isMC) {
    TFile*  fW      = TFile::Open("/Users/lal035/Documents/PWGDQ/data/calibration/MC/ratio_13TeV_two.root");
    TH1F*   weights = (TH1F*)fW->Get("weightsNew");
    jpsi2eeCorrAnalysis->SetMCJpsiPtWeights(weights);
  }

  // initialize analysis task
  //-----------------------------------------------------------------------------------------------------------
  jpsi2eeCorrAnalysis->Init();
  
  // define histograms histograms
  //-----------------------------------------------------------------------------------------------------------
  SetupHistogramManager(jpsi2eeCorrAnalysis, runningMode, runNumbers);

  // mixing handler
  //-----------------------------------------------------------------------------------------------------------
  SetupMixingHandlers(jpsi2eeCorrAnalysis, triggerChoice, systematicsMode);

  // initialize wrapper AliAnalysisTask
  // (in order to run AliReducedAnalysisJpsi2eeCorrelations in an aliroot analysis train )
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeCorrAnalysis);
  
  // intercept isAliRoot=kFALSE (nothing to be done yet)
  //-----------------------------------------------------------------------------------------------------------
  if (!isAliRoot) return 0;

  // get analysis manager
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_laltenka_jpsi2ee_correlations", "No analysis manager found.");
    return 0;
  }

  // get data container
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisDataContainer* cReducedEvent = NULL;
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
    printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): use on the fly events\n");
    cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
    if (!cReducedEvent) {
      printf("ERROR: In AddTask_laltenka_jpsi2ee_correlations(), ReducedEvent exchange container could not be found!\n");
      return 0x0;
    }
  }

  // add task to analysis manager
  //-----------------------------------------------------------------------------------------------------------
  mgr->AddTask(task);

  // connect input data containers
  //-----------------------------------------------------------------------------------------------------------
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)              mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  else if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)  mgr->ConnectInput(task, 0, cReducedEvent);
  else {
    printf("ERROR: In AddTask_laltenka_jpsi2ee_correlations(), runMode %d not defined!\n", runMode);
    return 0;
  }

  // connect output data containers
  //-----------------------------------------------------------------------------------------------------------
  TString                               outputName = "dstCorrelationsAnalysisHistograms.root";
  if (triggerString.CompareTo(""))      outputName = Form("dstCorrelationsAnalysisHistograms_%s.root", triggerString.Data());
  if (systematicsString.CompareTo(""))  outputName.ReplaceAll(".root", Form("_%s.root", systematicsString.Data()));
  if (runningString.CompareTo(""))      outputName.ReplaceAll(".root", Form("_%s.root", runningString.Data()));
  if (isMC)                             outputName.ReplaceAll(".root", "_MC.root");
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): output container name = %s\n", outputName.Data());
  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  mgr->ConnectOutput(task, 1, cOutputHist);
  
  // done
  //-----------------------------------------------------------------------------------------------------------
  return task;
}

//_______________________________________________________________________________________________________________
TString SetRunNumbers(TString prod /*=""*/) {
  
  TString runNumbers = "";
  if (!prod.CompareTo("")) {
    cout << "WARNING in AddTask_laltenka_unbiasedEvents::SetRunNumbers(), production " << prod.Data() << " empty, run numbers not defined!" << endl;
    AliReducedVarManager::SetRunNumbers(runNumbers);
    return runNumbers;
  }
  
  // 2016 data / MC
  //___________________________________________________________________________________________________________
  else if (prod.Contains("LHC16d") || prod.Contains("LHC17f6")    || prod.Contains("LHC17h2a"))   runNumbers = "252330;252326;252325;252322;252319;252317;252310;252271;252248;252235";
  else if (prod.Contains("LHC16e") || prod.Contains("LHC17f9")    || prod.Contains("LHC17h2b"))   runNumbers = "253591;253589;253563;253530;253529;253517;253488;253482;253481;253478;253437";
  else if (prod.Contains("LHC16g") || prod.Contains("LHC17d17")   || prod.Contains("LHC17h2d"))   runNumbers = "254332;254331;254330;254304;254303;254302;254293;254205;254204;254199;254193;254178;254175;254174;254149;254148;254147;254128";
  else if (prod.Contains("LHC16h") || prod.Contains("LHC17f5")    || prod.Contains("LHC17h2e"))   runNumbers = "255467;255466;255465;255463;255447;255442;255440;255421;255420;255419;255418;255415;255407;255402;255398;255352;255351;255350;255283;255280;255276;255275;255256;255255;255253;255252;255251;255249;255248;255247;255242;255240;255182;255181;255180;255177;255176;255174;255173;255171;255167;255162;255159;255154;255111;255091;255086;255085;255082;255079;254984;254983;254654;254653;254652;254651;254649;254648;254646;254644;254640;254632;254630;254629;254621;254608;254606;254604";
  else if (prod.Contains("LHC16i") || prod.Contains("LHC17d3")    || prod.Contains("LHC17h2f"))   runNumbers = "255618;255617;255616;255615;255614;255591;255583;255582;255577;255543;255542;255541;255540;255539";
  else if (prod.Contains("LHC16j") || prod.Contains("LHC17e5")    || prod.Contains("LHC17h2g"))   runNumbers = "256420;256418;256417;256415;256373;256372;256371;256368;256366;256365;256364;256363;256362;256361;256357;256356;256311;256309;256307;256302;256299;256298;256297;256295;256292;256290;256289;256287;256284;256283;256282;256281;256231;256228;256227;256223;256222;256219;256215;256213;256212;256210;256207;256204";
  else if (prod.Contains("LHC16k") || prod.Contains("LHC18f1")    || prod.Contains("LHC17h2h2"))  runNumbers = "258537;258499;258477;258456;258454;258452;258426;258393;258391;258387;258359;258336;258332;258307;258306;258303;258302;258301;258299;258278;258274;258273;258271;258270;258258;258257;258256;258204;258203;258202;258198;258197;258178;258117;258114;258113;258109;258108;258107;258063;258062;258060;258059;258053;258049;258048;258045;258042;258041;258039;258019;258017;258014;258012;258008;258003;257992;257989;257986;257979;257963;257960;257958;257957;257939;257937;257936;257893;257855;257853;257851;257850;257804;257803;257800;257799;257798;257797;257773;257765;257757;257754;257737;257735;257734;257733;257727;257725;257724;257697;257694;257692;257691;257689;257688;257687;257685;257684;257682;257644;257642;257636;257635;257632;257630;257606;257605;257604;257601;257595;257594;257592;257590;257588;257587;257566;257562;257561;257560;257541;257540;257539;257537;257531;257530;257492;257491;257490;257488;257487;257474;257468;257457;257433;257364;257358;257330;257322;257320;257318;257260;257224;257209;257206;257204;257145;257144;257142;257141;257140;257139;257138;257137;257136;257100;257095;257092;257086;257084;257083;257082;257080;257077;257028;257026;257021;257012;257011;256944;256942;256941";
  else if (prod.Contains("LHC16l") || prod.Contains("LHC18d8")    || prod.Contains("LHC17h2i2"))  runNumbers = "259888;259868;259867;259866;259860;259842;259841;259822;259789;259788;259781;259756;259752;259751;259750;259748;259747;259477;259473;259396;259395;259394;259389;259388;259382;259378;259342;259341;259340;259339;259336;259334;259307;259305;259303;259302;259274;259273;259272;259271;259270;259269;259264;259263;259261;259257;259204;259164;259162;259118;259117;259099;259096;259091;259090;259088;258964;258962";
  else if (prod.Contains("LHC16o") || prod.Contains("LHC17d16")   || prod.Contains("LHC17h2j"))   runNumbers = "264035;264033;263985;263984;263981;263978;263977;263923;263920;263917;263916;263905;263866;263863;263861;263810;263803;263793;263792;263790;263787;263786;263785;263784;263744;263743;263741;263739;263738;263737;263691;263690;263689;263682;263663;263662;263657;263654;263653;263652;263647;263529;263497;263496;263490;263487;263332;263331;262858;262855;262853;262849;262847;262844;262842;262841;262778;262777;262776;262768;262760;262727;262725;262723;262719;262717;262713;262708;262706;262705;262428;262426;262425;262424";
  else if (prod.Contains("LHC16p") || prod.Contains("LHC17d18")   || prod.Contains("LHC17h2k"))   runNumbers = "264347;264346;264345;264341;264336;264312;264306;264305;264281;264279;264277;264273;264267;264266;264265;264264;264262;264261;264260;264259;264238;264235;264233;264232;264198;264197;264194;264190;264188;264168;264164;264139;264138;264137;264129;264110;264109;264086;264085;264082;264078;264076";
  
  // 2017 data / MC
  //___________________________________________________________________________________________________________
  else if (prod.Contains("LHC17e") || prod.Contains("LHC17h1")    || prod.Contains("LHC18b1a-e")) runNumbers = "270830;270828;270827;270824;270822";
  else if (prod.Contains("LHC17f") || prod.Contains("LHC18d3")    || prod.Contains("LHC18b1a-f")) runNumbers = "270865;270861;270856;270855;270854";
  else if (prod.Contains("LHC17h") || prod.Contains("LHC18c12")   || prod.Contains("LHC18b1a-h")) runNumbers = "273103;273100;273099;273077;273010;273009;272985;272983;272976;272949;272947;272939;272935;272934;272933;272932;272905;272903;272880;272873;272871;272870;272836;272834;272833;272829;272828;272784;272783;272782;272764;272763;272760;272749;272747;272746;272712;272692;272691;272690;272620;272610;272608;272607;272585;272577;272575;272574;272521;272468;272466;272463;272462;272461;272413;272411;272400;272399;272395;272394;272389;272388;272360;272359;272340;272335;272194;272156;272155;272154;272153;272152;272151;272123;272101;272100;272076;272042;272041;272040;272039;272038;272036;272020;272018;271886;271881;271880;271874;271873;271871;271870;271868";
  else if (prod.Contains("LHC17i") || prod.Contains("LHC17k4")    || prod.Contains("LHC18b1a-i")) runNumbers = "274442;274390;274389;274388;274387;274386;274385;274364;274363;274360;274352;274351;274329;274283;274281;274280;274278;274276;274271;274270;274269;274268;274266;274264;274263;274259;274258;274232;274212;274174;274148;274147;274125;274094;274092;274064;274058;273986;273985;273946;273943;273942;273918;273889;273887;273886;273885;273825;273824;273719;273711;273709;273695;273690;273689;273687;273654;273653;273593;273592;273591";
  else if (prod.Contains("LHC17j") || prod.Contains("LHC17h11")   || prod.Contains("LHC18b1a-j")) runNumbers = "274671;274669;274667;274657;274653;274601;274596;274595;274594;274593";
  else if (prod.Contains("LHC17k") || prod.Contains("LHC18c13")   || prod.Contains("LHC18b1a-k")) runNumbers = "276508;276507;276506;276462;276439;276438;276437;276435;276351;276348;276312;276307;276302;276297;276294;276292;276291;276290;276259;276257;276230;276205;276178;276177;276170;276169;276166;276145;276140;276135;276104;276102;276099;276098;276097;276045;276041;276040;276020;276019;276017;276013;276012;275925;275924;275847;275664;275661;275657;275650;275648;275647;275624;275623;275622;275621;275617;275612;275559;275558;275515;275472;275471;275467;275459;275457;275456;275453;275452;275448;275443;275406;275404;275401;275395;275394;275372;275369;275361;275360;275333;275332;275328;275326;275324;275322;275314;275283;275247;275246;275245;275239;275188;275184;275180;275177;275174;275173;275151;275150;275149;275076;275075;275073;275068;275067;274979;274978;274886;274882;274878;274877;274822;274821;274817;274815;274811;274807;274806;274803;274802;274801;274708;274690";
  else if (prod.Contains("LHC17l") || prod.Contains("LHC18a8")    || prod.Contains("LHC18b1a-l")) runNumbers = "278216;278215;278191;278189;278167;278166;278165;278164;278163;278158;278127;278126;278123;278122;278121;277996;277991;277989;277987;277952;277930;277907;277904;277903;277901;277900;277899;277898;277897;277876;277870;277848;277847;277845;277842;277841;277836;277834;277805;277802;277801;277800;277799;277795;277794;277749;277747;277746;277745;277725;277723;277722;277721;277577;277576;277575;277574;277537;277536;277534;277531;277530;277479;277478;277477;277476;277473;277472;277418;277417;277416;277389;277386;277385;277384;277383;277360;277314;277312;277310;277293;277262;277257;277256;277197;277196;277194;277193;277189;277188;277184;277183;277182;277181;277180;277155;277121;277117;277091;277087;277082;277079;277076;277073;277037;277017;277016;277015;276972;276971;276970;276969;276967;276920;276917;276916;276762;276675;276674;276672;276671;276670;276644;276608;276557;276556;276553;276552;276551";
  else if (prod.Contains("LHC17m") || prod.Contains("LHC17l5")    || prod.Contains("LHC18b1a-m")) runNumbers = "280140;280135;280134;280131;280126;280118;280114;280111;280108;280107;280066;280052;280051;279880;279879;279855;279854;279853;279830;279827;279826;279773;279749;279747;279719;279718;279715;279689;279688;279687;279684;279683;279682;279679;279677;279676;279642;279641;279632;279630;279559;279550;279491;279488;279487;279483;279441;279439;279435;279410;279391;279355;279354;279349;279348;279344;279342;279312;279310;279309;279274;279273;279270;279268;279267;279265;279264;279242;279238;279235;279234;279232;279208;279207;279201;279199;279157;279155;279130;279123;279122;279118;279117;279107;279106;279075;279074;279073;279069;279068;279044;279043;279041;279036;279035;279008;279007;279005;279000;278999;278964;278963;278960;278959;278941;278939;278936;278915;278914";
  else if (prod.Contains("LHC17o") || prod.Contains("LHC18a9")    || prod.Contains("LHC18b1a-o")) runNumbers = "281961;281956;281953;281940;281939;281932;281931;281928;281920;281918;281916;281915;281895;281894;281893;281892;281633;281592;281583;281574;281569;281568;281563;281562;281557;281511;281509;281477;281475;281450;281449;281446;281444;281443;281441;281415;281321;281301;281277;281275;281273;281271;281244;281243;281242;281241;281240;281213;281212;281191;281190;281189;281181;281180;281179;281081;281080;281062;281061;281060;281036;281035;281033;281032;280999;280998;280997;280996;280994;280990;280947;280943;280940;280936;280897;280890;280881;280880;280856;280854;280849;280848;280847;280845;280844;280842;280793;280792;280787;280786;280768;280767;280766;280765;280764;280763;280762;280761;280757;280756;280755;280754;280753;280729;280706;280705;280681;280679;280671;280647;280645;280639;280637;280636;280634;280613;280583;280581;280576;280575;280574;280551;280550;280547;280546;280519;280518;280499;280490;280448;280447;280446;280445;280443;280419;280415;280413;280412;280406;280405;280403;280375;280374;280351;280350;280349;280348;280312;280310;280290;280286;280285;280284;280283;280282";
  else if (prod.Contains("LHC17r") || prod.Contains("LHC18a1")    || prod.Contains("LHC18b1a-r")) runNumbers = "282704;282703;282702;282700;282677;282676;282673;282671;282670;282668;282667;282666;282653;282651;282629;282622;282620;282618;282609;282608;282607;282606;282580;282579;282575;282573;282546;282545;282544;282528";

  // 2018 data
  //___________________________________________________________________________________________________________
  else if (prod.Contains("LHC18b") || prod.Contains("LHC18g4")    || prod.Contains("LHC19c6-b"))  runNumbers = "285396;285365;285364;285347;285328;285327;285224;285222;285203;285202;285200;285165;285127;285125;285108;285106;285066;285065;285064;285015;285014;285013;285012;285011;285009";
  else if (prod.Contains("LHC18d") || prod.Contains("LHC18g5")    || prod.Contains("LHC19c6-d"))  runNumbers = "286350;286349;286348;286345;286341;286340;286337;286336;286314;286313;286312;286311;286310;286309;286308;286289;286288;286287;286284;286282;286263;286261;286258;286257;286255;286254;286231;286230;286229;286203;286202;286201;286199;286198;286159;286130;286129;286127;286124;286064;286025;286014;285980;285979;285978";
  else if (prod.Contains("LHC18e") || prod.Contains("LHC18g6")    || prod.Contains("LHC19c6-e"))  runNumbers = "286937;286936;286933;286932;286931;286930;286911;286910;286908;286907;286877;286876;286874;286852;286850;286848;286846;286810;286809;286805;286801;286799;286731;286695;286661;286653;286633;286594;286592;286591;286569;286568;286567;286566;286511;286509;286508;286502;286501;286482;286455;286454;286428;286427;286426;286380";
  else if (prod.Contains("LHC18f") || prod.Contains("LHC18h2")    || prod.Contains("LHC19c6-f"))  runNumbers = "287658;287657;287656;287654;287578;287576;287575;287573;287524;287521;287520;287518;287517;287516;287513;287486;287484;287481;287480;287451;287413;287389;287388;287387;287385;287381;287380;287360;287356;287355;287353;287349;287347;287346;287344;287343;287325;287324;287323;287283;287254;287251;287250;287249;287248;287209;287208;287204;287203;287202;287201;287185;287155;287137;287077;287072;287071;287066;287064;287063;287021;287000";
  else if (prod.Contains("LHC18g") || prod.Contains("LHC18h4-g")  || prod.Contains("LHC19c6-g"))  runNumbers = "288619;288640;288642;288644;288650;288687;288689;288690;288743;288748;288750";
  else if (prod.Contains("LHC18h") || prod.Contains("LHC18h4-h")  || prod.Contains("LHC19c6-h"))  runNumbers = "288804;288806";
  else if (prod.Contains("LHC18i") || prod.Contains("LHC18h4-i")  || prod.Contains("LHC19c6-i"))  runNumbers = "288861;288862;288863;288864;288868;288902;288903;288908;288909";
  else if (prod.Contains("LHC18j") || prod.Contains("LHC18h4-j")  || prod.Contains("LHC19c6-j"))  runNumbers = "288943";
  else if (prod.Contains("LHC18k") || prod.Contains("LHC18h4-k")  || prod.Contains("LHC19c6-k"))  runNumbers = "289165;289166;289167;289169;289172;289175;289176;289177;289198;289199;289200;289201";
  else if (prod.Contains("LHC18l") || prod.Contains("LHC18j1")    || prod.Contains("LHC19c6-l"))  runNumbers = "289971;289966;289965;289943;289941;289940;289935;289931;289928;289884;289880;289879;289857;289856;289855;289854;289852;289849;289830;289818;289817;289816;289815;289814;289811;289808;289775;289757;289732;289731;289729;289724;289723;289721;289547;289521;289494;289493;289468;289466;289465;289463;289462;289444;289426;289374;289373;289370;289369;289368;289367;289366;289365;289356;289355;289354;289353;289309;289308;289306;289303;289300;289281;289280;289278;289277;289276;289275;289254;289253;289249;289247;289243;289242;289241;289240";
  else if (prod.Contains("LHC18m") || prod.Contains("LHC18j4")    || prod.Contains("LHC19c6-m"))  runNumbers = "292839;292836;292834;292832;292831;292811;292810;292809;292804;292803;292758;292754;292752;292750;292748;292747;292744;292739;292737;292704;292701;292698;292696;292695;292693;292586;292584;292563;292560;292559;292557;292554;292553;292526;292524;292523;292521;292500;292497;292496;292495;292461;292460;292457;292456;292434;292432;292430;292429;292428;292406;292405;292398;292397;292298;292274;292273;292270;292269;292265;292242;292241;292240;292218;292192;292168;292167;292166;292164;292163;292162;292161;292160;292140;292115;292114;292109;292108;292107;292106;292081;292080;292077;292075;292067;292062;292061;292060;292040;292012;291982;291977;291976;291953;291948;291946;291945;291944;291943;291942;291803;291796;291795;291769;291768;291766;291762;291760;291756;291755;291729;291706;291702;291698;291697;291694;291692;291690;291665;291661;291657;291626;291624;291622;291618;291615;291614;291590;291485;291484;291482;291481;291457;291456;291453;291451;291447;291446;291424;291420;291419;291417;291416;291402;291400;291399;291397;291377;291375;291363;291362;291361;291360;291286;291285;291284;291283;291282;291266;291265;291263;291262;291257;291240;291209;291188;291143;291116;291111;291110;291101;291100;291093;291069;291066;291065;291041;291037;291035;291006;291005;291004;291003;291002;290980;290979;290976;290975;290974;290948;290944;290943;290941;290935;290932;290895;290894;290892;290888;290887;290886;290862;290860;290853;290848;290846;290843;290841;290790;290787;290776;290774;290769;290766;290764;290721;290699;290696;290692;290689;290687;290665;290660;290658;290645;290632;290627;290615;290614;290613;290612;290590;290588;290553;290550;290549;290544;290540;290539;290538;290501;290500;290499;290469;290467;290459;290458;290456;290428;290427;290426;290425;290423;290412;290411;290404;290401;290399;290376;290375;290374;290350;290327;290324;290323";
  else if (prod.Contains("LHC18n") || prod.Contains("LHC18k1")    || prod.Contains("LHC19c6-n"))  runNumbers = "293357;293359";
  else if (prod.Contains("LHC18o") || prod.Contains("LHC18k2")    || prod.Contains("LHC19c6-o"))  runNumbers = "293898;293896;293893;293891;293886;293856;293831;293830;293829;293809;293807;293806;293805;293802;293799;293776;293774;293773;293770;293741;293740;293698;293696;293695;293692;293691;293690;293689;293686;293588;293587;293583;293582;293579;293578;293573;293571;293570;293475";
  else if (prod.Contains("LHC18p") || prod.Contains("LHC18k3")    || prod.Contains("LHC19c6-p"))  runNumbers = "294925;294916;294884;294883;294880;294877;294875;294852;294818;294817;294816;294815;294813;294809;294805;294775;294774;294772;294769;294749;294747;294746;294745;294744;294743;294742;294741;294722;294721;294718;294716;294715;294710;294703;294653;294636;294634;294633;294632;294593;294591;294590;294588;294587;294586;294563;294562;294558;294556;294553;294531;294530;294529;294527;294526;294525;294524;294310;294308;294307;294242;294241;294212;294210;294208;294205;294201;294200;294199;294156;294155;294154;294152;294131;294013;294012;294011;294010;294009";
  
  // run numbers undefined
  else {
    cout << "WARNING in AddTask_laltenka_assocHadronEfficiency::SetRunNumbers(), production " << prod.Data() << " not known, run numbers not defined!" << endl;
  }
  
  AliReducedVarManager::SetRunNumbers(runNumbers);
  return runNumbers;
}

//_______________________________________________________________________________________________________________
void SetJpsiEfficiencyMap(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/, Int_t systematicsMode /*=kStandard*/) {

  // get correct efficiency file
  TString                           triggerName = "";
  if      (triggerChoice==kMB)      triggerName = "MB";
  else if (triggerChoice==kV0HM)    triggerName = "HM";
  else if (triggerChoice==kEG1DG1)  triggerName = "EG1DG1";
  else if (triggerChoice==kEG2DG2)  triggerName = "EG2DG2";
  else {
    cout << "WARNING in AddTask_laltenka_jpsi2ee_correlations::SetJpsiEfficiencyMap(), trigger choice " << triggerChoice << " not known, using MB efficiencies!" << endl;
    triggerName = "MB";
  }

  TString                                       electronName = "";
  if      (systematicsMode==kStandard)          electronName = "standardElectron";
  else if (systematicsMode==kElectronPIDstrict) electronName = "PIDstrictElectron";
  else if (systematicsMode==kElectronPIDopen)   electronName = "PIDopenElectron";
  else if (systematicsMode==kHadronSPDstrict)   electronName = "standardElectron";
  else {
    cout << "WARNING in AddTask_laltenka_jpsi2ee_correlations::SetJpsiEfficiencyMap(), systematics choice " << systematicsMode << " not known, using standard efficiencies!" << endl;
    electronName = "standardElectron";
  }
  
  TString effFileName = Form("/Users/lal035/Documents/PWGDQ/analysis_MC/jpsi/efficiencies/efficiency_jpsi_%s_%s_average.root", electronName.Data(), triggerName.Data());
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): J/psi efficiencies loaded from %s\n", effFileName.Data());
  TFile* effFile = TFile::Open(effFileName);
  effFile->cd();
  
  // read efficiency from file
  AliReducedVarManager::SetPairEfficiencyMap((TH1F*)gDirectory->Get("jpsiEff_mcTruthJpsi_pT_rebin"), AliReducedVarManager::kPt);
  AliReducedVarManager::SetPairEfficiencyMapDependeciesCorrelation(AliReducedVarManager::kTriggerPt);

  // close file and set global flag
  effFile->Close();
  gJpsiEfficiencyUsed = kTRUE;
}

//_______________________________________________________________________________________________________________
void SetHadronEfficiencyMap(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/, Int_t systematicsMode /*=kStandard*/) {
  
  // get correct efficiency file
  TString                           triggerName = "";
  if      (triggerChoice==kMB)      triggerName = "MB";
  else if (triggerChoice==kV0HM)    triggerName = "HM";
  else if (triggerChoice==kEG1DG1)  triggerName = "EG1DG1";
  else if (triggerChoice==kEG2DG2)  triggerName = "EG2DG2";
  else {
    cout << "WARNING in AddTask_laltenka_jpsi2ee_correlations::SetHadronEfficiencyMap(), trigger choice " << triggerChoice << " not known, using MB efficiencies!" << endl;
    triggerName = "MB";
  }

  TString                                       hadronName = "";
  if      (systematicsMode==kStandard)          hadronName = "standardHadron";
  else if (systematicsMode==kElectronPIDstrict) hadronName = "standardHadron";
  else if (systematicsMode==kElectronPIDopen)   hadronName = "standardHadron";
  else if (systematicsMode==kHadronSPDstrict)   hadronName = "SPDstrictHadron";
  else {
    cout << "WARNING in AddTask_laltenka_jpsi2ee_correlations::SetHadronEfficiencyMap(), systematics choice " << systematicsMode << " not known, using standard efficiencies!" << endl;
    hadronName = "standardHadron";
  }
  
  TString effFileName = Form("/Users/lal035/Documents/PWGDQ/analysis_MC/hadron/efficiencies/efficiency_hadron_%s_%s_average.root", hadronName.Data(), triggerName.Data());
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): Hadron efficiencies loaded from %s\n", effFileName.Data());
  TFile* effFile = TFile::Open(effFileName);
  effFile->cd();

  // read efficiency from file
  AliReducedVarManager::SetAssociatedHadronEfficiencyMap((TH1F*)gDirectory->Get("hadronEff_proxy_pT_rebin_purity"), AliReducedVarManager::kAssociatedPt);

  // close file and set global flag
  effFile->Close();
  gHadronEfficiencyUsed = kTRUE;
}

//_______________________________________________________________________________________________________________
void SetEventCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/) {
  
  // default selection (vertex z and PS)
  AliReducedEventCut* defaultCut = new AliReducedEventCut("DefaultCut", "Default selection");
  defaultCut->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  if (!task->GetRunOverMC()) defaultCut->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);

  // trigger selection
  AliReducedEventCut* triggerCut = new AliReducedEventCut("TriggerCut", "Trigger selection");
  if (      triggerChoice==kMB)     triggerCut->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  else if ( triggerChoice==kV0HM)   triggerCut->AddEventTriggerFilter(AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD);
  else if ( triggerChoice==kSPDHM)  triggerCut->AddEventTriggerFilter(AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD);
  else if ( triggerChoice==kEMC7)   triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMC7);
  else if ( triggerChoice==kEG1)    triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else if ( triggerChoice==kEG2)    triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else if ( triggerChoice==kDG1)    triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else if ( triggerChoice==kDG2)    triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else if ( triggerChoice==kEG1DG1) triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else if ( triggerChoice==kEG2DG2) triggerCut->AddEventTriggerFilter(AliReducedVarManager::kEMCEGA);
  else {
    printf("WARNING: In AddTask_laltenka_unbiasedEvents(), trigger choice %d not defined, using INT7!\n", triggerChoice);
    triggerCut->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  }
  if (triggerChoice==kV0HM) triggerCut->AddEventTriggerClassFilter("CVHMV0M");
  if (triggerChoice==kSPDHM) {
    triggerCut->AddEventTriggerClassFilter("CVHMSH2");
    triggerCut->AddEventTriggerClassFilter("CSHM7");
    triggerCut->EventTriggerClassFilterUseOr();
  }
  if (triggerChoice==kEG1) triggerCut->AddEventTriggerClassFilter("EG1");
  if (triggerChoice==kEG2) triggerCut->AddEventTriggerClassFilter("EG2");
  if (triggerChoice==kDG1) triggerCut->AddEventTriggerClassFilter("DG1");
  if (triggerChoice==kDG2) triggerCut->AddEventTriggerClassFilter("DG2");
  if (triggerChoice==kEG1DG1) {
    triggerCut->AddEventTriggerClassFilter("EG1");
    triggerCut->AddEventTriggerClassFilter("DG1");
    triggerCut->EventTriggerClassFilterUseOr(); // AliReducedVarManager::kEMCEGA, EG1 || DG1
  }
  if (triggerChoice==kEG2DG2) {
    triggerCut->AddEventTriggerClassFilter("EG2");
    triggerCut->AddEventTriggerClassFilter("DG2");
    triggerCut->EventTriggerClassFilterUseOr(); // AliReducedVarManager::kEMCEGA, EG2 || DG2
  }

  // SPD pileup cut (https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup)
  AliReducedEventCut* spdPileupCut = new AliReducedEventCut("SpdPileupCut", "SPD pileup rejection");
  if (triggerChoice==kV0HM || triggerChoice==kSPDHM)  spdPileupCut->AddEventTagFilterBit(11,  kTRUE); // event->IsPileupFromSPD(5,0.8,3.,2.,5.)
  else                                                spdPileupCut->AddEventTagFilterBit(9,   kTRUE); // event->IsPileupFromSPD(3,0.8,3.,2.,5.)

  // MV pileup cut (https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup)
  AliReducedEventCut* mvPileupCut = new AliReducedEventCut("MvPileupCut", "MV pileup rejection");
  mvPileupCut->AddEventTagFilterBit(1, kTRUE); // reject multi-vertexer pileup events, min. weight. dist. 15

  // combine all cuts
  AliReducedCompositeCut* eventCut = new AliReducedCompositeCut("eventCut", "", kTRUE);
  eventCut->AddCut(defaultCut);
  eventCut->AddCut(triggerCut);
  eventCut->AddCut(spdPileupCut);
  eventCut->AddCut(mvPileupCut);
  task->AddEventCut(eventCut);
  
  // set additional MB event cut required for event mixing with EMCal triggered events
  if (triggerChoice>=kEMC7) {
    AliReducedEventCut* triggerCutMB = new AliReducedEventCut("TriggerCutMB", "Trigger selection MB");
    triggerCutMB->AddEventTriggerFilter(AliReducedVarManager::kINT7);
    
    AliReducedEventCut* spdPileupCutMB = new AliReducedEventCut("spdPileupCutMB", "SPD pileup rejection MB");
    spdPileupCutMB->AddEventTagFilterBit(9, kTRUE);
    
    AliReducedCompositeCut* eventCutMB = new AliReducedCompositeCut("eventCutMB", "", kTRUE);
    eventCutMB->AddCut(defaultCut);
    eventCutMB->AddCut(triggerCutMB);
    eventCutMB->AddCut(spdPileupCutMB);
    eventCutMB->AddCut(mvPileupCut);
    task->AddMBEventCut(eventCutMB);
  }
}

//_______________________________________________________________________________________________________________
void SetClusterCuts(AliReducedAnalysisJpsi2eeCorrelations* task) {
  AliReducedCaloClusterCut* clusterCut = new AliReducedCaloClusterCut("clusters","cluster selection");
  clusterCut->AddCut(AliReducedVarManager::kEMCALdetector, AliReducedCaloClusterInfo::kEMCAL, AliReducedCaloClusterInfo::kEMCAL);
  clusterCut->AddCut(AliReducedVarManager::kEMCALm02, 0.01, 0.70);
  task->AddClusterCut(clusterCut);
}

//_______________________________________________________________________________________________________________
void SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/, Int_t systematicsMode /*=kStandard*/) {

  // default cuts to be used in either cut variation
  //_____________________________________________________________________________________________________________
  AliReducedTrackCut* defaultCut = new AliReducedTrackCut("defaultCut", "Electron default selection");
  defaultCut->SetTrackFilterBit(0); // fQualityFlag bit 32 (see AddTask_laltenka_dst_correlations.C)
  defaultCut->AddCut(AliReducedVarManager::kPt,                 1.0,   1.0e+30);
  defaultCut->AddCut(AliReducedVarManager::kEta,               -0.9,   0.9);
  defaultCut->AddCut(AliReducedVarManager::kDcaXY,             -1.0,   1.0);
  defaultCut->AddCut(AliReducedVarManager::kDcaZ,              -3.0,   3.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCncls,           70.0, 161.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCchi2,            0.0,   4.0);
  defaultCut->AddCut(AliReducedVarManager::kITSchi2,            0.0,  30.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired,  6.0,   9.0);
  defaultCut->SetRejectKinks();
  defaultCut->SetRequestITSrefit();
  defaultCut->SetRequestTPCrefit();

  // standard electron cut
  //_____________________________________________________________________________________________________________
  if (systematicsMode==kStandard || systematicsMode==kSystematicsFull ||
      systematicsMode==kHadronSPDstrict) {
    AliReducedTrackCut* standardCut = (AliReducedTrackCut*)defaultCut->Clone("standardElectron");
    standardCut->SetRequestSPDany();
    standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -1.5, 3.0);
    if (triggerChoice>=kEMC7 || gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
      // EMCal/DCal trigger -> good EMCal runs
      standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.5, 1.e8);
      standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.5, 1.e8, kFALSE,
                          AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
      standardCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.8,  1.3, kFALSE,
                          AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
    } else {
      // MB/HM trigger
      standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.5, 1.e8);
      standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.5, 1.e8);
    }
    task->AddTrackCut(standardCut);
  }
  
  // PID strict electron cut
  //_____________________________________________________________________________________________________________
  if (systematicsMode==kElectronPIDstrict || systematicsMode==kSystematicsFull) {
    AliReducedTrackCut* PIDstrictCut = (AliReducedTrackCut*)defaultCut->Clone("PIDstrictElectron");
    PIDstrictCut->SetRequestSPDany();
    PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -1.0, 2.5);
    if (triggerChoice>=kEMC7 || gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
      // EMCal/DCal trigger -> good EMCal runs
      PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  4.0, 1.e8);
      PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    4.0, 1.e8, kFALSE,
                           AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
      PIDstrictCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.9, 1.2, kFALSE,
                           AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
    } else {
      // MB/HM trigger
      PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  4.0, 1.e8);
      PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    4.0, 1.e8);
    }
    task->AddTrackCut(PIDstrictCut);
  }
  
  // PID open electron cut
  //_____________________________________________________________________________________________________________
  if (systematicsMode==kElectronPIDopen || systematicsMode==kSystematicsFull) {
    AliReducedTrackCut* PIDopenCut = (AliReducedTrackCut*)defaultCut->Clone("PIDopenElectron");
    PIDopenCut->SetRequestSPDany();
    PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -2.0, 3.5);
    if (triggerChoice>=kEMC7 || gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
      // EMCal/DCal trigger -> good EMCal runs
      PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.0, 1.e8);
      PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.0, 1.e8, kFALSE,
                         AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
      PIDopenCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.7, 1.4, kFALSE,
                         AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
    } else {
      // MB/HM trigger
      PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.0, 1.e8);
      PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.0, 1.e8);
    }
    task->AddTrackCut(PIDopenCut);
  }
}

//_______________________________________________________________________________________________________________
void SetAssociatedTrackCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t systematicsMode /*=kStandard*/) {
  // standard hadron cut
  if (systematicsMode==kStandard || systematicsMode==kSystematicsFull ||
      systematicsMode==kElectronPIDstrict || systematicsMode==kElectronPIDopen) {
    AliReducedTrackCut* standardCut = new AliReducedTrackCut("standardHadron","Associated hadron selection standard");
    standardCut->SetTrackFilterBit(4); // fQualityFlag bit 36 (see AddTask_laltenka_dst_correlations.C) -> associated hadrons, ITS + TPC refit, no kinks
    task->AddAssociatedTrackCut(standardCut);
  }

  // SPD strict hadron cut
  if (systematicsMode==kHadronSPDstrict) {
    AliReducedTrackCut* SPDstrictCut = new AliReducedTrackCut("SPDstrictHadron","Associated hadron selection SPD any");
    SPDstrictCut->SetTrackFilterBit(5); // fQualityFlag bit 37 (see AddTask_laltenka_dst_correlations.C) -> associated hadrons, ITS + TPC refit, no kinks, SPD any
    task->AddAssociatedTrackCut(SPDstrictCut);
  }
}

//_______________________________________________________________________________________________________________
void SetTrackPrefilterCuts(AliReducedAnalysisJpsi2eeCorrelations* task) {
  AliReducedTrackCut* trackPrefilterCut = new AliReducedTrackCut("TrackPrefilterCut","Track prefilter selection");
  trackPrefilterCut->SetTrackFilterBit(1); // fQualityFlag bit 33 (see AddTask_laltenka_dst_correlations.C)
  trackPrefilterCut->AddCut(AliReducedVarManager::kPt, 0.7, 1.0e+30);
  task->AddPrefilterTrackCut(trackPrefilterCut);
}

//_______________________________________________________________________________________________________________
void SetPairCuts(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/, Int_t systematicsMode /*=kStandard*/) {
  
  // default pair cut
  AliReducedTrackCut* defaultPairCut = new AliReducedTrackCut("defaultPairCut", "Pair default selection");
  defaultPairCut->AddCut(AliReducedVarManager::kRap,  -0.9, 0.9);

  // standard pair cut
  if (triggerChoice==kEG1 || triggerChoice==kDG1 || triggerChoice==kEG1DG1 || gUseEG1DG1PidMC) {
    AliReducedTrackCut* standardPairCutLeg1 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg1");
    AliReducedTrackCut* standardPairCutLeg2 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg2");

    standardPairCutLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
    standardPairCutLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

    AliReducedCompositeCut* standardPairCut = new AliReducedCompositeCut("standardPair", "", kFALSE);
    standardPairCut->AddCut(standardPairCutLeg1);
    standardPairCut->AddCut(standardPairCutLeg2);
    task->AddPairCut(standardPairCut);
  } else if (triggerChoice==kEG2 || triggerChoice==kDG2 || triggerChoice==kEG2DG2 || gUseEG2DG2PidMC) {
    AliReducedTrackCut* standardPairCutLeg1 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg1");
    AliReducedTrackCut* standardPairCutLeg2 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg2");

    standardPairCutLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
    standardPairCutLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

    AliReducedCompositeCut* standardPairCut = new AliReducedCompositeCut("standardPair", "", kFALSE);
    standardPairCut->AddCut(standardPairCutLeg1);
    standardPairCut->AddCut(standardPairCutLeg2);
    task->AddPairCut(standardPairCut);
  } else {
    AliReducedTrackCut* standardPairCut = (AliReducedTrackCut*)defaultPairCut->Clone("standardPair");
    task->AddPairCut(standardPairCut);
  }

  // non-prompt pair cuts
  if (systematicsMode==kStandard) {
    TFile* cutFile = TFile::Open("/Users/lal035/Documents/PWGDQ/analysis_MC_decayLengthCut/decayLengthResolution_MB_MC_standardElectron_paper.root");

    // default cut with decay-length cut applied - 5% prompt contamination
    AliReducedTrackCut* defaultPairCut5perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut5perc");
    defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_5percent_fit"), 1.0e+30, kFALSE,
                                AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
    defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_5percent_fit"), 1.0e+30, kFALSE,
                                AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
    defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_5percent_fit"), 1.0e+30, kFALSE,
                                AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);
    
    // default cut with decay-length cut applied - 10% prompt contamination
    AliReducedTrackCut* defaultPairCut10perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut10perc");
    defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_10percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
    defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_10percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
    defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_10percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);

    // default cut with decay-length cut applied - 20% prompt contamination
    AliReducedTrackCut* defaultPairCut20perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut20perc");
    defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_20percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
    defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_20percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
    defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_20percent_fit"), 1.0e+30, kFALSE,
                                 AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                                 AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);
    
    // set pair cuts
    if (triggerChoice==kEG1 || triggerChoice==kDG1 || triggerChoice==kEG1DG1 || gUseEG1DG1PidMC) {
      // non-prompt, 5% prompt contamination
      AliReducedTrackCut* pairCutNonPromptLeg1 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg1");
      AliReducedTrackCut* pairCutNonPromptLeg2 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg2");

      pairCutNonPromptLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
      pairCutNonPromptLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt = new AliReducedCompositeCut("nonPromptPair5perc", "", kFALSE);
      pairCutNonPrompt->AddCut(pairCutNonPromptLeg1);
      pairCutNonPrompt->AddCut(pairCutNonPromptLeg2);
      task->AddPairCut(pairCutNonPrompt);

      // non-prompt, 10% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt1Leg1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg1");
      AliReducedTrackCut* pairCutNonPrompt1Leg2 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg2");

      pairCutNonPrompt1Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
      pairCutNonPrompt1Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt1 = new AliReducedCompositeCut("nonPromptPair10perc", "", kFALSE);
      pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg1);
      pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg2);
      task->AddPairCut(pairCutNonPrompt1);

      // non-prompt, 20% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt2Leg1 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg1");
      AliReducedTrackCut* pairCutNonPrompt2Leg2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg2");

      pairCutNonPrompt2Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
      pairCutNonPrompt2Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt2 = new AliReducedCompositeCut("nonPromptPair20perc", "", kFALSE);
      pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg1);
      pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg2);
      task->AddPairCut(pairCutNonPrompt2);
    } else if (triggerChoice==kEG2 || triggerChoice==kDG2 || triggerChoice==kEG2DG2 || gUseEG2DG2PidMC) {
      // non-prompt, 5% prompt contamination
      AliReducedTrackCut* pairCutNonPromptLeg1 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg1");
      AliReducedTrackCut* pairCutNonPromptLeg2 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg2");

      pairCutNonPromptLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
      pairCutNonPromptLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt = new AliReducedCompositeCut("nonPromptPair5perc", "", kFALSE);
      pairCutNonPrompt->AddCut(pairCutNonPromptLeg1);
      pairCutNonPrompt->AddCut(pairCutNonPromptLeg2);
      task->AddPairCut(pairCutNonPrompt);

      // non-prompt, 10% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt1Leg1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg1");
      AliReducedTrackCut* pairCutNonPrompt1Leg2 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg2");

      pairCutNonPrompt1Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
      pairCutNonPrompt1Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt1 = new AliReducedCompositeCut("nonPromptPair10perc", "", kFALSE);
      pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg1);
      pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg2);
      task->AddPairCut(pairCutNonPrompt1);

      // non-prompt, 20% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt2Leg1 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg1");
      AliReducedTrackCut* pairCutNonPrompt2Leg2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg2");

      pairCutNonPrompt2Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
      pairCutNonPrompt2Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

      AliReducedCompositeCut* pairCutNonPrompt2 = new AliReducedCompositeCut("nonPromptPair20perc", "", kFALSE);
      pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg1);
      pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg2);
      task->AddPairCut(pairCutNonPrompt2);
    } else {
      // non-prompt, 5% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt = (AliReducedTrackCut*)defaultPairCut5perc->Clone("nonPromptPair5perc");
      task->AddPairCut(pairCutNonPrompt);

      // non-prompt, 10% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("nonPromptPair10perc");
      task->AddPairCut(pairCutNonPrompt1);

      // non-prompt, 20% prompt contamination
      AliReducedTrackCut* pairCutNonPrompt2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("nonPromptPair20perc");
      task->AddPairCut(pairCutNonPrompt2);
    }
  }
}

//_______________________________________________________________________________________________________________
void SetPairPrefilterCuts(AliReducedAnalysisJpsi2eeCorrelations* task) {
  AliReducedVarCut* pairPrefilterCut = new AliReducedVarCut("PairPrefilterCut","Pair prefilter selection");
  pairPrefilterCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  task->AddPrefilterPairCut(pairPrefilterCut);
}

//_______________________________________________________________________________________________________________
void SetClusterTrackMatcher(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/) {
  AliReducedCaloClusterTrackMatcher* matcher = new AliReducedCaloClusterTrackMatcher();
  matcher->SetMaximumMatchingDistance(0.005);
  task->SetClusterTrackMatcher(matcher);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsLeg(AliReducedAnalysisJpsi2eeCorrelations* task) {
  // electron from inclusive J/psi
  AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons from inclusive Jpsi with MC truth");
  trueElectron->SetMCFilterBit(kJpsiDecayElectron);
  task->AddLegCandidateMCcut(trueElectron);
  
  // electron from non-prompt J/psi
  AliReducedTrackCut* trueElectronNonPrompt = new AliReducedTrackCut("TrueElectronNonPrompt", "reconstructed electrons from non-prompt Jpsi with MC truth");
  trueElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronNonPrompt);
  
  // electron from prompt J/psi
  AliReducedTrackCut* trueElectronPrompt = new AliReducedTrackCut("TrueElectronPrompt", "reconstructed electrons from prompt Jpsi with MC truth");
  trueElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronPrompt);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2eeCorrelations* task) {
  // inclusive J/psi
  AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth inclusive J/psi");
  mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
  mcTruthJpsi->SetApplyReweightMCpt(kFALSE); //reweight is applied only for the prompt (kTRUE will reweight twice)
  mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from inclusive J/psi");
  mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta,   -0.9, 0.9);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt,    1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);
  
  // non-prompt J/psi
  AliReducedTrackCut* mcTruthJpsiNonPrompt = new AliReducedTrackCut("mcTruthJpsiNonPrompt", "Pure MC truth non-prompt J/psi");
  mcTruthJpsiNonPrompt->SetMCFilterBit(kJpsiNonPrompt);
  mcTruthJpsiNonPrompt->SetApplyReweightMCpt(kFALSE); //no reweight
  mcTruthJpsiNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronNonPrompt = new AliReducedTrackCut("mcTruthJpsiElectronNonPrompt", "Pure MC truth electron from non-prompt J/psi");
  mcTruthJpsiElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kEta,  -0.9, 0.9);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kPt,   1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsiNonPrompt, mcTruthJpsiElectronNonPrompt);

  // prompt J/psi
  AliReducedTrackCut* mcTruthJpsiPrompt = new AliReducedTrackCut("mcTruthJpsiPrompt", "Pure MC truth prompt J/psi");
  mcTruthJpsiPrompt->SetMCFilterBit(kJpsiPrompt);
  mcTruthJpsiPrompt->SetApplyReweightMCpt(kTRUE); //Apply reweighting only for prompt J/psi
  mcTruthJpsiPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronPrompt = new AliReducedTrackCut("mcTruthJpsiElectronPrompt", "Pure MC truth electron from prompt J/psi");
  mcTruthJpsiElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kPt,  1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsiPrompt, mcTruthJpsiElectronPrompt);
}

//_______________________________________________________________________________________________________________
void SetupMixingHandlers(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t triggerChoice /*=0*/, Int_t systematicsMode /*=kStandard*/) {
  
  // jpsi signal extraction mixing handler
  AliMixingHandler* handler = task->GetMixingHandler();
  handler->SetPoolDepth(100);
  handler->SetMixingThreshold(1.0);
  handler->SetDownscaleEvents(1);
  handler->SetDownscaleTracks(1);
  handler->AddMixingVariable(AliReducedVarManager::kVtxZ, gNVtxZLimits, gVtxZLimits);
  
  // provide pair cuts to mixing handler
  AliReducedInfoCut*          pairCutME[10];
  for (Int_t i=0; i<10; i++)  pairCutME[i] = NULL;
  for (Int_t iCut=0; iCut<task->GetNPairCuts(); iCut++) {
    pairCutME[iCut] = (AliReducedInfoCut*)(task->GetPairCut(iCut))->Clone(task->GetPairCutName(iCut));
    handler->AddPairsCut(pairCutME[iCut]);
  }

  // correlation mixing handler
  AliMixingHandler* handler2 = task->GetCorrelationMixingHandler();
  if (!task->GetUseLikeSignPairs()) handler2->SetMixLikeSign(kFALSE);
  handler2->SetPoolDepth(100);
  handler2->SetMixingThreshold(1.0);
  handler2->SetDownscaleEvents(1);
  handler2->SetDownscaleTracks(1);
  handler2->AddMixingVariable(AliReducedVarManager::kVtxZ, gNVtxZLimits, gVtxZLimits);
}

//_______________________________________________________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t runningMode, TString runNumbers /*=""*/) {
  AliReducedVarManager::SetDefaultVarNames();
  DefineHistograms(task, runningMode, runNumbers);
  AliReducedVarManager::SetUseVars((Bool_t*)task->GetHistogramManager()->GetUsedVars());
}

//_______________________________________________________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2eeCorrelations* task, Int_t runningMode, TString runNumbers /*=""*/) {
  AliHistogramManager* man = task->GetHistogramManager();

  //-------------------------------------------------------------------------------------------------------------
  // binning
  //-------------------------------------------------------------------------------------------------------------

  // mass binning
  const Int_t kNMassBins = 71;      // [1.4, 4.2]
  Double_t    massBins[kNMassBins];
  for (Int_t i=0; i<kNMassBins; i++) massBins[i] = 1.4+i*0.04;

  const Int_t kNMassBinsCorr = 31;  // [1.4, 4.2]
  Double_t    massBinsCorr[kNMassBinsCorr];
  for (Int_t i=0; i<6; i++)   massBinsCorr[i]     = 1.4+i*0.2;  // 1.40 - 2.20
  for (Int_t i=0; i<20; i++)  massBinsCorr[i+6]   = 2.6+i*0.04; // 2.60 - 3.36
  for (Int_t i=0; i<5; i++)   massBinsCorr[i+26]  = 3.4+i*0.2;  // 3.40 - 4.20

  // pt binning
  const Int_t kNPtBins          = 11;             // [0., 40.]
  Double_t    ptBins[kNPtBins]  = {0.0, 1.0, 3.0, 5.0, 7.0, 8.0, 12.0, 15.0, 20.0, 30.0, 40.0};
  
  const Int_t kNPtBinsEff = 201;                  // [0., 40.]
  Double_t    ptBinsEff[kNPtBinsEff];
  for (Int_t i=0; i<kNPtBinsEff; i++) ptBinsEff[i] = 0. + i*0.2;
  
  const Int_t kNPtAssocBins               = 5;    // [0.15, 10.]
  Double_t    ptAssocBins[kNPtAssocBins]  = {0.15, 1.0, 3.0, 5.0, 10.0};

  // phi binning
  const Int_t kNPhiBins = 41;                     // [0., 2*pi]
  Double_t    phiBins[kNPhiBins];
  for (Int_t i=0; i<kNPhiBins; i++) phiBins[i] = i*TMath::Pi()/20;

  const Int_t kNDeltaPhiBins = 13;                // [0, pi]
  Double_t    deltaPhiBins[kNDeltaPhiBins];
  for (Int_t i=0; i<kNDeltaPhiBins; i++) deltaPhiBins[i] = i*TMath::Pi()/12.;
  
  // eta binning
  const Int_t kNEtaBins = 21;                     // [-1., 1.]
  Double_t    etaBins[kNEtaBins];
  for (Int_t i=0; i<kNEtaBins; i++) etaBins[i] = -1+i*0.1;
  
  const Int_t kNDeltaEtaBins                = 11; // [-3., 3.]
  Double_t    deltaEtaBins[kNDeltaEtaBins]  = {-3., -2., -1., -0.5, -0.25, 0., 0.25, 0.5, 1., 2., 3.};
  
  // SPD pair type
  const Int_t kNPairTypeSPDBins                   = 4;
  Double_t    pairTypeSPDBins[kNPairTypeSPDBins]  = {-0.5, 0.5, 1.5, 2.5};
  
  //-------------------------------------------------------------------------------------------------------------
  // run numbers
  //-------------------------------------------------------------------------------------------------------------
  Int_t     runNBins        = 0;
  Double_t  runHistRange[2] = {0.0,0.0};
  if (runNumbers.CompareTo("")) {
    TObjArray* runArr = runNumbers.Tokenize(";");
    runNBins = runArr->GetEntries();
    runHistRange[1] = (Double_t)runNBins;
  } else {
    cout << "WARNING in AddTask_laltenka_jpsi2ee_correlations::DefineHistograms(): Run ranges not defined!" << endl;
  }
  
  //-------------------------------------------------------------------------------------------------------------
  // histogram classes
  //-------------------------------------------------------------------------------------------------------------
  
  // event quantities
  TString histClasses = "";
  if (runningMode==kFull || runningMode==kQA) {
    histClasses += "Event_BeforeCuts;";
    histClasses += "Event_AfterCuts;";
  }
  
  // calorimeter quantities
  if (task->GetFillCaloClusterHistograms() && (runningMode==kFull || runningMode==kQA)) {
    histClasses += "CaloCluster_BeforeCuts;";
    for (Int_t i=0; i<task->GetNClusterCuts(); ++i) {
      TString cutName = task->GetClusterCutName(i);
      histClasses += Form("CaloCluster_%s;", cutName.Data());
    }
  }

  // J/psi MC gen. tracks
  if (task->GetRunOverMC()) {
    for (Int_t mcSel=0; mcSel<task->GetNJpsiMotherMCCuts(); ++mcSel) {
      histClasses += Form("%s_PureMCTruth_BeforeSelection;", task->GetJpsiMotherMCcutName(mcSel));
      histClasses += Form("%s_PureMCTruth_AfterSelection;", task->GetJpsiMotherMCcutName(mcSel));
    }
  }
  
  // electron track quantities
  if (task->GetLoopOverTracks() && (runningMode==kFull || runningMode==kQA)) {
    histClasses += "Track_BeforeCuts;";
    for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
      if (task->GetRunOverMC()) {
        for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) {
          histClasses += Form("Track_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
        }
      }
    }
  }
  
  // associated track quantities
  if (task->GetRunCorrelation() && (runningMode==kFull || runningMode==kQA)) {
    histClasses += "AssociatedTrack_BeforeCuts;";
    for (Int_t i=0; i<task->GetNAssociatedTrackCuts(); ++i) {
      TString assocCutName = task->GetAssociatedTrackCutName(i);
      histClasses += Form("AssociatedTrack_%s;", assocCutName.Data());
    }
  }
   
  // pair quantities
  if (task->GetRunPairing() && (runningMode==kFull || runningMode==kCorr)) {
    if (task->GetNPairCuts()>1) {
      for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
        for (Int_t j=0; j<task->GetNPairCuts(); ++j) {
          TString trackCutName  = task->GetTrackCutName(i);
          TString pairCutName   = task->GetPairCutName(j);
          if (!task->GetRunLikeSignPairing()) histClasses += Form("PairSEPM_%s_%s;",
                                                                  trackCutName.Data(), pairCutName.Data());
          else                                histClasses += Form("PairSEPP_%s_%s;PairSEPM_%s_%s;PairSEMM_%s_%s;",
                                                                  trackCutName.Data(), pairCutName.Data(),
                                                                  trackCutName.Data(), pairCutName.Data(),
                                                                  trackCutName.Data(), pairCutName.Data());
          if (!task->GetRunOverMC() && !pairCutName.Contains("nonPrompt"))  histClasses += Form("PairMEPP_%s_%s;PairMEPM_%s_%s;PairMEMM_%s_%s;",
                                                                                                trackCutName.Data(), pairCutName.Data(),
                                                                                                trackCutName.Data(), pairCutName.Data(),
                                                                                                trackCutName.Data(), pairCutName.Data());
          if (task->GetRunOverMC()) {
            for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) histClasses += Form("PairSEPM_%s_%s_%s;",
                                                                                                   trackCutName.Data(), pairCutName.Data(),
                                                                                                   task->GetLegCandidateMCcutName(mcSel));
          }
        }
      }
    } else {
      for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
        TString trackCutName  = task->GetTrackCutName(i);
        if (!task->GetRunLikeSignPairing()) histClasses += Form("PairSEPM_%s;", trackCutName.Data());
        else                                histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", trackCutName.Data(), trackCutName.Data(), trackCutName.Data());
        if (!task->GetRunOverMC())          histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", trackCutName.Data(), trackCutName.Data(), trackCutName.Data());
        if (task->GetRunOverMC()) {
          for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) histClasses += Form("PairSEPM_%s_%s;", trackCutName.Data(),
                                                                                                 task->GetLegCandidateMCcutName(mcSel));
        }
      }
    }
  }
  
  // correlation quantities
  if (task->GetRunCorrelation() && (runningMode==kFull || runningMode==kCorr)) {
    if (task->GetNPairCuts()>1) {
      for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
        for (Int_t j=0; j<task->GetNPairCuts(); j++) {
          TString elecCutName   = task->GetTrackCutName(i);
          TString assocCutName  = task->GetAssociatedTrackCutName(i);
          TString pairCutName   = task->GetPairCutName(j);
          if (!task->GetUseLikeSignPairs()) histClasses += Form("CorrSEPM_%s_%s_%s;",
                                                                elecCutName.Data(), assocCutName.Data(), pairCutName.Data());
          else                              histClasses += Form("CorrSEPP_%s_%s_%s;CorrSEPM_%s_%s_%s;CorrSEMM_%s_%s_%s;",
                                                                elecCutName.Data(), assocCutName.Data(), pairCutName.Data(),
                                                                elecCutName.Data(), assocCutName.Data(), pairCutName.Data(),
                                                                elecCutName.Data(), assocCutName.Data(), pairCutName.Data());
          if(task->GetRunCorrelationMixing()) {
            if (!task->GetUseLikeSignPairs()) histClasses += Form("CorrMEPM_%s_%s_%s;",
                                                                  elecCutName.Data(), assocCutName.Data(), pairCutName.Data());
            else                              histClasses += Form("CorrMEPP_%s_%s_%s;CorrMEPM_%s_%s_%s;CorrMEMM_%s_%s_%s;",
                                                                  elecCutName.Data(), assocCutName.Data(), pairCutName.Data(),
                                                                  elecCutName.Data(), assocCutName.Data(), pairCutName.Data(),
                                                                  elecCutName.Data(), assocCutName.Data(), pairCutName.Data());
          }
        }
      }
    } else {
      for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
        TString elecCutName   = task->GetTrackCutName(i);
        TString assocCutName  = task->GetAssociatedTrackCutName(i);
        if (!task->GetUseLikeSignPairs()) histClasses += Form("CorrSEPM_%s_%s;",
                                                              elecCutName.Data(), assocCutName.Data());
        else                              histClasses += Form("CorrSEPP_%s_%s;CorrSEPM_%s_%s;CorrSEMM_%s_%s;",
                                                              elecCutName.Data(), assocCutName.Data(),
                                                              elecCutName.Data(), assocCutName.Data(),
                                                              elecCutName.Data(), assocCutName.Data());
        if(task->GetRunCorrelationMixing()) {
          if (!task->GetUseLikeSignPairs()) histClasses += Form("CorrMEPM_%s_%s;",
                                                                elecCutName.Data(), assocCutName.Data());
          else                              histClasses += Form("CorrMEPP_%s_%s;CorrMEPM_%s_%s;CorrMEMM_%s_%s;",
                                                                elecCutName.Data(), assocCutName.Data(),
                                                                elecCutName.Data(), assocCutName.Data(),
                                                                elecCutName.Data(), assocCutName.Data());
        }
      }
    }
  }
  
  // add histograms according to class to histogram manager
  TString classesStr(histClasses);
  TObjArray* arr = classesStr.Tokenize(";");
  
  //-------------------------------------------------------------------------------------------------------------
  // histograms
  //-------------------------------------------------------------------------------------------------------------

  // loop over histogram classes and add histograms
  printf("INFO on AddTask_laltenka_jpsi2ee_correlations(): Histgram classes included in histogram manager\n");
  for (Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
      
    //-----------------------------------------------------------------------------------------------------------
    // event histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // run numbers
      man->AddHistogram(classStr.Data(), "RunNo", "Run numbers", kFALSE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // multiplicity
      man->AddHistogram(classStr.Data(), "SPDnTracklets", "SPD #tracklets in |#eta|<1.0", kFALSE,
                        200, 0.0, 200.,                               AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(), "SPDtracklets_RunNo_QA", "<SPD no. tracklets> vs. Run ID", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, 0.0, 200.,                               AliReducedVarManager::kSPDntracklets,
                        0,0,0,                                        AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(),"MultiplicityV0", "", kFALSE,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"MultiplicityV0_RunNo_QA", "<mult. V0> vs. Run ID", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROTotalMult,
                        0,0,0,                                        AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(),"MultiplicityV0A", "", kFALSE,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROATotalMult);
      man->AddHistogram(classStr.Data(),"MultiplicityV0A_RunNo_QA", "<mult. V0A> vs. Run ID", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROATotalMult,
                        0,0,0,                                        AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(),"MultiplicityV0AC", "", kFALSE,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROCTotalMult);
      man->AddHistogram(classStr.Data(),"MultiplicityV0C_RunNo_QA", "<mult. V0C> vs. Run ID", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROCTotalMult,
                        0,0,0,                                        AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // vertex histograms
      man->AddHistogram(classStr.Data(), "VtxX", "Vtx X", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(), "VtxX_RunNo_QA", "<Vtx X> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxX,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "VtxY", "Vtx Y", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(), "VtxY_RunNo_QA", "<Vtx Y> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxY,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "VtxZ", "Vtx Z", kFALSE,
                        150, -15.0, 15.0,                             AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(), "VtxZ_RunNo_QA", "<Vtx Z> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        150, -15.0, 15.0,                             AliReducedVarManager::kVtxZ,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "NVtxContributors_RunNo_QA", "number of vertex contributors vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        100, 0., 100.,                                AliReducedVarManager::kNVtxContributors,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // SPD information
      for (Int_t il=0; il<2; ++il) {
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), kFALSE,
                          100, 0., 500.,                              AliReducedVarManager::kSPDFiredChips+il);
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d_RunNo_QA",il+1), Form("SPD <#fired chips> in layer %d vs run number",il+1), kTRUE,
                          runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunID,
                          100, 0., 500.,                              AliReducedVarManager::kSPDFiredChips+il,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          runNumbers.Data());
        man->AddHistogram(classStr.Data(), Form("SPDfirecChips_layer%d_Vs_SPDnTracklets", il+1), Form("SPD fired chips in layer %d vs SPD tracklets",il+1), kFALSE,
                          200, 0., 200.,                              AliReducedVarManager::kSPDntracklets,
                          250, 0., 500.,                              AliReducedVarManager::kSPDFiredChips+il);
      }

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // calorimeter cluster histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("CaloCluster_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      man->AddHistogram(classStr.Data(), "Energy",          "Cluster energy",                 kFALSE, 500,  0.0,  50.0, AliReducedVarManager::kEMCALclusterEnergy);
      man->AddHistogram(classStr.Data(), "Dx",              "Cluster dx",                     kFALSE, 200, -0.2,   0.2, AliReducedVarManager::kEMCALclusterDx);
      man->AddHistogram(classStr.Data(), "Dz",              "Cluster dz",                     kFALSE, 200, -0.2,   0.2, AliReducedVarManager::kEMCALclusterDz);
      man->AddHistogram(classStr.Data(), "M02",             "Cluster M02",                    kFALSE, 100,  0.0,   1.0, AliReducedVarManager::kEMCALm02);
      man->AddHistogram(classStr.Data(), "M20",             "Cluster M20",                    kFALSE, 100,  0.0,   1.0, AliReducedVarManager::kEMCALm20);
      man->AddHistogram(classStr.Data(), "Detector",        "Detector",                       kFALSE,   3, -0.5,   2.5, AliReducedVarManager::kEMCALdetector);
      man->AddHistogram(classStr.Data(), "NCells",          "No. cells per cluster",          kFALSE, 100,  0.0, 100.0, AliReducedVarManager::kEMCALnCells);
      man->AddHistogram(classStr.Data(), "NMatchedTracks",  "No. matched tracks per cluster", kFALSE, 100,  0.0, 100.0, AliReducedVarManager::kEMCALnMatchedTracks);
      
      continue;
    }

    //-----------------------------------------------------------------------------------------------------------
    // track and associated track histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // momenta
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE,
                        200, 0., 100.,                                AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_RunNo_QA", "<p_{T}> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, 0., 100.,                                AliReducedVarManager::kPt,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // angles
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE,
                        90, -0.9, 0.9,                                AliReducedVarManager::kEta,
                        100, 0., 6.3,                                 AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_RunNo_QA", "<#eta> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        90, -0.9, 0.9,                                AliReducedVarManager::kEta,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "Phi_RunNo_QA", "<#varphi> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        100, 0., 6.3,                                 AliReducedVarManager::kPhi,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // skip rest for associated (i.e. base) tracks
      if (classStr.Contains("Associated")) continue;
      
      // DCA
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy vs. pT", kFALSE,
                        100, -2.0, 2.0,                               AliReducedVarManager::kDcaXY,
                        100, 0.0, 100.,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_RunNo_QA", "<DCAxy> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaXY,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz vs. pT", kFALSE,
                        100, -3.0, 3.0,                               AliReducedVarManager::kDcaZ,
                        100, 0.0, 100.,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_RunNo_QA", "<DCAz> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaZ,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // ITS
      man->AddHistogram(classStr.Data(), "ITSncls_RunNo_QA", "<ITS nclusters> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        7, -0.5, 6.5,                                 AliReducedVarManager::kITSncls,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "ITSchi2_RunNo_QA", "<ITS #chi^{2}> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, 0.0, 20.0,                               AliReducedVarManager::kITSchi2,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // TPC
      man->AddHistogram(classStr.Data(),  "TPCncls_RunNo_QA", "<TPC #cls> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        10, 0., 10.,                                  AliReducedVarManager::kTPCncls,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(),  "TPCnclsBitsFired_RunNo_QA", "<TPC #cls> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        10, 0., 10.,                                  AliReducedVarManager::kTPCNclusBitsFired,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "TPCchi2_RunNo_QA", "TPC <#chi^{2}> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        10, 0., 10.,                                  AliReducedVarManager::kTPCchi2,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // TPC PID
      man->AddHistogram(classStr.Data(), "TPCsignal_Pin", "TPC dE/dx vs. inner param P", kFALSE,
                        200, 0., 100.,                                AliReducedVarManager::kPin,
                        200, -0.5, 199.5,                             AliReducedVarManager::kTPCsignal);
      man->AddHistogram(classStr.Data(), "TPCsignal_Pin_Profile", "TPC dE/dx vs. inner param P", kTRUE,
                        200, 0., 100.,                                AliReducedVarManager::kPin,
                        200, -0.5, 199.5,                             AliReducedVarManager::kTPCsignal);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_Pin", "TPC N_{#sigma} electron vs. inner param P", kFALSE,
                        200, 0., 100.,                                AliReducedVarManager::kPin,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_Pin_Profile", "TPC N_{#sigma} electron vs. inner param P", kTRUE,
                        200, 0., 100.,                                AliReducedVarManager::kPin,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_Eta", "TPC N_{#sigma} electron vs. #eta", kFALSE,
                        36, -0.9, 0.9,                                AliReducedVarManager::kEta,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_Eta_Profile", "TPC N_{#sigma} electron vs. #eta", kTRUE,
                        36, -0.9, 0.9,                                AliReducedVarManager::kEta,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_RunNo_QA", "TPC N_{#sigma} electron vs. run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "TPCnsigPion_RunNo_QA", "TPC N_{#sigma} pion vs. run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -10.0, 10.0,                             AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "TPCnsigProton_RunNo_QA", "TPC N_{#sigma} proton vs. run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -10.0, 10.0,                             AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // matched calorimeter cluster
      man->AddHistogram(classStr.Data(), "EMCMatchedDeltaPhi", "#Delta#varphi track-cluster", kFALSE,
                        100,  -0.05,   0.05,                          AliReducedVarManager::kEMCALmatchedDeltaPhi);
      man->AddHistogram(classStr.Data(), "EMCMatchedDeltaEta", "#Delta#eta track-cluster", kFALSE,
                        100,  -0.05,   0.05,                          AliReducedVarManager::kEMCALmatchedDeltaEta);
      man->AddHistogram(classStr.Data(), "EMCMatchedDeltaPhi_EMCMatchedDeltaEta", "#Delta#eta vs. #Delta#varphi track-cluster", kFALSE,
                        100,  -0.05,   0.05,                          AliReducedVarManager::kEMCALmatchedDeltaPhi,
                        100,  -0.05,   0.05,                          AliReducedVarManager::kEMCALmatchedDeltaEta);
      man->AddHistogram(classStr.Data(), "EMCMatchedDistance", "rack-cluster distance", kFALSE,
                        50,  0.0,   0.05,                             AliReducedVarManager::kEMCALmatchedDistance);
      man->AddHistogram(classStr.Data(), "EMCMatchedNCells", "No. cells per cluster", kFALSE,
                        100,  0.,   100.,                             AliReducedVarManager::kEMCALmatchedNCells);
      man->AddHistogram(classStr.Data(), "EMCMatchedEnergy", "Energy from the calorimeter matched cluster", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kEMCALmatchedEnergy);
      man->AddHistogram(classStr.Data(), "EMCMatchedEnergy_Pt", "Energy from the calorimeter matched cluster vs pT", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt,
                        200, 0.0, 100.,                               AliReducedVarManager::kEMCALmatchedEnergy);
      
      // matched calorimeter cluster PID
      man->AddHistogram(classStr.Data(), "EMCMatchedEOverP", "E/P from the calorimeter matched cluster", kFALSE,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "EMCMatchedEOverP_Pt", "E/P from the calorimeter matched cluster vs. pT", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "EMCMatchedM02", "M02 from the calorimeter matched cluster", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM02);
      man->AddHistogram(classStr.Data(), "EMCMatchedM02_Pt", "M02 from the calorimeter matched cluster vs. pT", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM02);
      man->AddHistogram(classStr.Data(), "EMCMatchedEOverP_EMCMatchedM02", "E/P from the calorimeter matched cluster vs. M02", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM02,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "EMCMatchedM20", "M20 from the calorimeter matched cluster", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM20);
      man->AddHistogram(classStr.Data(), "EMCMatchedM20_Pt", "M20 from the calorimeter matched cluster vs. pT", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM20);
      man->AddHistogram(classStr.Data(), "EMCMatchedEOverP_EMCMatchedM20", "E/P from the calorimeter matched cluster vs. M20", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM20,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "EMCnSigmaElectron_Pt", "nSigma electron from the calorimeter matched cluster vs. pT", kFALSE,
                        100, 0.0, 50.,                                AliReducedVarManager::kPt,
                        100, -5.0, 5.0,                               AliReducedVarManager::kEMCALmatchedNSigmaElectron);
      man->AddHistogram(classStr.Data(), "EMCnSigmaElectron_EMCMatchedM02", "nSigma electron from the calorimeter matched cluster vs. M02", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM02,
                        100, -5.0, 5.0,                               AliReducedVarManager::kEMCALmatchedNSigmaElectron);
      man->AddHistogram(classStr.Data(), "EMCnSigmaElectron_EMCMatchedM20", "nSigma electron from the calorimeter matched cluster vs. M20", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM20,
                        100, -5.0, 5.0,                               AliReducedVarManager::kEMCALmatchedNSigmaElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_EMCMatchedM02", "nSigma electron from TPC vs. M02", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM02,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_EMCMatchedM20", "nSigma electron from TPC vs. M20", kFALSE,
                        100, 0.0, 1.0,                                AliReducedVarManager::kEMCALmatchedM20,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_EMCnSigmaElectron", "nSigma electron from TPC vs. nSigma electron from the calorimeter matched cluster", kFALSE,
                        100, -5.0, 5.0,                               AliReducedVarManager::kEMCALmatchedNSigmaElectron,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "EMCnSigmaElectron_EMCMatchedEOverP", "nSigma electron from the calorimeter matched cluster vs. E/p", kFALSE,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP,
                        100, -5.0, 5.0,                               AliReducedVarManager::kEMCALmatchedNSigmaElectron);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_EMCMatchedEOverP", "nSigma electron from TPC vs. E/p", kFALSE,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // MC histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("MCTruth")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      man->AddHistogram(classStr.Data(), "MassMC",      "MC mass",      kFALSE,  70,  1.4,   4.2, AliReducedVarManager::kMassMC);
      man->AddHistogram(classStr.Data(), "RapidityMC",  "MC rapidity",  kFALSE,  40, -1.0,   1.0, AliReducedVarManager::kRapMC);
      man->AddHistogram(classStr.Data(), "PtMC",        "p_{T} MC",     kFALSE, 200,  0.0, 100.0, AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "PhiMC",       "MC #varphi",   kFALSE, 100,  0.0,   6.3, AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "EtaMC",       "MC #eta",      kFALSE, 100, -1.5,   1.5, AliReducedVarManager::kEtaMC);

      // pair histogram for efficiency evaluation
      const Int_t kNVarsEff = 4;
      Int_t varsEff[kNVarsEff] = {AliReducedVarManager::kMassMC, AliReducedVarManager::kPtMC, AliReducedVarManager::kEtaMC, AliReducedVarManager::kPhiMC};
      TArrayD pairHistEffBinLimits[kNVarsEff];
      pairHistEffBinLimits[0] = TArrayD(kNMassBins,   massBins);
      pairHistEffBinLimits[1] = TArrayD(kNPtBinsEff,  ptBinsEff);
      pairHistEffBinLimits[2] = TArrayD(kNEtaBins,    etaBins);
      pairHistEffBinLimits[3] = TArrayD(kNPhiBins,    phiBins);
      man->AddHistogram(classStr.Data(), "PairInvMass_Eff", "Differential pair inv. mass", kNVarsEff, varsEff, pairHistEffBinLimits, 0x0, -1, kTRUE); // kTRUE = THnSparseF

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // pair histograms
    //-----------------------------------------------------------------------------------------------------------
    if (classStr.Contains("PairSE") || classStr.Contains("PairME")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // pair histogram
      const Int_t kNVars = 2;
      Int_t vars[kNVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt};
      TArrayD pairHistBinLimits[kNVars];
      pairHistBinLimits[0] = TArrayD(kNMassBins,          massBins);
      pairHistBinLimits[1] = TArrayD(kNPtBins,            ptBins);
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", kNVars, vars, pairHistBinLimits, 0x0, -1, kFALSE); // kTRUE = THnSparseF
      
      // pair histogram for efficiency evaluation
      if (task->GetRunOverMC()) {
        const Int_t kNVarsEff = 4;
        Int_t varsEff[kNVarsEff] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kEta, AliReducedVarManager::kPhi};
        TArrayD pairHistEffBinLimits[kNVarsEff];
        pairHistEffBinLimits[0] = TArrayD(kNMassBins,   massBins);
        pairHistEffBinLimits[1] = TArrayD(kNPtBinsEff,  ptBinsEff);
        pairHistEffBinLimits[2] = TArrayD(kNEtaBins,    etaBins);
        pairHistEffBinLimits[3] = TArrayD(kNPhiBins,    phiBins);
        man->AddHistogram(classStr.Data(), "PairInvMass_Eff", "Differential pair inv. mass", kNVarsEff, varsEff, pairHistEffBinLimits, 0x0, -1, kTRUE); // kTRUE = THnSparseF
      }

      // kinematic quantities
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt);
      if (gJpsiEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "Pt_efficiencyCorrected", "", kFALSE,
                          200, 0.0, 100.,                               AliReducedVarManager::kPt,
                          0, 0, 0,                                      AliReducedVarManager::kNothing,
                          0, 0, 0,                                      AliReducedVarManager::kNothing,
                          "", "", "",
                          AliReducedVarManager::kNothing, AliReducedVarManager::kOneOverPairEff);
      }
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE,
                        100, -1.0, 1.0,                               AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE,
                        100, 0.0, 6.3,                                AliReducedVarManager::kPhi);

      // cluster energy for leg
      man->AddHistogram(classStr.Data(), "Leg1EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 1)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+0);
      man->AddHistogram(classStr.Data(), "Leg2EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 2)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+1);
      man->AddHistogram(classStr.Data(), "Leg1EMCMatchedEnergy_Leg2EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 1 vs. leg 2)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+0,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+1);

      // skip rest for other than SE-PM
      if (!classStr.Contains("SEPM")) continue;

      // mean pT
      man->AddHistogram(classStr.Data(), "MeanPt_Mass", "<#it{p}_{T}> vs inv. mass", kTRUE,
                        70, 1.4, 4.2,                                 AliReducedVarManager::kMass,
                        400, 0.0, 40.0,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_Mass", "#it{p}_{T} vs inv. mass", kFALSE,
                        70, 1.4, 4.2,                                 AliReducedVarManager::kMass,
                        400, 0.0, 40.0,                               AliReducedVarManager::kPt);
      if (gJpsiEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "Pt_Mass_efficiencyCorrected", "#it{p}_{T} vs inv. mass", kFALSE,
                          70, 1.4, 4.2,                               AliReducedVarManager::kMass,
                          400, 0.0, 40.0,                             AliReducedVarManager::kPt,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          "", "", "",
                          AliReducedVarManager::kNothing, AliReducedVarManager::kOneOverPairEff);
      }

      // mass
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE,
                        125, 0.0, 5.0,                                AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Mass_RunNo_QA", "<Invariant mass> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        125, 0.0, 5.0,                                AliReducedVarManager::kMass,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // SPD pair type
      man->AddHistogram(classStr.Data(), "PairTypeSPD", "Pair type SPD", kFALSE,
                        kNPairTypeSPDBins-1, pairTypeSPDBins,         AliReducedVarManager::kPairTypeSPD);
      man->AddHistogram(classStr.Data(), "PairTypeSPD_Pt", "Pair type SPD vs. pt", kFALSE,
                        kNPtBins-1, ptBins,                           AliReducedVarManager::kPt,
                        kNPairTypeSPDBins-1, pairTypeSPDBins,         AliReducedVarManager::kPairTypeSPD);
      
      // pseudo-proper decay length
      man->AddHistogram(classStr.Data(), "DecayLength", "pseudo-proper decay length", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt", "pseudo-proper decay length vs. pt", kFALSE,
                        80, 0.0, 40.0,                                AliReducedVarManager::kPt,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_PairTypeSPD", "pseudo-proper decay length vs. pair type SPD", kFALSE,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kPairTypeSPD,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt_PairTypeSPD", "pseudo-proper decay length vs. pt vs. SPD pair type", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime,
                        80, 0.0, 40.0,                                AliReducedVarManager::kPt,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kPairTypeSPD);

      // chi2
      man->AddHistogram(classStr.Data(), "Leg1TPCchi2_Leg2TPCchi2", "", kFALSE,
                        100, 0.0, 6.0,                                AliReducedVarManager::kPairLegTPCchi2,
                        100, 0.0, 6.0,                                AliReducedVarManager::kPairLegTPCchi2+1);
      man->AddHistogram(classStr.Data(), "Leg1ITSchi2_Leg2ITSchi2", "", kFALSE,
                        100, 0.0, 6.0,                                AliReducedVarManager::kPairLegITSchi2,
                        100, 0.0, 6.0,                                AliReducedVarManager::kPairLegITSchi2+1);

      continue;
    }

    //-----------------------------------------------------------------------------------------------------------
    // correlation histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("CorrSE") || classStr.Contains("CorrME")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // correlation histogram
      const Int_t kNVarsCorr            = 5;
      Int_t       varsCorr[kNVarsCorr]  = {AliReducedVarManager::kMass, AliReducedVarManager::kTriggerPt, AliReducedVarManager::kAssociatedPt, AliReducedVarManager::kDeltaPhiSym, AliReducedVarManager::kDeltaEta};
      TArrayD corrHistBinLimits[kNVarsCorr];
      corrHistBinLimits[0] = TArrayD(kNMassBinsCorr,      massBinsCorr);
      corrHistBinLimits[1] = TArrayD(kNPtBins,            ptBins);
      corrHistBinLimits[2] = TArrayD(kNPtAssocBins,       ptAssocBins);
      corrHistBinLimits[3] = TArrayD(kNDeltaPhiBins,      deltaPhiBins);
      corrHistBinLimits[4] = TArrayD(kNDeltaEtaBins,      deltaEtaBins);
      man->AddHistogram(classStr.Data(), "CorrHist", "Differential correlation histogram", kNVarsCorr, varsCorr, corrHistBinLimits, 0x0, -1, kTRUE); // kTRUE = THnSparseF

      // correlation histogram, using hadron efficiency from map
      if (gJpsiEfficiencyUsed && gHadronEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "CorrHist_efficiencyCorrected", "Differential correlation histogram, efficiency corrected",
                          kNVarsCorr, varsCorr, corrHistBinLimits, 0x0, AliReducedVarManager::kOneOverTriggerEffTimesAssocHadronEff, kTRUE); // kTRUE = THnSparseF
      }

      // skip rest for SE/ME-PP/MM
      if (classStr.Contains("SEPP") || classStr.Contains("SEMM") ||
          classStr.Contains("MEPP") || classStr.Contains("MEMM")) continue;

      // J/psi candidates
      man->AddHistogram(classStr.Data(), "NJpsiCandidates", "number of J/psi candidates per event", kFALSE,
                        30, 0., 30,                                   AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(), "NJpsiCandidates_RunNo_QA", "<number of J/psi candidates> vs. run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        30, 0., 30.,                                  AliReducedVarManager::kNpairsSelected,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // associated tracks
      man->AddHistogram(classStr.Data(), "NAssociatedTracks", "number of associated tracks per event", kFALSE,
                        300, 0., 300.,                                AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(), "NAssociatedTracks_RunNo_QA", "<number of associated tracks> vs. run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        300, 0., 300.,                                AliReducedVarManager::kNtracksAnalyzed,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // pseudo-proper decay length of pair
      man->AddHistogram(classStr.Data(), "DecayLength", "pseudo-proper decay length", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kTriggerPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt", "pseudo-proper decay length vs. pt", kFALSE,
                        80, 0.0, 40.0,                                AliReducedVarManager::kTriggerPt,
                        800, -2.0, 2.0,                               AliReducedVarManager::kTriggerPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_PairTypeSPD", "pseudo-proper decay length vs. pair type SPD", kFALSE,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kTriggerPairTypeSPD,
                        800, -2.0, 2.0,                               AliReducedVarManager::kTriggerPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt_PairTypeSPD", "pseudo-proper decay length vs. pt vs. SPD pair type", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kTriggerPseudoProperDecayTime,
                        80, 0.0, 40.0,                                AliReducedVarManager::kTriggerPt,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kTriggerPairTypeSPD);
      
      // correlation QA histograms
      man->AddHistogram(classStr.Data(), "PtJpsi", "", kFALSE,
                          200, 0.0, 100.,                             AliReducedVarManager::kTriggerPt);
      if (gJpsiEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "PtJpsi_efficiencyCorrected", "", kFALSE,
                          200, 0.0, 100.,                             AliReducedVarManager::kTriggerPt,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          "", "", "",
                          AliReducedVarManager::kNothing, AliReducedVarManager::kOneOverTriggerEff);
      }
      man->AddHistogram(classStr.Data(), "PtAssoc", "", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kAssociatedPt);
      if (gHadronEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "PtAssoc_efficiencyCorrected", "", kFALSE,
                          200, 0.0, 100.,                             AliReducedVarManager::kAssociatedPt,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          "", "", "",
                          AliReducedVarManager::kNothing, AliReducedVarManager::kOneOverAssocHadronEff);
      }
      man->AddHistogram(classStr.Data(), "PtJpsi_PtAssoc", "", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kAssociatedPt,
                        200, 0.0, 100.,                               AliReducedVarManager::kTriggerPt);
      if (gJpsiEfficiencyUsed && gHadronEfficiencyUsed) {
        man->AddHistogram(classStr.Data(), "PtJpsi_PtAssoc_efficiencyCorrected", "", kFALSE,
                          200, 0.0, 100.,                             AliReducedVarManager::kAssociatedPt,
                          200, 0.0, 100.,                             AliReducedVarManager::kTriggerPt,
                          0, 0, 0,                                    AliReducedVarManager::kNothing,
                          "", "", "",
                          AliReducedVarManager::kNothing, AliReducedVarManager::kOneOverTriggerEffTimesAssocHadronEff);
      }
      man->AddHistogram(classStr.Data(), "DeltaPhi_PtAssoc", "", kFALSE,
                        kNPtAssocBins-1, ptAssocBins,                 AliReducedVarManager::kAssociatedPt,
                        kNDeltaPhiBins-1, deltaPhiBins,               AliReducedVarManager::kDeltaPhiSym);
      man->AddHistogram(classStr.Data(), "DeltaEta_PtAssoc", "", kFALSE,
                        kNPtAssocBins-1, ptAssocBins,                 AliReducedVarManager::kAssociatedPt,
                        kNDeltaEtaBins-1, deltaEtaBins,               AliReducedVarManager::kDeltaEta);

      continue;
    }
  }
}
