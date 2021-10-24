//===============================================================================================================
// addtask for associated hadron efficiency extraction in pp 13TeV (last updated: 05/11/2018)
//
// NOTE:  All settings (cuts, binning, etc.) must be identical to AddTask_laltenka_jpsi2ee_correlations
//        in order to represent the efficiency correctly!
//
//===============================================================================================================

// MC signals, in order of appearance in AddTask_laltenka_dst_assocHadronEff.C
enum MCFilters {
  kPionPrimMC=0,
  kPionSecMC,
  kPionSecMatMC,
  kKaonPrimMC,
  kKaonSecMC,
  kKaonSecMatMC,
  kProtonPrimMC,
  kProtonSecMC,
  kProtonSecMatMC,
  kElectronPrimMC,
  kElectronSecMC,
  kElectronSecMatMC,
  kMuonPrimMC,
  kMuonSecMC,
  kMuonSecMatMC
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
  kEGDGmerged,
  kNTriggerModes=11
};

// forward declaration
TString SetRunNumbers(TString prod="");
void    SetEventCuts(AliReducedAnalysisSingleTrack* task, Int_t triggerChoice=0);
void    SetAssociatedTrackCuts(AliReducedAnalysisSingleTrack* task);
void    SetMCSignalCuts(AliReducedAnalysisSingleTrack* task);
void    SetupHistogramManager(AliReducedAnalysisSingleTrack* task, TString prod="", TString runNumbers="");
void    DefineHistograms(AliReducedAnalysisSingleTrack* task, TString prod="", TString runNumbers="");

//_______________________________________________________________________________________________________________
AliAnalysisTask* AddTask_laltenka_assocHadronEfficiency(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prodString="") {
  //
  // arguments:
  //  - isAliroot = kTRUE for ESD/AOD analysis, = KFALSE for reduced trees
  //  - runMode = 1 for AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents, = 2 for AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree
  //
  printf("INFO on AddTask_laltenka_assocHadronEfficiency(): (isAliRoot, runMode) :: (%d,%d)\n", isAliRoot, runMode);

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
    if (!triggerString.CompareTo("MB"))         triggerChoice = kMB;          // AliReducedVarManager::kINT7
    if (!triggerString.CompareTo("V0HM"))       triggerChoice = kV0HM;        // AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD, V0  (trigger class)
    if (!triggerString.CompareTo("SPDHM"))      triggerChoice = kSPDHM;       // AliReducedVarManager::kHighMultV0+AliReducedVarManager::kHighMultSPD, SPD (trigger class)
    if (!triggerString.CompareTo("EMC7"))       triggerChoice = kEMC7;        // AliReducedVarManager::kEMC7
    if (!triggerString.CompareTo("EG1"))        triggerChoice = kEG1;         // AliReducedVarManager::kEMCEGA, EG1
    if (!triggerString.CompareTo("EG2"))        triggerChoice = kEG2;         // AliReducedVarManager::kEMCEGA, EG2
    if (!triggerString.CompareTo("DG1"))        triggerChoice = kDG1;         // AliReducedVarManager::kEMCEGA, DG1
    if (!triggerString.CompareTo("DG2"))        triggerChoice = kDG2;         // AliReducedVarManager::kEMCEGA, DG2
    if (!triggerString.CompareTo("EG1DG1"))     triggerChoice = kEG1DG1;      // AliReducedVarManager::kEMCEGA, EG1 | DG1
    if (!triggerString.CompareTo("EG2DG2"))     triggerChoice = kEG2DG2;      // AliReducedVarManager::kEMCEGA, EG2 | DG2
    if (!triggerString.CompareTo("EGDGmerged")) triggerChoice = kEGDGmerged;  // AliReducedVarManager::kEMCEGA, (EG1 | DG1) | (EG2 | DG2)
  }
  printf("INFO on AddTask_laltenka_assocHadronEfficiency(): trigger choice = %d (%s)\n", triggerChoice, triggerString.Data());

  // set trigger to MB for MC
  if (setMC && triggerChoice!=kMB) {
    printf("INFO on AddTask_laltenka_assocHadronEfficiency(): trigger choice (%s) not MB but must be MB for MC, setting trigger to MB... \n", triggerString.Data());
    triggerChoice = kMB;
  }
  
  // initialize analysis task
  //-----------------------------------------------------------------------------------------------------------
  AliReducedAnalysisSingleTrack* assocHadEffAnalysis = new AliReducedAnalysisSingleTrack("assocHadEffAnalysis", "Associated hadron efficiency analysis");
  assocHadEffAnalysis->Init();
  assocHadEffAnalysis->SetRunOverMC(setMC);

  // set run numbers for runwise histograms
  //-----------------------------------------------------------------------------------------------------------
  TString runNumbers = SetRunNumbers(prod);

  // set analysis and prefilter cuts
  //-----------------------------------------------------------------------------------------------------------
  SetEventCuts(assocHadEffAnalysis, triggerChoice);
  SetAssociatedTrackCuts(assocHadEffAnalysis);

  // set MC signal selection
  //-----------------------------------------------------------------------------------------------------------
  SetMCSignalCuts(assocHadEffAnalysis);

  // define histograms histograms
  //-----------------------------------------------------------------------------------------------------------
  SetupHistogramManager(assocHadEffAnalysis, prod, runNumbers);

  // initialize wrapper AliAnalysisTask
  // (in order to run AliReducedAnalysisJpsi2ee in an aliroot analysis train )
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(assocHadEffAnalysis);

  // intercept isAliRoot=kFALSE (nothing to be done yet)
  //-----------------------------------------------------------------------------------------------------------
  if (!isAliRoot) return 0;
  
  // get analysis manager
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_laltenka_assocHadronEfficiency", "No analysis manager found.");
    return 0;
  }

  // get data container
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisDataContainer* cReducedEvent = NULL;
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
    printf("INFO on AddTask_laltenka_assocHadronEfficiency(): use on the fly events\n");
    cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
    if (!cReducedEvent) {
      printf("ERROR: In AddTask_laltenka_assocHadronEfficiency(), ReducedEvent exchange container could not be found!\n");
      return 0;
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
    printf("ERROR: In AddTask_laltenka_assocHadronEfficiency(), runMode %d not defined!\n", runMode);
    return 0;
  }
  
  // connect output data containers
  //-----------------------------------------------------------------------------------------------------------
  TString                                   outputName = "assocHadronEfficiencyAnalysisHistograms.root";
  if (triggerString.CompareTo(""))          outputName = Form("assocHadronEfficiencyAnalysisHistograms_%s.root", triggerString.Data());
  if (assocHadEffAnalysis->GetRunOverMC())  outputName.ReplaceAll(".root", "_MC.root");
  printf("INFO on AddTask_laltenka_assocHadronEfficiency(): output container name = %s\n", outputName.Data());
  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("histos", THashList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  mgr->ConnectOutput(task, 1, cOutputHist);
  
  // done
  //-----------------------------------------------------------------------------------------------------------
  return task;
}

//_____________________________________________________________________________________________________________
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
void SetEventCuts(AliReducedAnalysisSingleTrack* task, Int_t triggerChoice /*=0*/) {
  
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
}

//_______________________________________________________________________________________________________________
void SetAssociatedTrackCuts(AliReducedAnalysisSingleTrack* task) {
  // standard cut (ITS + TPC refit + no kinks)
  AliReducedTrackCut* hadronCut = new AliReducedTrackCut("standardHadron", "Associated hadron selection standard");
  hadronCut->SetTrackFilterBit(2); // NOTE: fQualityFlag bit 34 (see AddTask_laltenka_dst_assocHadronEff.C)
  task->AddTrackCut(hadronCut);

  // standard cut + SPD any
  AliReducedTrackCut* SPDstrictCut = new AliReducedTrackCut("SPDstrictHadron", "Associated hadron selection SPD any");
  SPDstrictCut->SetTrackFilterBit(3); // NOTE: fQualityFlag bit 35 (see AddTask_laltenka_dst_assocHadronEff.C)
  task->AddTrackCut(SPDstrictCut);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCuts(AliReducedAnalysisSingleTrack* task) {
  // pions
  AliReducedTrackCut* truePionPrim = new AliReducedTrackCut("TruePionPrimary", "primary pions MC truth");
  truePionPrim->SetMCFilterBit(kPionPrimMC);
  truePionPrim->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  truePionPrim->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(truePionPrim);

  AliReducedTrackCut* truePionSec = new AliReducedTrackCut("TruePionSecondary", "secondary pions MC truth");
  truePionSec->SetMCFilterBit(kPionSecMC);
  truePionSec->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  truePionSec->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(truePionSec);

  AliReducedTrackCut* truePionSecMat = new AliReducedTrackCut("TruePionSecondaryFromMaterial", "secondary pions from material MC truth");
  truePionSecMat->SetMCFilterBit(kPionSecMatMC);
  truePionSecMat->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  truePionSecMat->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(truePionSecMat);

  // kaons
  AliReducedTrackCut* trueKaonPrim = new AliReducedTrackCut("TrueKaonPrimary", "primary kaons MC truth");
  trueKaonPrim->SetMCFilterBit(kKaonPrimMC);
  trueKaonPrim->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueKaonPrim->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueKaonPrim);

  AliReducedTrackCut* trueKaonSec = new AliReducedTrackCut("TrueKaonSecondary", "secondary kaons MC truth");
  trueKaonSec->SetMCFilterBit(kKaonSecMC);
  trueKaonSec->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueKaonSec->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueKaonSec);

  AliReducedTrackCut* trueKaonSecMat = new AliReducedTrackCut("TrueKaonSecondaryFromMaterial", "secondary kaons from material MC truth");
  trueKaonSecMat->SetMCFilterBit(kKaonSecMatMC);
  trueKaonSecMat->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueKaonSecMat->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueKaonSecMat);

  // protons
  AliReducedTrackCut* trueProtonPrim = new AliReducedTrackCut("TrueProtonPrimary", "primary protons MC truth");
  trueProtonPrim->SetMCFilterBit(kProtonPrimMC);
  trueProtonPrim->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueProtonPrim->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueProtonPrim);

  AliReducedTrackCut* trueProtonSec = new AliReducedTrackCut("TrueProtonSecondary", "secondary protons MC truth");
  trueProtonSec->SetMCFilterBit(kProtonSecMC);
  trueProtonSec->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueProtonSec->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueProtonSec);

  AliReducedTrackCut* trueProtonSecMat = new AliReducedTrackCut("TrueProtonSecondaryFromMaterial", "secondary protons from material MC truth");
  trueProtonSecMat->SetMCFilterBit(kProtonSecMatMC);
  trueProtonSecMat->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueProtonSecMat->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueProtonSecMat);

  // electrons
  AliReducedTrackCut* trueElectronPrim = new AliReducedTrackCut("TrueElectronPrimary", "primary electrons MC truth");
  trueElectronPrim->SetMCFilterBit(kElectronPrimMC);
  trueElectronPrim->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueElectronPrim->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueElectronPrim);

  AliReducedTrackCut* trueElectronSec = new AliReducedTrackCut("TrueElectronSecondary", "secondary electrons MC truth");
  trueElectronSec->SetMCFilterBit(kElectronSecMC);
  trueElectronSec->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueElectronSec->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueElectronSec);

  AliReducedTrackCut* trueElectronSecMat = new AliReducedTrackCut("TrueElectronSecondaryFromMaterial", "secondary electrons from material MC truth");
  trueElectronSecMat->SetMCFilterBit(kElectronSecMatMC);
  trueElectronSecMat->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueElectronSecMat->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueElectronSecMat);

  // muons
  AliReducedTrackCut* trueMuonPrim = new AliReducedTrackCut("TrueMuonPrimary", "primary muons MC truth");
  trueMuonPrim->SetMCFilterBit(kMuonPrimMC);
  trueMuonPrim->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueMuonPrim->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueMuonPrim);

  AliReducedTrackCut* trueMuonSec = new AliReducedTrackCut("TrueMuonSecondary", "secondary muons MC truth");
  trueMuonSec->SetMCFilterBit(kMuonSecMC);
  trueMuonSec->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueMuonSec->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueMuonSec);

  AliReducedTrackCut* trueMuonSecMat = new AliReducedTrackCut("TrueMuonSecondaryFromMaterial", "secondary muons from material MC truth");
  trueMuonSecMat->SetMCFilterBit(kMuonSecMatMC);
  trueMuonSecMat->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trueMuonSecMat->AddCut(AliReducedVarManager::kEta, 0.0, 0.0, kTRUE);
  task->AddMCSignalCut(trueMuonSecMat);
}

//_______________________________________________________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisSingleTrack* task, TString prod /*=""*/, TString runNumbers /*=""*/) {
  AliReducedVarManager::SetDefaultVarNames();
  DefineHistograms(task, prod, runNumbers);
  AliReducedVarManager::SetUseVars((Bool_t*)task->GetHistogramManager()->GetUsedVars());
}

//_______________________________________________________________________________________________________________
void DefineHistograms(AliReducedAnalysisSingleTrack* task, TString prod /*=""*/, TString runNumbers /*=""*/) {
  AliHistogramManager* man = task->GetHistogramManager();
  
  //-------------------------------------------------------------------------------------------------------------
  // binning
  //-------------------------------------------------------------------------------------------------------------

  // pT binning
  const Int_t nPtBins = 401; // [0., 20.]
  Double_t    ptBins[nPtBins];
  for (Int_t i=0; i<nPtBins; i++) ptBins[i] = i*0.05;

  // eta binning
  const Int_t nEtaBins = 21;
  Double_t    etaBins[nEtaBins];
  for (Int_t i=0; i<nEtaBins; i++) etaBins[i] = -1.+i*0.1;

  // phi binning
  const Int_t nPhiBins = 41;                     // [0., 2*pi]
  Double_t    phiBins[nPhiBins];
  for (Int_t i=0; i<nPhiBins; i++) phiBins[i] = i*TMath::Pi()/20;
  
  // dca binning
  const Int_t nDcaBins = 201;
  Double_t    dcaBins[nDcaBins];
  for (Int_t i=0; i<nDcaBins; i++) dcaBins[i] = -4.+i*0.04;
  
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
    cout << "WARNING in AddTask_laltenka_assocHadronEfficiency::DefineHistograms(), production " << prod.Data() << " not known, run ranges not defined!" << endl;
  }

  //-------------------------------------------------------------------------------------------------------------
  // histogram classes
  //-------------------------------------------------------------------------------------------------------------
  TString histClasses = "";

  // event quantities
  histClasses += "Event_BeforeCuts;";
  histClasses += "Event_AfterCuts;";

  // track quantites
  histClasses += "Track_BeforeCuts;";
  for (Int_t iCut=0; iCut<task->GetNTrackCuts(); ++iCut) {
    TString cutName = task->GetTrackCutName(iCut);
    histClasses += Form("Track_%s;", cutName.Data());
    if (task->GetRunOverMC()) {
      for (Int_t iMC=0; iMC<task->GetNMCSignalCuts(); ++iMC) {
        TString cutNameMC = task->GetMCSignalCutName(iMC);
        histClasses += Form("Track_%s_%s;", cutName.Data(), cutNameMC.Data());
      }
    }
  }
  
  // pure MC truth quantities
  if (task->GetRunOverMC()) {
    for (Int_t iMC=0; iMC<task->GetNMCSignalCuts(); ++iMC) {
      TString cutNameMC = task->GetMCSignalCutName(iMC);
      histClasses += Form("%s_PureMCTruth;", cutNameMC.Data());
    }
  }
  
  //-------------------------------------------------------------------------------------------------------------
  // histograms
  //-------------------------------------------------------------------------------------------------------------

  // add histograms according to class to histogram manager
  TString classesStr(histClasses);
  TObjArray* arr = classesStr.Tokenize(";");
  
  // loop over histogram classes and add histograms
  printf("INFO on AddTask_laltenka_assocHadronEfficiency(): Histgram classes included in histogram manager\n");
  for (Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();

    //-------------------------------------------------------------------------------------------------------------
    // event histograms
    //-------------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // run numbers
      man->AddHistogram(classStr.Data(), "RunNo", "Run numbers", kFALSE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // vertex histograms
      man->AddHistogram(classStr.Data(), "VtxX", "Vtx X", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(), "VtxY", "Vtx Y", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(), "VtxZ", "Vtx Z", kFALSE,
                        150, -15.0, 15.0,                             AliReducedVarManager::kVtxZ);

      continue;
    }
    
    //-------------------------------------------------------------------------------------------------------------
    // track histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // kinematical quantities
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE,
                        nPtBins-1, ptBins,                            AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_RunNo_QA", "<p_{T}> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, 0., 100.,                                AliReducedVarManager::kPt,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE,
                        nEtaBins-1,  etaBins,                         AliReducedVarManager::kEta,
                        nPhiBins-1,  phiBins,                         AliReducedVarManager::kPhi);
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

      // DCA
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE,
                        nDcaBins-1, dcaBins,                          AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy vs. pT", kFALSE,
                        nDcaBins-1, dcaBins,                          AliReducedVarManager::kDcaXY,
                        nPtBins-1, ptBins,                            AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_RunNo_QA", "<DCAxy> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaXY,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE,
                        nDcaBins-1, dcaBins,                          AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz vs. pT", kFALSE,
                        nDcaBins-1, dcaBins,                          AliReducedVarManager::kDcaZ,
                        nPtBins-1, ptBins,                            AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_RunNo_QA", "<DCAz> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaZ,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      
      // histogram for efficiency calculation (pT, eta, phi)
      const Int_t nVars = 3;
      Int_t vars[nVars] = {AliReducedVarManager::kPt, AliReducedVarManager::kEta, AliReducedVarManager::kPhi};
      TArrayD histBinLimits[nVars];
      histBinLimits[0] = TArrayD(nPtBins,   ptBins);
      histBinLimits[1] = TArrayD(nEtaBins,  etaBins);
      histBinLimits[2] = TArrayD(nPhiBins,  phiBins);
      man->AddHistogram(classStr.Data(), "TrackHist", "Differential track histogram", nVars, vars, histBinLimits, 0x0, -1, kTRUE);

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // MC gen. histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("PureMCTruth")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // kinematical quantities
      man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE,
                        nPtBins-1, ptBins,                            AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "PtMC_RunNo_QA", "<p_{T}> MC vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        200, 0., 100.,                                AliReducedVarManager::kPtMC,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "EtaMC_PhiMC", "", kFALSE,
                        nEtaBins-1,  etaBins,                         AliReducedVarManager::kEtaMC,
                        nPhiBins-1,  phiBins,                         AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "EtaMC_RunNo_QA", "<#eta> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        90, -0.9, 0.9,                                AliReducedVarManager::kEtaMC,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());
      man->AddHistogram(classStr.Data(), "PhiMC_RunNo_QA", "<#varphi> vs run number", kTRUE,
                        runNBins, runHistRange[0], runHistRange[1],   AliReducedVarManager::kRunID,
                        100, 0., 6.3,                                 AliReducedVarManager::kPhiMC,
                        0, 0, 0,                                      AliReducedVarManager::kNothing,
                        runNumbers.Data());

      // histogram for efficiency calculation (pT, eta, phi)
      const Int_t nVarsMC = 3;
      Int_t varsMC[nVarsMC] = {AliReducedVarManager::kPtMC, AliReducedVarManager::kEtaMC, AliReducedVarManager::kPhiMC};
      TArrayD histBinLimitsMC[nVarsMC];
      histBinLimitsMC[0] = TArrayD(nPtBins,   ptBins);
      histBinLimitsMC[1] = TArrayD(nEtaBins,  etaBins);
      histBinLimitsMC[2] = TArrayD(nPhiBins,  phiBins);
      man->AddHistogram(classStr.Data(), "TrackHist", "Differential track histogram", nVarsMC, varsMC, histBinLimitsMC, 0x0, -1, kTRUE);

      continue;
    }
  }
}
