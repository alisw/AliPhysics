{
  /*

    //  < 1 are private
    // run >= 100 < 1000 are real data
    // runs > 1000 are sim data
    1 = LHC09a1 pt hard 15 - 50 GeV 
    2 = LHC09a2 pt hard 50 -100 GeV
    3 = LHC09a3 pt hard > 100
    100 = LHC09d 0.9 GeV real with TPC
    200 = LHC10b 7 TeV real with TPC  
    202 = LHC10b 7 TeV real with TPC  (pass2) 
    300 = LHC10c 900 GeV real with TPC
    302 = LHC10c 900 GeV  pass2 real with TPC
    400 = LHC10c 7000 GeV real with TPC
    402 = LHC10c 7000 GeV pass2  real with TPC
    500 = LHC10d 7000 GeV Real
    502 = LHC10d 7000 GeV Real pass2
    801 = LHC10h Pb+Pb real
    802 = LHC10h Pb+Pb real pass2
   1012 = LHC10a12 0.9 GeV MC Pythiy D6T  
   1014 = LHC10a14 0.9 GeV MC ojet D6T  
   2001 = LHC10b1 7 TeV Phojet with real OCDB
   2002 = LHC10b2 7 TeV Pythia Pergiua-0
   4001 = LHC10d1 pp, Pythia6 Perugia-0, 0.5T, 7000GeV, // misses some runs since i a different direcotory LHC10c9
   4002 = LHC10d2 pp, Phojet, 0.5T, 7000GeV
   4004 = LHC10d4 pp, Perugis, 0.5T, 7000GeV
   5014 = LHC10e14 pp, Jet+Jet different p_T hard bins 0.5T, 7000GeV
   6006 = LHC10f6 pp, Phojet Nachor Runs LHC10d
   8012 = 50 80 GeV (LHC10h12)
   8013 = 80-inf GeV (LHC10h13)
   8102 = LHC11a2XX (where x = a - j)
  */
  Bool_t debugConfig =false;
  Bool_t ckbConfig = false;
  Bool_t productionConfig = false; // make this true for productions mode
  Int_t  iVersion = 1; // this is the train version for one datapass

  // bextra == 0 4 plus
  // bextra == 1 large pass1 split..
  // bextra == 2 3 plus
  Int_t bRun = 802; Int_t bExtra = 0;  char* cDate = "110908a";
  //  Int_t bRun = 8102; Int_t bExtra = 1;  char* cDate = "110725a";
  iAODanalysis = 0; 
  // 1 == Read Jets and tracks form the input AOD
  // needs the jet branchnames set explicitly
  // 2 == Read tracks from input and jets from output
  
  // 1 =  write the Full AOD for all events 
  // 2 =  write the Full AOD for triggered events
  // 3 =  write the deltaAOD for all events
  //  iFilterAnalysis = 2;kJetTriggerPtCut = 40; 
  iFilterAnalysis = 3;
  // iFilterAnalysis = 2;
  
  if (kPluginMode.Contains("merge")){
    // currently merging this one...
    //       cDate = "110717a";
    //    bRun = 802; Int_t bExtra = 0; cDate = "110812a";
  }
  kUseDebug = kFALSE;
  // this is for testing just one run...
  //  kGridMaxRunsFromList = 1;kUseDebug = kTRUE;

  kUseSysInfo = 100;
  
  kFillAOD                = kTRUE; // fill the aod
  kGridMergeExclude = ""; kSaveAOD = (1<<0)|(1<<1)|(1<<3); // 
  if(!productionConfig){
    kUsePAR              = kTRUE; // 
    kUseCPAR            = kTRUE;
    kGridFilesPerJob       = 20;
  }

  iPhysicsSelection = 1;

  if(iAODanalysis){
    // iAODAn
    if(iAODanalysis == 1){
      iJETAN = 0;
      iPWG4Cluster = 0;
      iJETSUBTRACT = 0;
    }
  }

  iJETAN = 3;
  iDIJETAN = 1; // only run on one JetFinder for the moment
  iPWG1QASym = 0; // excluded since also on QA train         
  iPWG4TmpSourceSara = 0; 
  iPWG4JetServices = 1; 
  iPWG4Fragmentation = 0; 
  iPWG4JetSpectrum = 7; 
  iPWG4UE = 0; // tmp off awating updates
  iPWG4LeadingUE = 1; 
  iPWG4CorrectionsUE = 0; // 19.07. OFF awaiting changes by Sara
  iPWG4PtQAMC       = 1;
  iPWG4PtSpectra    = 1;
  iPWG4PtQATPC      = 3;
  iPWG4PtTrackQA    = 1;
  iPWG4Cosmics      = 0; 
  iPWG4JetChem      = 1;
  iPWG4QGSep  = 1;
  iPWG4Minijet  = 1;
  iPWG4ThreeJets    = 0; // tmp off mem leak
  iPWG4KMeans       = 1; // Off no feedback 
  iPWG4Cluster      = 5;
  iPWG4Tagged       = 1; 
  iPWG4PartCorr     = 1;
  iPWG4CaloQA       = 1;
  iPWG4JetCorr      = 0; 
  iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
  iPWG4GammaConv    = 1; // TMP OFF for merging
  iPWG4CaloConv    = 0;  // 25.08. off: Output size 450MB in memory 03.09 crashes
  iPWG4omega3pi     = 1; //
  kDeltaAODJetName   = ""; // for OTF we do not need to write a delta/saves some extension gymnastics in the tasks       
  kDeltaAODPartCorrName = "";
  kPluginExecutableCommand = "cat wn.xml; root -b -q "; // dump the file list to stdout for debugging
  kPluginAliRootVersion    = ""; 
  kGridRunsPerMaster = 1; // To have the separate configs for eacj run in separate dirs

  if(bRun<100){ // private MC
    iPWG1QASym = 0;
    iPWG4TmpSourceSara = 0;
    iPWG4JetChem = 0;
    iPWG4UE = 0;
    iPWG4Cluster      = 0;
    iPWG4PtQAMC       = 0;
    iPWG4PtSpectra    = 0;
    iPWG4PtQATPC      = 0;
    iPWG4Cosmics      = 0; // tmp on
    iPWG4ThreeJets    = 0; 
    iPWG4KMeans       = 0;
    iPWG4PartCorr     = 0;
    iPWG4CaloQA       = 0; 
    iPWG4CaloConv     = 0; 
    iPWG4Tagged    = 0; 
    iPWG4JetCorr      = 0;     
    iPWG4GammaConv    = 0; 
    iPWG4omega3pi      = 0;     
    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE;  
    kUseTR                 = kTRUE; 
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    //   kErrorIgnoreLevel = 3001;
    //    kPluginFastReadOption = kTRUE;
  }

  if(bRun>=1000){
    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE; 
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    kUseTR                 = kTRUE; 
    iPhysicsSelection      = 1;
  }


  if (bRun == 1){
    kGridRunRange[0]       =  0;  
    kGridRunRange[1]       =  -1; 
    kGridDatadir           = "/alice/sim/PDC_08b/LHC09a1/"; 
    kGridDataSet           = "LHC09a1"; 
    kGridRunPattern        = "%05d"; 
    kGridLocalRunList      = "fp_lhc09a1.txt";
    kTrainName             = Form("pwg4train_LHC09a1_%s",cDate);
    iPhysicsSelection      = 0;
    kHighPtFilterMask      = 32;     
    kGridFilesPerJob       = 200;
  }
  else if (bRun == 2){
    kGridRunRange[0]       =  0;  
    kGridRunRange[1]       =  -1; 
    kGridDatadir           = "/alice/sim/PDC_08b/LHC09a2/"; 
    kGridDataSet           = "LHC09a2"; 
    kGridRunPattern        = "%05d"; 
    kGridLocalRunList      = "fp_lhc09a2.txt";
    kTrainName             = Form("pwg4train_LHC09a2_%s",cDate);
    iPhysicsSelection      = 0;
    kHighPtFilterMask      = 32;     
    kGridFilesPerJob       = 80;
  }
  else if (bRun == 3){
    kGridRunRange[0]       =  0;  
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/PDC_08b/LHC09a3/"; 
    kGridDataSet           = "LHC09a3"; 
    kGridRunPattern        = "%05d"; 
    kGridLocalRunList      = "fp_lhc09a3.txt";
    kTrainName             = Form("pwg4train_LHC09a3_%s",cDate);
    iPhysicsSelection      = 0;
    kHighPtFilterMask      = 32;     
    kGridFilesPerJob       = 80;
  }
  else if (bRun == 100){
    kGridRunRange[0]       =  0; 
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2009/LHC09d/"; 
    kGridDataSet           = "LHC09d"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    //    kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_runlist_pass1.txt";    kTrainName  = "pwg4train_LHC09d_pass1_100104";
    kGridPassPattern       = "pass6";    kGridLocalRunList      = "fp_runlist_pass6.txt";    kTrainName  = Form("pwg4train_LHC09d_pass6_%s",cDate);
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
  }
  else if (bRun == 200){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10b/"; 
    kGridDataSet           = "LHC10b"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_lhc10b_runlist_pass1.txt";    kTrainName  = Form("pwg4train_LHC10b_pass1_%s",cDate);
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 202){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10b/"; 
    kGridDataSet           = "LHC10b"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass2";    kGridLocalRunList      = "fp_lhc10b_runlist_pass2.txt";    kTrainName  = Form("pwg4train_LHC10b_pass2_%s",cDate);

    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks

    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback

    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 300){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10c/"; 
    kGridDataSet           = "LHC10c"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_lhc10c_runlist_900_pass1.txt";    kTrainName  = Form("pwg4train_LHC10c_900_pass1_%s",cDate);
    
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 302){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10c/"; 
    kGridDataSet           = "LHC10c"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass2";    kGridLocalRunList      = "fp_lhc10c_runlist_900_pass2.txt";    kTrainName  = Form("pwg4train_LHC10c_900_pass2_%s",cDate);

    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback

    
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 400){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10c/"; 
    kGridDataSet           = "LHC10c"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_lhc10c_runlist_7000_pass1.txt";    kTrainName  = Form("pwg4train_LHC10c_7000_pass1_%s",cDate);

    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 402){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridFilesPerJob       = 10;
    kGridDatadir           = "/alice/data/2010/LHC10c/"; 
    kGridDataSet           = "LHC10c"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass2";    kGridLocalRunList      = "fp_lhc10c_runlist_7000_pass2.txt";    kTrainName  = Form("pwg4train_LHC10c_7000_pass2_%s",cDate);

    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback

    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 500){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10d/"; 
    kGridDataSet           = "LHC10d"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_lhc10d_runlist_7000_pass1.txt";    kTrainName  = Form("pwg4train_LHC10d_7000_pass1_%s",cDate);

    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 502){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10d/"; 
    kGridDataSet           = "LHC10d"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 

    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback


    // swtich off most tasks for Pb+Pb
    iDIJETAN = 0; // only run on one JetFinder for the moment
    iPWG4Fragmentation = 0; 
    iPWG4LeadingUE = 0; 
    iPWG4JetChem      = 0;
    iPWG4QGSep  = 0;
    iPWG4Minijet  = 0;
    iPWG4PartCorr     = 0;
    iPWG4GammaConv    = 0; 

    // running as light a possible 

    iPWG4PtQAMC     = 1;
    iPWG4PtSpectra   = 1;
    iPWG4PtQATPC   = 1;
    iPWG4JetSpectrum = 1;
    iPWG4JetServices  = 1; // !!!!!!!!!!! 
    iPWG4Cluster      = 1;// not 5....
    kHighPtFilterMask = 1<<8; // 256 TPC related to SPD     

    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass2";    kGridLocalRunList      = "fp_lhc10d_runlist_7000_pass2.txt";    kTrainName  = Form("p4t_10d_7000_p2_%s",cDate);
    //    kDeltaAODJetName   = "AliAOD.Jets.root"; 
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 801){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10h/"; 
    kGridDataSet           = "LHC10h"; 
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 

    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback

    // swtich off most tasks for Pb+Pb
    iDIJETAN = 0; // only run on one JetFinder for the moment
    iPWG4Fragmentation = 0; 
    iPWG4LeadingUE = 0; 
    iPWG4JetChem      = 0;
    iPWG4QGSep  = 0;
    iPWG4Minijet  = 0;
    iPWG4PartCorr     = 0;
    iPWG4GammaConv    = 0; 

    // running as light a possible 

    iPWG4PtQAMC     = 1;
    iPWG4PtSpectra   = 1;
    iPWG4PtQATPC   = 1;
    iPWG4JetSpectrum = 1;
    iPWG4JetServices  = 1; // !!!!!!!!!!! 
    iPWG4Cluster      = 1;// not 5....
    kHighPtFilterMask = 1<<8; // 256 TPC related to SPD     


    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass1_4plus";    kGridLocalRunList      = "fp_lhc10h_pass1_4plus.txt";    
    kTrainName  = Form("t_10h_1_4p_%s",cDate);
    if(bExtra==1){
      kGridPassPattern       = "pass1";    kGridLocalRunList      = "fp_lhc10h_pass1.txt";   
      //      kTrainName  = Form("train_pass1_%s",cDate);
      kTrainName  = Form("t_10h_p1_%s",cDate);
    }
    else if(bExtra==2){
      kGridPassPattern       = "pass1_plusplusplus";    kGridLocalRunList      = "fp_lhc10h_pass1_3plus.txt";
      //      kTrainName  = Form("train_LHC10h_pass1_3p_%s",cDate);      
      kTrainName  = Form("t_10h_p1_3p_%s",cDate);
    }
    kSaveAOD            = 1;  
    kDeltaAODJetName   = "AliAOD.Jets.root";kSaveAOD = 2;

    if (kPluginMode.Contains("merge")){
      kSaveAOD = 0; // 
    }
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    //    gROOT->LoadMacro("cleanXML.C");
    // cleanXML();
  }
  else if (bRun == 802){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/data/2010/LHC10h/"; 
    kGridDataSet           = "LHC10h"; 
    kGridExtraAliendirLevel  = Form("/%s",cDate);kGridOutdir = "output";
    kGridRunPattern        = "%09d"; 
    kUseKinefilter         = kFALSE;
    kIsMC                  = kFALSE;  // is real! 
    kUseMC                 = kFALSE; 
    kUseAODMC              = kFALSE; 

    if(iVersion==1){
    // switch off tasks with no feedback...
      iPWG4ThreeJets    = 0; // tmp off mem leak
      iPWG4KMeans       = 0;  // OFF no FEEDBACK
      iPWG4Tagged       = 0; // OFF crashes on MC 
      iPWG4CaloQA       = 0; // OFF not needed on MC   
      iPWG4JetCorr      = 0; 
      iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
      iPWG4omega3pi     = 0; // OFF no feedback
      
      // swtich off most tasks for Pb+Pb
      iDIJETAN = 0; // only run on one JetFinder for the moment
      iPWG4LeadingUE = 0; 
      iPWG4JetChem      = 0;
      iPWG4QGSep  = 0;
      iPWG4Minijet  = 0;
      iPWG4PartCorr     = 0;
      iPWG4GammaConv    = 0; 
      
      // running as light a possible 
      
      iPWG4PtQAMC     = 0;  // 
      iPWG4PtQATPC   = 0;  // 
      iPWG4PtSpectra   = 0;  //  
      iPWG4PtTrackQA    = 0; // currently not needed 08.09.2011
      iPWG4JetSpectrum = 1; 
      iPWG4JetServices  = 1; // !!!!!!!!!!! 
      iPWG4Cluster      = 1;// not 5....
      kHighPtFilterMask = 1<<4|1<<8; // Global tracks with SPD requirment global constraitn for the rest
      iPWG4Fragmentation = 1;
    //
    }// version1


    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "pass2";    // reset for AODs below

    kGridLocalRunList      = "fp_lhc10h_pass2c.txt";    
    if(bExtra==1)kGridLocalRunList      = "fp_lhc10h_pass2.txt";    
    kTrainName  = Form("p4t_10h_pass2_%s",cDate);
    if (kPluginMode.Contains("merge")){
      //      kGridLocalRunList      = "fp_lhc10h_pass2b.txt";    
      kGridLocalRunList      = "out_merge_100_1.txt";    
    }
    if(iAODanalysis)kTrainName  = Form("tAOD_10h_pass2_%s",cDate);

    if(iFilterAnalysis==1){
      kSaveAOD            = 1;
      kGridFilesPerJob       = 5;
      kDeltaAODJetName   = "";
      kFilterAOD = false;
    }
    else if(iFilterAnalysis==2){
      kSaveAOD            = 1;
      kGridFilesPerJob    = 20;
      kDeltaAODJetName   = "";
      kFilterAOD = true;
    }
    else if(iFilterAnalysis == 3){
      kSaveAOD            = 2;
      kGridFilesPerJob    = 20;
      kDeltaAODJetName   = "AliAOD.Jets.root";
      kFilterAOD = true;
    }

    if(iAODanalysis==0){
      // select what is running in the spectrum task, list of jet names is filles automatically
      kGridFilesPerJob       = 30;
      /*
      // Thnsparses... 

      ############# Possible jet branches ###################
      1: jetsAOD_UA104_B0_Filter00016_Cut01000
      2: jetsAOD_UA104_B0_Filter00016_Cut02000
      3: jetsAOD_UA104_B2_Filter00016_Cut01000
      4: jetsAOD_UA104_B2_Filter00016_Cut02000
      5: clustersAOD_KT04_B0_Filter00016_Cut00150_Skip00
      6: clustersAOD_KT04_B0_Filter00016_Cut00150_Skip00RandomConeSkip00
      7: clustersAOD_KT04_B0_Filter00016_Cut00150_Skip00RandomCone_random
      8: clustersAOD_KT04_B0_Filter00016_Cut02000_Skip00
      9: clustersAOD_KT04_B0_Filter00016_Cut02000_Skip00RandomConeSkip00
      10: clustersAOD_KT04_B0_Filter00016_Cut02000_Skip00RandomCone_random
      11: clustersAOD_KT02_B0_Filter00016_Cut00150_Skip00
      12: clustersAOD_ANTIKT04_B0_Filter00016_Cut00150_Skip02
      13: clustersAOD_ANTIKT04_B0_Filter00016_Cut00150_Skip02RandomConeSkip02
      14: clustersAOD_ANTIKT04_B0_Filter00016_Cut00150_Skip02RandomCone_random
      15: clustersAOD_ANTIKT04_B0_Filter00016_Cut02000_Skip02
      16: clustersAOD_ANTIKT02_B0_Filter00016_Cut00150_Skip00
      17: clustersAOD_KT04_B1_Filter00016_Cut00150_Skip00RandomConeSkip00
      18: clustersAOD_KT04_B1_Filter00016_Cut00150_Skip00RandomCone_random
      19: clustersAOD_ANTIKT04_B1_Filter00016_Cut00150_Skip02
      20: clustersAOD_ANTIKT04_B1_Filter00016_Cut00150_Skip02RandomConeSkip02
      21: clustersAOD_ANTIKT04_B1_Filter00016_Cut00150_Skip02RandomCone_random
      22: clustersAOD_ANTIKT02_B1_Filter00016_Cut00150_Skip00
      23: clustersAOD_KT04_B2_Filter00016_Cut00150_Skip00RandomConeSkip00
      24: clustersAOD_KT04_B2_Filter00016_Cut00150_Skip00RandomCone_random
      25: clustersAOD_ANTIKT04_B2_Filter00016_Cut00150_Skip02
      26: clustersAOD_ANTIKT04_B2_Filter00016_Cut00150_Skip02RandomConeSkip02
      27: clustersAOD_ANTIKT04_B2_Filter00016_Cut00150_Skip02RandomCone_random
      28: clustersAOD_ANTIKT02_B2_Filter00016_Cut00150_Skip00
      29: clustersAOD_KT04_B2_Filter00016_Cut02000_Skip00RandomConeSkip00
      30: clustersAOD_KT04_B2_Filter00016_Cut02000_Skip00RandomCone_random
      31: clustersAOD_ANTIKT04_B2_Filter00016_Cut02000_Skip02
            */
      /*
	############# Possible jet branches ###################
	1: jetsAOD_UA104_B0_Filter00272_Cut01000
	2: jetsAOD_UA104_B0_Filter00272_Cut02000
	3: jetsAOD_UA104_B2_Filter00272_Cut01000
	4: jetsAOD_UA104_B2_Filter00272_Cut02000
	5: clustersAOD_KT04_B0_Filter00272_Cut00150_Skip00
	6: clustersAOD_KT04_B0_Filter00272_Cut00150_Skip00RandomConeSkip00
	7: clustersAOD_KT04_B0_Filter00272_Cut00150_Skip00RandomCone_random
	8: clustersAOD_KT04_B0_Filter00272_Cut01000_Skip00
	9: clustersAOD_KT04_B0_Filter00272_Cut01000_Skip00RandomConeSkip00
	10: clustersAOD_KT04_B0_Filter00272_Cut01000_Skip00RandomCone_random
	11: clustersAOD_KT04_B0_Filter00272_Cut02000_Skip00
	12: clustersAOD_KT04_B0_Filter00272_Cut02000_Skip00RandomConeSkip00
	13: clustersAOD_KT04_B0_Filter00272_Cut02000_Skip00RandomCone_random
	14: clustersAOD_KT02_B0_Filter00272_Cut00150_Skip00
	15: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02
	16: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02RandomConeSkip02
	17: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02RandomCone_random
	18: clustersAOD_ANTIKT04_B0_Filter00272_Cut01000_Skip02
	19: clustersAOD_ANTIKT04_B0_Filter00272_Cut02000_Skip02
	20: clustersAOD_ANTIKT02_B0_Filter00272_Cut00150_Skip00
	21: clustersAOD_KT04_B1_Filter00272_Cut00150_Skip00RandomConeSkip00
	22: clustersAOD_KT04_B1_Filter00272_Cut00150_Skip00RandomCone_random
	23: clustersAOD_ANTIKT04_B1_Filter00272_Cut00150_Skip02
	24: clustersAOD_ANTIKT04_B1_Filter00272_Cut00150_Skip02RandomConeSkip02
	25: clustersAOD_ANTIKT04_B1_Filter00272_Cut00150_Skip02RandomCone_random
	26: clustersAOD_ANTIKT02_B1_Filter00272_Cut00150_Skip00
	27: clustersAOD_KT04_B2_Filter00272_Cut00150_Skip00RandomConeSkip00
	28: clustersAOD_KT04_B2_Filter00272_Cut00150_Skip00RandomCone_random
	29: clustersAOD_ANTIKT04_B2_Filter00272_Cut00150_Skip02
	30: clustersAOD_ANTIKT04_B2_Filter00272_Cut00150_Skip02RandomConeSkip02
	31: clustersAOD_ANTIKT04_B2_Filter00272_Cut00150_Skip02RandomCone_random
	32: clustersAOD_ANTIKT02_B2_Filter00272_Cut00150_Skip00
	33: clustersAOD_KT04_B2_Filter00272_Cut01000_Skip00RandomConeSkip00
	34: clustersAOD_KT04_B2_Filter00272_Cut01000_Skip00RandomCone_random
	35: clustersAOD_ANTIKT04_B2_Filter00272_Cut01000_Skip02
	36: clustersAOD_KT04_B2_Filter00272_Cut02000_Skip00RandomConeSkip00
	37: clustersAOD_KT04_B2_Filter00272_Cut02000_Skip00RandomCone_random
	38: clustersAOD_ANTIKT04_B2_Filter00272_Cut02000_Skip02
      */
      
      // in the first map we fill the correlations we want to plot
      // in the jet back map we associated the branche used for background calculation
      // to fetch the multiplicity

      // UA1 
      kJetMapSpectrum.Add(4,2);
      kJetBackMapSpectrum.Add(4,8);

      // anti kT 150 MeV
      kJetMapSpectrum.Add(29,15);
      kJetBackMapSpectrum.Add(29,5);
      kJetBackMapSpectrum.Add(15,5);

      // anti kT 1000 MeV
      kJetMapSpectrum.Add(35,18);
      kJetBackMapSpectrum.Add(35,8);      
      kJetBackMapSpectrum.Add(35,8);

      // anti kT 2000 MeV
      kJetMapSpectrum.Add(38,19);
      kJetBackMapSpectrum.Add(35,11);      
      kJetBackMapSpectrum.Add(19,11);


      // anti kT 0.2
      kJetMapSpectrum.Add(32,14);
      kJetBackMapSpectrum.Add(32,5);      
      kJetBackMapSpectrum.Add(14,5);

      // random cones
      kJetMapSpectrum.Add(6,7);
      kJetBackMapSpectrum.Add(6,5);
      kJetBackMapSpectrum.Add(7,5);

      kJetMapSpectrum.Add(9,10);
      kJetBackMapSpectrum.Add(9,8);
      kJetBackMapSpectrum.Add(10,8);
      
    }
    else if (iAODanalysis){
      kGridFilesPerJob       = 20;
      if(iAODanalysis==1){
	kGridPassPattern = "pass2/AOD049";
	iJETAN = 0;
	iPWG4Cluster = 0;
	iJETSUBTRACT = 0;


	/* 
	   reading from AOD043 input
	   1 jetsAOD_UA104_B0_Filter00128_Cut01000 1 0x102c068d0
	   2 jetsAOD_UA104_B2_Filter00128_Cut01000 1 0x102c068d0
	   3 jetsAOD_SISCONE04_B0_Filter00128_Cut00150 1 0x102c068d0
	   4 clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00 1 0x102c068d0
	   5 clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00_random 1 0x102c068d0
	   6 clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomConeSkip00 1 0x102c068d0
	   7 clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomCone_random 1 0x102c068d0
	   8 clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02 1 0x102c068d0
	   9 clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02_random 1 0x102c068d0
	   10 clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomConeSkip02 1 0x102c068d0
	   11 clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomCone_random 1 0x102c068d0
	   12 clustersAOD_ANTIKT02_B0_Filter00128_Cut00150_Skip00 1 0x102c068d0
	   13 clustersAOD_ANTIKT02_B0_Filter00128_Cut00150_Skip00_random 1 0x102c068d0
	   14 jetsAOD_SISCONE04_B1_Filter00128_Cut00150 1 0x102c068d0
	   15 clustersAOD_ANTIKT04_B1_Filter00128_Cut00150_Skip02 1 0x102c068d0
	   16 clustersAOD_ANTIKT02_B1_Filter00128_Cut00150_Skip00 1 0x102c068d0
	*/

	kJetListSpectrum.Add(new TObjString("jetsAOD_UA104_B0_Filter00128_Cut01000"));
	kJetListSpectrum.Add(new TObjString("jetsAOD_UA104_B2_Filter00128_Cut01000"));
	kJetListSpectrum.Add(new TObjString("jetsAOD_SISCONE04_B0_Filter00128_Cut00150"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00_random"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomConeSkip00"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomCone_random"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02_random"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomConeSkip02"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomCone_random"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT02_B0_Filter00128_Cut00150_Skip00"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT02_B0_Filter00128_Cut00150_Skip00_random"));
	kJetListSpectrum.Add(new TObjString("jetsAOD_SISCONE04_B1_Filter00128_Cut00150"));
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT04_B1_Filter00128_Cut00150_Skip02"));	
	kJetListSpectrum.Add(new TObjString("clustersAOD_ANTIKT02_B1_Filter00128_Cut00150_Skip00"));

	// anti kT 150 MeV
	kJetMapSpectrum.Add(15,8);
	kJetBackMapSpectrum.Add(8,4);
	kJetBackMapSpectrum.Add(15,4);

	// anti kT R  = 0.2
	kJetMapSpectrum.Add(17,12);
	kJetBackMapSpectrum.Add(17,4);      
	kJetBackMapSpectrum.Add(12,4);
	
	// random cones
	kJetMapSpectrum.Add(6,7);
	kJetBackMapSpectrum.Add(6,4);
	kJetBackMapSpectrum.Add(7,4);
	
	kJetMapSpectrum.Add(10,11);
	kJetBackMapSpectrum.Add(10,4);
	kJetBackMapSpectrum.Add(11,4);
	kDeltaAODJetNameInput = "AliAOD.Jets.root";
	kDeltaAODJetName = "";
	kHighPtFilterMask      = 128; // centrally produced AOD     
      }
      else if(iAODanalysis==2){
	/*
	  ############# Possible jet branches ###################
	  1: jetsAOD_UA104_B0_Filter00128_Cut01000
	  2: jetsAOD_UA104_B0_Filter00128_Cut02000
	  3: jetsAOD_UA104_B2_Filter00128_Cut01000
	  4: jetsAOD_UA104_B2_Filter00128_Cut02000
	  5: clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00
	  6: clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomConeSkip00
	  7: clustersAOD_KT04_B0_Filter00128_Cut00150_Skip00RandomCone_random
	  8: clustersAOD_KT04_B0_Filter00128_Cut02000_Skip00
	  9: clustersAOD_KT04_B0_Filter00128_Cut02000_Skip00RandomConeSkip00
	  10: clustersAOD_KT04_B0_Filter00128_Cut02000_Skip00RandomCone_random
	  11: clustersAOD_KT02_B0_Filter00128_Cut00150_Skip00
	  12: clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02
	  13: clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomConeSkip02
	  14: clustersAOD_ANTIKT04_B0_Filter00128_Cut00150_Skip02RandomCone_random
	  15: clustersAOD_ANTIKT04_B0_Filter00128_Cut02000_Skip02
	  16: clustersAOD_ANTIKT02_B0_Filter00128_Cut00150_Skip00
	  17: clustersAOD_KT04_B1_Filter00128_Cut00150_Skip00Skip00RandomConeSkip00
	  18: clustersAOD_KT04_B1_Filter00128_Cut00150_Skip00Skip00RandomCone_random
	  19: clustersAOD_ANTIKT04_B1_Filter00128_Cut00150_Skip02
	  20: clustersAOD_ANTIKT04_B1_Filter00128_Cut00150_Skip02Skip02RandomConeSkip02
	  21: clustersAOD_ANTIKT04_B1_Filter00128_Cut00150_Skip02Skip02RandomCone_random
	  22: clustersAOD_ANTIKT02_B1_Filter00128_Cut00150_Skip00
	  23: clustersAOD_KT04_B2_Filter00128_Cut00150_Skip00Skip00RandomConeSkip00
	  24: clustersAOD_KT04_B2_Filter00128_Cut00150_Skip00Skip00RandomCone_random
	  25: clustersAOD_ANTIKT04_B2_Filter00128_Cut00150_Skip02
	  26: clustersAOD_ANTIKT04_B2_Filter00128_Cut00150_Skip02Skip02RandomConeSkip02
	  27: clustersAOD_ANTIKT04_B2_Filter00128_Cut00150_Skip02Skip02RandomCone_random
	  28: clustersAOD_ANTIKT02_B2_Filter00128_Cut00150_Skip00
	  29: clustersAOD_KT04_B2_Filter00128_Cut02000_Skip00Skip00RandomConeSkip00
	  30: clustersAOD_KT04_B2_Filter00128_Cut02000_Skip00Skip00RandomCone_random
	  31: clustersAOD_ANTIKT04_B2_Filter00128_Cut02000_Skip02
	 */
	
	// anti kT 150 MeV
	kJetMapSpectrum.Add(25,12);
	kJetBackMapSpectrum.Add(25,5);
	kJetBackMapSpectrum.Add(12,5);

	// anti kT R  = 0.2
	kJetMapSpectrum.Add(28,16);
	kJetBackMapSpectrum.Add(28,5);      
	kJetBackMapSpectrum.Add(16,5);

	// anti kT 2 GeV
	kJetMapSpectrum.Add(31,15);
	kJetBackMapSpectrum.Add(31,8);      
	kJetBackMapSpectrum.Add(15,8);

	
	// random cones
	kJetMapSpectrum.Add(9,10);
	kJetBackMapSpectrum.Add(10,5);
	kJetBackMapSpectrum.Add(9,5);
	
	kJetMapSpectrum.Add(13,14);
	kJetBackMapSpectrum.Add(14,5);
	kJetBackMapSpectrum.Add(13,5);
	kDeltaAODJetNameInput = "";


	kDeltaAODJetName   = "AliAOD.Jets.root";kSaveAOD = 2;
	kHighPtFilterMask      = 128; // centrally produced AOD     
	iPWG4Fragmentation = 0; // off for a pass 

	
      }
    }
    if (kPluginMode.Contains("merge")){
      kSaveAOD = 0; // 
    }
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    //    gROOT->LoadMacro("cleanXML.C");
    // cleanXML();
  }
  else if (bRun == 1012){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10a12/"; 
    kGridDataSet           = "LHC10a12"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10a12.txt";    kTrainName  = Form("pwg4train_LHC10a12_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 250; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 1014){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10a14/"; 
    kGridDataSet           = "LHC10a14"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10a14.txt";    kTrainName  = Form("pwg4train_LHC10a14_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 250; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 2001){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10b1/"; 
    kGridDataSet           = "LHC10b1"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10b1.txt";    kTrainName  = Form("pwg4train_LHC10b1_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 250; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 2002){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10b2/"; 
    kGridDataSet           = "LHC10b2"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10b2.txt";    kTrainName  = Form("pwg4train_LHC10b2_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 250; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 4001){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10d1/"; 
    kGridDataSet           = "LHC10d1"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10d1.txt";    kTrainName  = Form("pwg4train_LHC10d1_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 50; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 4002){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10d2/"; 
    kGridDataSet           = "LHC10d2"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10d2.txt";    kTrainName  = Form("pwg4train_LHC10d2_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 100; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 4004){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10d4/"; 
    kGridDataSet           = "LHC10d4"; 
    kGridRunPattern        = "%06d"; 
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10d4.txt";    kTrainName  = Form("pwg4train_LHC10d4_%s",cDate);
    kGridRunsPerMaster     = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 100; // only few events in a sim file
    // stich of tasks not from PWG4JetTasks
    iPWG4TmpSourceSara = 0; 
    iPWG4UE = 0; //
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4CaloConv    = 0;  // 25.08. off: Output size 03.09 crashes 
    iPWG4omega3pi     = 0; // OFF no feedback
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if(bRun==5014){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10e14/120823/"; 
    kGridRunPattern        = "%d"; 
    kGridDataSet           = "LHC10e14"; 
    kGridLocalRunList      = "fp_lhc10e14.txt";   kTrainName  = Form("pwg4train_LHC10e14_%s",cDate);
    kGridPassPattern       = "";

    iPWG4TmpSourceSara = 0; 
    iPWG4UE = 0; // off not needed on Jet+Jet
    iPWG4LeadingUE = 0; // off not needed on Jet +Jet
    iPWG4CorrectionsUE = 0; // off not needed on Jet +Jet
    iPWG4Cosmics      = 0; // off not needed on Jet +Jet MC
    iPWG4JetChem      = 0; // OFF no FEEDBACK
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4GammaConv    = 0; // TMP OFF cuts not updated not so important for jet+jet
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4CaloConv    = 0;  // 25.08. off: Output size 03.09 crashes 
    iPWG4omega3pi     = 0; // OFF no feedback
    iPWG1QASym = 0; // excluded since also on QA train         

    kDeltaAODJetName   = ""; // for OTF we do not need to write a delta/saves some extension gymnastics in the tasks       

    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE;  
    kUseTR                 = kTRUE; 
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    //    gSystem->Exec("cp wn_lhc10b.xml wn.xml");
    kGridFilesPerJob       = 50; // only few events in a sim file
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 6006){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10f6/"; 
    kGridDataSet           = "LHC10f6"; 
    kGridRunPattern        = "%06d"; 
    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE;  
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    kGridFilesPerJob       = 100;
    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback
    iPWG4TmpSourceSara = 0; 
    iPWG4UE = 0; //
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4CaloConv    = 0;  // 25.08. off: Output size 03.09 crashes 
    iPWG4PartCorr    = 0;  // OFF cjecked back with Gustavo
    iPWG4omega3pi     = 0; // OFF no feedback
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10f6.txt";    kTrainName  = Form("pwg4train_LHC10f6_%s",cDate);
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    //    gROOT->LoadMacro("cleanXML.C");
    //    cleanXML();
  }
  else if (bRun == 8012){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 
    kGridDatadir           = "/alice/sim/LHC10h12/"; 
    kGridDataSet           = "LHC10h12"; 
    kGridRunPattern        = "%06d"; 
    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE;  
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    kGridFilesPerJob       = 50;
    // switch off tasks with no feedback...


    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback


    // swtich off most tasks for Pb+Pb
    iDIJETAN = 0; // only run on one JetFinder for the moment
    iPWG4Fragmentation = 0; 
    iPWG4LeadingUE = 0; 
    iPWG4JetChem      = 0;
    iPWG4QGSep  = 0;
    iPWG4Minijet  = 0;
    iPWG4PartCorr     = 0;
    iPWG4GammaConv    = 0; 

    // running as light a possible 

    iPWG4PtQAMC     = 0;
    iPWG4PtSpectra   = 0;
    iPWG4PtQATPC   = 0;
    iPWG4JetSpectrum = 0;
    iPWG4JetServices  = 0; // !!!!!!!!!!! 
    iPWG4Cluster      = 1;// not 5....
    kHighPtFilterMask = 256;     

    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_lhc10h12.txt";    kTrainName  = Form("pwg4train_LHC10h12_%s",cDate);
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }
  else if (bRun == 8102){
    kGridRunRange[0]       =  0;  // 0 is a bad number :(, needs a manual commit in alien...
    kGridRunRange[1]       =  -1; // 


    char a;
    if(bExtra<10)a = bExtra+97;

    kGridDatadir           = Form("/alice/sim/LHC11a2%c/",a); 
    kGridDataSet           = Form("LHC11a2%c/",a); 
    kGridRunPattern        = "%06d"; 
    kUseKinefilter         = kTRUE;
    kIsMC                  = kTRUE;  
    kUseMC                 = kTRUE; 
    kUseAODMC              = kTRUE; 
    kGridFilesPerJob       = 50;
    // switch off tasks with no feedback...


    // switch off tasks with no feedback...
    iPWG4ThreeJets    = 0; // tmp off mem leak
    iPWG4KMeans       = 0;  // OFF no FEEDBACK
    iPWG4Tagged       = 0; // OFF crashes on MC 
    iPWG4CaloQA       = 0; // OFF not needed on MC   
    iPWG4JetCorr      = 0; 
    iPWG4JCORRAN      = 0;  // TMP OFF: Writes a Tree
    iPWG4omega3pi     = 0; // OFF no feedback


    // swtich off most tasks for Pb+Pb
    iDIJETAN = 0; // only run on one JetFinder for the moment
    iPWG4Fragmentation = 0; 
    iPWG4LeadingUE = 0; 
    iPWG4JetChem      = 0;
    iPWG4QGSep  = 0;
    iPWG4Minijet  = 0;
    iPWG4PartCorr     = 0;
    iPWG4GammaConv    = 0; 

    // running as light a possible 
    iJETSUBTRACT = 0; // no subtraction


    iPWG4PtQAMC     = 1;
    iPWG4PtSpectra   = 1;
    iPWG4PtQATPC   = 1;
    iPWG4JetSpectrum = 1;
    iPWG4JetServices  = 1; // !!!!!!!!!!! 
    iPWG4Cluster      = 1;// not 5....
    kHighPtFilterMask = 1<<4|1<<8;     



    if(iFilterAnalysis==1){
      kSaveAOD            = 1;
      kDeltaAODJetName   = "";
      kFilterAOD = false;
    }
    else if(iFilterAnalysis==2){
      kJetTriggerPtCut = 20; //pt 
      kSaveAOD            = 1;
      kDeltaAODJetName   = "";
      kFilterAOD = true;
    }
    else if(iFilterAnalysis == 3){
      kSaveAOD            = 2;
      kDeltaAODJetName   = "AliAOD.Jets.root";
      kFilterAOD = true;
    }
    kGridFilesPerJob       = 100;
    /*
############# Possible jet branches ###################
  1: jetsAOD_UA104_B0_Filter00272_Cut01000
      2: jetsAOD_UA104_B0_Filter00272_Cut02000
      3: jetsAODMC_UA104_B0_Filter00272_Cut01000
      4: jetsAODMC2_UA104_B0_Filter00272_Cut01000
      5: clustersAOD_KT06_B0_Filter00272_Cut00150_Skip00
      6: clustersAOD_KT04_B0_Filter00272_Cut00150_Skip00
      7: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02
      8: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02RandomConeSkip02
      9: clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02RandomCone_random
      10: clustersAOD_ANTIKT04_B0_Filter00272_Cut02000_Skip02
      11: clustersAOD_ANTIKT02_B0_Filter00272_Cut00150_Skip00
      12: clustersAODMC_KT06_B0_Filter00272_Cut00150_Skip00
      13: clustersAODMC2_KT06_B0_Filter00272_Cut00150_Skip00
      14: clustersAODMC_ANTIKT04_B0_Filter00272_Cut00150_Skip00
      15: clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00
*/
    // CLEAN XML FILES LOCALLY AND ON ALIEN WHEN STARTING A NEW PASS!
    kGridPassPattern       = "";    kGridLocalRunList      = "fp_runlist_lhc11a2.txt";    kTrainName  = Form("pwg4train_LHC11a2%c_%s",a,cDate);
    kGridRunsPerMaster = 1; // Physcicsselection does not support more than on run per job
    kGridFilesPerJob       = 100;
    // stich of tasks not from PWG4JetTasks
    gROOT->LoadMacro("cleanXML.C");
    cleanXML();
  }



  if(ckbConfig){
    // switch off all but my analyses
    iPWG4KMeans     = 0;
    iPWG1QASym = 0;
    iPWG4TmpSourceSara = 0;
    iPWG4UE = 0;
    iPWG4PtQAMC        = 0;
    iPWG4PtSpectra     = 0;
    iPWG4PtQATPC        = 0;
    iPWG4ThreeJets        = 0;
    iPWG4KMeans     = 0;
    iPWG4Tagged    = 0; 
    iPWG4PartCorr     = 0;
    iPWG4CaloQA     = 0;
    iPWG4CaloConv     = 0; 
    iPWG4JetCorr      = 0;     
    iPWG4GammaConv    = 0; 
    iPWG4JetChem = 0; // tmp on
    iPWG4omega3pi      = 0;     
    kDeltaAODJetName   = ""; // for OTF we do not need to write a delta/saves some extension gymnastics in the tasks       
    kDeltaAODPartCorrName = "";
    kUseDebug = kFALSE;
    kPluginAliRootVersion    = ""; 
    kGridFilesPerJob       = 60;
    kTrainName             = Form("%s_ckb",kTrainName.Data());
  }

  if(debugConfig){
    // debug mode 
    //    kUsePAR                = kFALSE; // cannot patch the macro for local test, need to laod FASTjet libs before loading FASTJETA.so
    //    kUseCPAR               = kFALSE;
    kUseSysInfo = 100;
    kUseDebug = kTRUE;
    kGridLocalRunList      = "fp_runlist_pass4_debug.txt";
    kTrainName             = Form("pwg4train_LHC09d_debug",cDate);
    
    // kPluginExecutableCommand = "root -b -q ";
    kPluginExecutableCommand = "cat wn.xml; echo \"Root.Stacktrace:         yes\" > .rootrc; root -b -q ";
    //    kPluginExecutableCommand = "echo \"run -b -q  pwg4train_LHC09d_debug.C\" > gdb.cmd; echo where >> gdb.cmd; echo quit >> gdb.cmd; echo y >> gdb.cmd; gdb -x gdb.cmd -batch root.exe";
  }

  if (kPluginMode.Contains("test")){
    kJetTriggerPtCut = 0.01; 
    if(kAnalysisMode.Contains("grid")){
      //      kPluginExecutableCommand = "source ~/setup_root.txt; alienroot -b -q";      
      kPluginExecutableCommand = "root -b -q";      
      //      kPluginExecutableCommand = "valgrind --tool=memcheck --error-limit=no --max-stackframe=3060888 --suppressions=$ROOTSYS/etc/valgrind-root.supp --leak-check=full --num-callers=15 --log-file=valgrind_memcheck.log root.exe -b -q";
      kPluginExecutableCommand = "export ALICE_ROOT=./ROOTFILES/;" + kPluginExecutableCommand; 
      kUseSysInfo = 1;
      kUseDebug = kTRUE;
      if(bRun==802){
	kGridLocalRunList      = "fp_lhc10h_anchor.txt";
      }
      kTrainName             = Form("pwg4train_test");
      if(iAODanalysis)kNumberOfEvents     = 500;
    }
    else{
      // local
      if(iAODanalysis)kNumberOfEvents     = 2000;
      kUseSysInfo = 1;
	kUseDebug = kTRUE;
	kTrainName             = Form("pwg4train_test_local");
	kLocalDataList = "local_esd_lhc10d_pass2.txt";
	kUsePAR              = kFALSE; // cannot patch the macro for local test, need to laod FASTjet libs before loading FASTJETA.so
	kUseCPAR            = kFALSE;
       	if(bRun==802){
	  kLocalDataList = "local_esd_lhc10h.txt";
	  if(iAODanalysis)	  kLocalDataList = "local_aod_lhc10h.txt";
	}
	//	iPWG4PtTrackQA    = 0;
	//	iPWG4PtQAMC       = 0;
	//	iPWG4PtSpectra    = 0;
	//	iPWG4PtQATPC      = 0;
	//	iPWG4PtTrackQA    = 0;
	//	iPWG4Cluster      = 0;
	kUseCPAR            = kFALSE;
	kUsePAR            = kFALSE;
	//	kNumberOfEvents     = 200;
	// all OFF
    }
  }
  if(kPluginAliRootVersion.Length()==0){
    kPluginExecutableCommand = "export ALICE_ROOT=./ROOTFILES/;" + kPluginExecutableCommand; 
  } 

  if (kPluginMode.Contains("merge")){
    // swtich of task where macros changed in the meantime
  }


}
