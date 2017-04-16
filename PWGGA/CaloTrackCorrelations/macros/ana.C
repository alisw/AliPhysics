/// \file ana.C
/// \ingroup CaloTrackCorrMacros
/// \brief Example of execution macro
///
/// Example macro to do analysis with the
/// analysis classes in CaloTrackCorrelations,
/// in local, grid or plugin modes.
///
/// Pay attention to the options and definitions
/// set in the lines below
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//---------------------------------------------------------------------------
/// Different analysis modes
enum anaModes
{
    mLocal  = 0, /// Analyze locally files in your computer.
    mPROOF  = 1, /// Analyze files on GRID with Plugin
    mPlugin = 2, /// Analyze files on GRID with Plugin
    mGRID   = 3  /// Analyze files on GRID, jobs launched from aliensh
};

//---------------------------------------------------------------------------
// Settings to read locally several files, only for "mLocal" mode
// The different values are default, they can be set with environmental 
// variables: INDIR, PATTERN, NFILES, respectivelly

char * kInDir   = "/user/data/files/";  /// Global,  path to data files
char * kPattern = ""; /// Data are in files kInDir/kPattern+i
Int_t  kFile    = 2;  /// Number of files to analyze in local mode.

//---------------------------------------------------------------------------
// Dataset for proof analysis, mode=mPROOF
// char * kDataset = "/alice/vernet/PbPb_LHC10h_ESD";

char *  kDatasetPROOF     = "/alice/vernet/LHC11b_149646";
Int_t   kDatasetNMaxFiles = 20;
TString ccin2p3UserName   = "arbor" ;
TString alienUserName     = "narbor" ;

//---------------------------------------------------------------------------
// Collection file for grid analysis

char * kXML = "collection.xml"; /// Global name for the xml collection file with data on grid

//---------------------------------------------------------------------------
// Scale histograms from file. Change to kTRUE when xsection file exists
// Put name of file containing xsection
// Put number of events per ESD file
// This is an specific case for normalization of Pythia files.
const char * kXSFileName = "pyxsec.root"; /// Name of file with pT-hard cross sections

// Container of xs if xs in file pyxsec_hist.root
TArrayF* xsArr;
TArrayI* trArr;

//---------------------------------------------------------------------------

// Set some default values, but used values are set in the code!

Bool_t  kMC        = kFALSE; /// With real data kMC = kFALSE
TString kInputData = "ESD"; /// ESD, AOD, MC, deltaAOD
Int_t   kYear      = 2011;
TString kCollision = "pp";
Bool_t  outAOD     = kFALSE; /// Some tasks doesnt need it.
TString kTreeName;
TString kPass      = "";
char    kTrigger[1024];
Int_t   kRun       = 0;

//___________________________
/// Main execution method. It:
/// * 1) loads the needed libraries in method LoadLibraries
/// * 2) depending on the files path, run etc, the variables year, collision type, data type, are obtained in methods CheckInputData and CheckEnvironmentVariables
/// * 3) put the data files in a list to be passed to the analysis frame in method CreateChain
/// * 4) In case of MC pt-Hard bin simulations, the file containing the cross sections is read and scaling parameter is obtained via the method GetAverageXsection
/// * 5) The analysis frame is initialized via de analysis manager
/// * 6) Different general analysis are initialized: Physics selection, centrality etc.
/// * 7) Specialized analysis are initialized: AliAnalysistaskCounter, AliAnalysisTaskEMCALClusterizer, AliAnalysisTaskCaloTrackCorrelations and executed for different settings.
/// * 8) The output/input containers are passed to the analysis manager
/// * 9) The analysis is executed
///
/// \param mode: analysis mode defined in enum anaModes
//___________________________
void ana(Int_t mode=mGRID)
{
  //--------------------------------------------------------------------
  // Load analysis libraries

  LoadLibraries(mode) ;
  //gSystem->ListLibraries();
  
  //-----------------------------------------------------------------------------
  // Create chain from ESD and from cross sections files, look below for options.
  
  // Set kInputData and kTreeName looking to the kINDIR
  
  CheckInputData(mode);
  
  // Check global analysis settings  
  
  CheckEnvironmentVariables();
  
  printf("*********************************************\n");
  printf("*** Input data < %s >, pass %s, tree < %s >, MC?  < %d > ***\n",kInputData.Data(),kPass.Data(),kTreeName.Data(),kMC);
  printf("*********************************************\n");
  
  TChain * chain   = new TChain(kTreeName) ;
  TChain * chainxs = new TChain("Xsection") ;
  CreateChain(mode, chain, chainxs); 
  
  Double_t scale  = -1;
  printf("===== kMC %d, chainxs %p\n",kMC,chainxs);
  
  if(kMC)
  {
    //Get the cross section
    Double_t xsection = 0;
    Float_t  ntrials  = 0;
    Int_t    nfiles =  0;
    
    Bool_t ok = GetAverageXsection(chainxs, xsection, ntrials, nfiles);
    
    printf("n xs files %d",nfiles);
    
    if(ok)
    {
      Int_t  nEventsPerFile = chain->GetEntries() / nfiles;
      
      Double_t trials = ntrials / nEventsPerFile ;
      
      scale = xsection / trials;
      
      printf("Get Cross section : nfiles  %d, nevents %d, nevents per file %d \n",nfiles, chain->GetEntries(),nEventsPerFile);
      printf("                    ntrials %d, trials %2.2f, xs %2.2e, scale factor %2.2e\n", ntrials,trials,xsection,scale);
      
      if(chainxs->GetEntries()!=chain->GetEntries()) printf("CAREFUL: Number of files in data chain %d, in cross section chain %d \n",
                                                            chainxs->GetEntries(),chain->GetEntries());
    } // ok
    
    // comment out this line in case the simulation did not have the cross section files produced in the directory
    if( scale <= 0  || !ok)
    { printf( "STOP, cross section not available! nfiles %d \n", chainxs->GetEntries() ) ; return ; }
    
  }

  printf("*********************************************\n");
  printf("number of entries # %lld, skipped %d\n", chain->GetEntries()) ; 	
  printf("*********************************************\n");
  
  if(!chain)
  { 
    printf("STOP, no chain available\n"); 
    return;
  }
  
  AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
  
  //------------------------------------------
  //  Alien handler part
  //------------------------------------------
  AliAnalysisGrid *alienHandler=0x0;
  if(mode==mPlugin)
  {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    alienHandler = CreateAlienHandler();
    if (!alienHandler) return;
  }  
  
  //--------------------------------------
  // Make the analysis manager
  //-------------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  //AliAnalysisManager::SetUseProgressBar(kTRUE);
  //mgr->SetSkipTerminate(kTRUE);
  //mgr->SetNSysInfo(1);
  
  if(mode==mPlugin)
  {
    // Connect plugin to the analysis manager
    mgr->SetGridHandler(alienHandler);
  }
  
  // MC handler
  if((kMC || kInputData == "MC") && !kInputData.Contains("AOD"))
  {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
    if( kInputData == "MC") 
    {
      cout<<"MC INPUT EVENT HANDLER"<<endl;
      mgr->SetInputEventHandler(NULL);
    }
  }
  
  
  // AOD output handler
  if(kInputData!="deltaAOD" && outAOD)
  {
    cout<<"Init output handler"<<endl;
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("aod.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
  }
  
  //input
  
  if(kInputData == "ESD")
  {
    // ESD handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdHandler);
    cout<<"ESD handler "<<mgr->GetInputEventHandler()<<endl;
  }
  else if(kInputData.Contains("AOD"))
  {
    // AOD handler
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
    if(kInputData == "deltaAOD") aodHandler->AddFriend("deltaAODCaloTrackCorr.root");
    cout<<"AOD handler "<<mgr->GetInputEventHandler()<<endl;
  }
  
  //mgr->RegisterExternalFile("deltaAODCaloTrackCorr.root");
  //mgr->SetDebugLevel(1); // For debugging, do not uncomment if you want no messages.
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
  
  //-------------------------------------------------------------------------
  // Define task, put here any other task that you want to use.
  //-------------------------------------------------------------------------
  
  // Physics selection
  if(kInputData=="ESD" && !kMC)
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kMC); 
  }
  
  // Centrality, valid for Run1, but superseeded by new task below
  if(kCollision=="PbPb" && kInputData=="ESD" && kRun < 200000)
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  }

  {
    // New centrality/multiplicity selector
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask * task = AddTaskMultSelection(kFALSE); // user mode:
    
    //use the default calibration for runs which have not yet been calibrated
    task->SetUseDefaultCalib(kTRUE); // data
    task->SetUseDefaultMCCalib(kTRUE); // MC
  }
  
  if(kCollision=="PbPb")
  {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
    AliVZEROEPSelectionTask  * EPV0 = AddTaskVZEROEPSelection();  
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
    AliEPSelectionTask * EP = AddTaskEventplane();
  }
  
  
  
  // Simple event counting tasks
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/AddTaskCounter.C");   

  AliAnalysisTask* count    = AddTaskCounter("",kMC);   // All, fill histo with cross section and trials if kMC is true
  AliAnalysisTask* countmb  = AddTaskCounter("MB"); // Min Bias
  AliAnalysisTask* countany = AddTaskCounter("Any"); 
  AliAnalysisTask* countint = AddTaskCounter("AnyINT");// Min Bias
  
  if(!kMC)
  {
    AliAnalysisTaskCounter* countemg = AddTaskCounter("EMCEGA"); 
    AliAnalysisTaskCounter* countemj = AddTaskCounter("EMCEJE"); 
    if(kCollision=="PbPb")
    {
      AliAnalysisTaskCounter* countcen = AddTaskCounter("Central"); 
      AliAnalysisTaskCounter* countsce = AddTaskCounter("SemiCentral"); 
      AliAnalysisTaskCounter* countssce= AddTaskCounter("SemiOrCentral"); 
      AliAnalysisTaskCounter* countphP = AddTaskCounter("PHOSPb"); 
    }
    else
    {
      AliAnalysisTaskCounter* countem1 = AddTaskCounter("EMC1"); // Trig Th > 1.5 GeV approx
      AliAnalysisTaskCounter* countem7 = AddTaskCounter("EMC7"); // Trig Th > 4-5 GeV 
      AliAnalysisTaskCounter* countphp = AddTaskCounter("PHOS"); 
    }
  }  
  // -----------------
  // Photon conversion
  // ----------------- 
/*  
  if(kInputData=="ESD"){
    printf("* Configure photon conversion analysis in macro \n");
    TString arguments = "-run-on-train -use-own-xyz  -force-aod -mc-off ";
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConversion/macros/ConfigGammaConversion.C");
    AliAnalysisTaskGammaConversion * taskGammaConversion = 
    ConfigGammaConversion(arguments,mgr->GetCommonInputContainer());
    taskGammaConversion->SelectCollisionCandidates();
    
    // Gamma Conversion AOD to AODPWG4Particle
    AliAnalysisTaskGCPartToPWG4Part * taskGCToPC = new AliAnalysisTaskGCPartToPWG4Part("GCPartToPWG4Part");
    taskGCToPC->SetGammaCutId("90035620401003321022000000090");
    mgr->AddTask(taskGCToPC);
    mgr->ConnectInput  (taskGCToPC, 0, mgr->GetCommonInputContainer() );
    mgr->ConnectOutput (taskGCToPC, 0, mgr->GetCommonOutputContainer()); 
  }
*/  
  
  Bool_t kPrint   = kFALSE;
  Bool_t deltaAOD = kFALSE;
  gROOT->LoadMacro("AddTaskCaloTrackCorr.C");   // $ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/AddTaskEMCALClusterize.C"); 
  
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");  
  //AliAnalysisTaskCaloTrackCorrelation * qatask = AddTaskCalorimeterQA(kInputData,kYear,kPrint,kMC); 
  
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/AddTaskEMCALTriggerQA.C");  
  //AliAnalysisTaskEMCALTriggerQA * qatrigtask = AddTaskEMCALTriggerQA(); 
  
  // Calibration, bad map ...
  
  Bool_t calibEE = kTRUE; // It is set automatically, but here we force to use it or not in any case
  Bool_t calibTT = kTRUE; // It is set automatically, but here we force to use it or not in any case
  if(kRun < 122195 || (kRun > 126437 && kRun < 136851) || kMC) calibTT=kFALSE ; // Recalibration parameters not available for LHC10a,b,c,e,f,g
  Bool_t badMap  = kTRUE; // It is set automatically, but here we force to use it or not in any case  
  
  if(kCollision=="pp")
  {
    printf("==================================== \n");
    printf("CONFIGURE ANALYSIS FOR PP COLLISIONS \n");
    printf("==================================== \n");
    
    Bool_t  mixing    = kTRUE; // mixing in reader and hadron correlation, and pi0
    Bool_t  clTM      = kTRUE;
    Bool_t  reTM      = kFALSE; // Recalculate matches if not already done in clusterizer
    Bool_t  anTM      = kTRUE;  // Remove matched
    Bool_t  exo       = kTRUE;  // Remove exotic cells
    Bool_t  clnonlin  = kTRUE;  // Apply non linearity (clusterization)
    Bool_t  annonlin  = kFALSE; // Apply non linearity (analysis)
    Int_t   minEcell  = 50;     // 50  MeV (10 MeV used in reconstruction)
    Int_t   minEseed  = 100;    // 100 MeV
    Int_t   dTime     = 0;      // default, 250 ns
    Int_t   wTime     = 30;     // default 425 < T < 825 ns, careful if time calibration is on
    Int_t   unfMinE   = 15;     // Remove cells with less than 15 MeV from cluster after unfolding
    Int_t   unfFrac   = 1;      // Remove cells with less than 1% of cluster energy after unfolding

    //Trigger
    TString clTrigger   = "";   
    TString anTrigger   = "EMC7";  
    
    if(kMC) 
    {
      clTrigger   = "";
      anTrigger   = "";
      dTime       = 0;      
      wTime       = 0;    
    }
    
    Bool_t  selectEvents = kFALSE; // Select events depending on V0, pile-up and vertex quality
    Bool_t  qa     = kTRUE; // Do besides calorimeter QA analysis
    Bool_t  hadron = kTRUE; // Do besides charged track correlations analysis    
    
    //Analysis with clusterizer V1
    
    TString arrayNameV1 = "";
    
    AliAnalysisTaskEMCALClusterize * clv1 = AddTaskEMCALClusterize(arrayNameV1,outAOD,kMC,exo,"V1",clTrigger, clTM,
                                                                   minEcell,minEseed,dTime,wTime,unfMinE,unfFrac,
                                                                   calibEE,badMap,calibTT,clnonlin);    
    
    printf("Name of clusterizer1 array: %s\n",arrayNameV1.Data());
    
    if(!kMC)
    {
      
      AliAnalysisTaskCaloTrackCorrelation *anav1trig   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                              kYear,kCollision,anTrigger,arrayNameV1,mixing,reTM,anTM,
                                                                              -1,-1, qa, hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    }
    
    AliAnalysisTaskCaloTrackCorrelation *anav1mb     = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,"AnyINT",arrayNameV1,mixing,reTM,anTM,
                                                                            -1,-1, qa, hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    
    
    
    //Analysis with clusterizer V2
    TString arrayNameV2 = "";
    AliAnalysisTaskEMCALClusterize * clv2 = AddTaskEMCALClusterize(arrayNameV2,outAOD,kMC,exo,"V2",clTrigger, clTM,
                                                                   minEcell,minEseed,dTime,wTime,
                                                                   calibEE,badMap,calibTT,clnonlin);    
    
    printf("Name of clusterizer2 array: %s\n",arrayNameV2.Data());
    
    hadron = kFALSE;
    if(!kMC)
    {
      
      
      AliAnalysisTaskCaloTrackCorrelation *anav2tr   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,anTrigger,arrayNameV2,mixing,reTM,anTM, 
                                                                            -1,-1,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    }
    
    
    AliAnalysisTaskCaloTrackCorrelation *anav2mb     = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,"AnyINT",arrayNameV2,mixing,reTM,anTM,
                                                                            -1,-1, qa, hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
  }
  
  if(kCollision=="PbPb")
  {
    printf("====================================== \n");
    printf("CONFIGURE ANALYSIS FOR PbPb COLLISIONS \n");
    printf("====================================== \n");
    Bool_t  mixing    = kTRUE; 
    Bool_t  clTM      = kTRUE;
    Bool_t  reTM      = kFALSE; // Recalculate matches if not already done in clusterizer
    Bool_t  anTM      = kTRUE;  // Remove matched
    Bool_t  exo       = kTRUE;  // Remove exotic cells
    Bool_t  clnonlin  = kTRUE;  // Apply non linearity (clusterization)
    Bool_t  annonlin  = kFALSE; // Apply non linearity (analysis)
    Int_t   minEcell  = 150;    // 50  MeV (10 MeV used in reconstruction)
    Int_t   minEseed  = 300;    // 100 MeV
    Int_t   dTime     = 0;      // default, 250 ns
    Int_t   wTime     = 0;      // default 425 < T < 825 ns
    Int_t   unfMinE   = 15;     // Remove cells with less than 15 MeV from cluster after unfolding
    Int_t   unfFrac   = 1;      // Remove cells with less than 1% of cluster energy after unfolding
    
    // Trigger
    TString clTrigger = ""; 
    TString anTrigger = "EMCGA";  
    if(kMC) 
    {
      clTrigger = "";
      anTrigger = "";
      dTime       = 0;      
      wTime       = 0;    
    }
    
    Bool_t  selectEvents = kFALSE; // Select events depending on V0, pile-up and vertex quality
    Bool_t  qa     = kTRUE; // Do besides calorimeter QA analysis
    Bool_t  hadron = kTRUE; // Do besides charged track correlations analysis    
    
    //Analysis with clusterizer V1
    
    TString arrayNameV1 = "";
    AliAnalysisTaskEMCALClusterize * clv1 = AddTaskEMCALClusterize(arrayNameV1,outAOD,kMC,exo,"V1",clTrigger, clTM,
                                                                   minEcell,minEseed,dTime,wTime,unfMinE,unfFrac,
                                                                   calibEE,badMap,calibTT,clnonlin);    
    
    printf("Name of clusterizer1 array: %s\n",arrayNameV1.Data());
    
    if(!kMC)
    {
      AliAnalysisTaskCaloTrackCorrelation *anav1c   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,anTrigger,arrayNameV1,mixing,reTM,anTM, 
                                                                           0,20,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
      AliAnalysisTaskCaloTrackCorrelation *anav1m   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,anTrigger,arrayNameV1,mixing,reTM,anTM,
                                                                           20,40,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
      AliAnalysisTaskCaloTrackCorrelation *anav1p   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC,  selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,anTrigger,arrayNameV1,mixing,reTM,anTM,
                                                                           60,80,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    }
    
    //Analysis with clusterizer V2
    
    TString arrayNameV2 = "";
    AliAnalysisTaskEMCALClusterize * clv2 = AddTaskEMCALClusterize(arrayNameV2,outAOD,kMC,exo,"V2",clTrigger, clTM,
                                                                   minEcell,minEseed,dTime,wTime,unfMinE,unfFrac,
                                                                   calibEE,badMap,calibTT,clnonlin);        
    
    printf("Name of clusterizer2 array: %s\n",arrayNameV2.Data());
    
    hadron = kFALSE;
    
    if(!kMC)
    {
      
      AliAnalysisTaskCaloTrackCorrelation *anav2cT   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,anTrigger,arrayNameV2,mixing,reTM,anTM, 
                                                                            0,20,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
      AliAnalysisTaskCaloTrackCorrelation *anav2mT   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,anTrigger,arrayNameV2,mixing,reTM,anTM,
                                                                            20,40,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
      AliAnalysisTaskCaloTrackCorrelation *anav2pT   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                            kYear,kCollision,anTrigger,arrayNameV2,mixing,reTM,anTM,
                                                                            60,80,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    }
    
    AliAnalysisTaskCaloTrackCorrelation *anav2cMB   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,"AnyINT",arrayNameV2,mixing,reTM,anTM, 
                                                                           0,20,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    AliAnalysisTaskCaloTrackCorrelation *anav2mMB   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,"AnyINT",arrayNameV2,mixing,reTM,anTM,
                                                                           20,40,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    AliAnalysisTaskCaloTrackCorrelation *anav2pMB   = AddTaskCaloTrackCorr(kInputData, "EMCAL", kMC, selectEvents, exo, annonlin, outputFile.Data(), 
                                                                           kYear,kCollision,"AnyINT",arrayNameV2,mixing,reTM,anTM,
                                                                           60,80,qa,hadron,calibEE,badMap,calibTT,deltaAOD,kPrint,scale,kRun);
    
    
  }
  
  //-----------------------
  // Run the analysis
  //-----------------------    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  if      (mode == mPlugin) mgr->StartAnalysis("grid");
  else if (mode == mPROOF ) mgr->StartAnalysis("proof",chain);
  else                      mgr->StartAnalysis("local",chain);
  
  cout <<" Analysis ended sucessfully "<< endl ;
}

//_____________________________
/// Load analysis libraries.
//_____________________________
void  LoadLibraries(Int_t mode)
{
  if (mode == mPROOF)
  {
    //TProof::Mgr("ccalpmaster")->SetROOTVersion("ALICE_v5-27-06b");
    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/EnableAliRootForLAF.C");
    TProof* proof = EnableAliRootForLAF("ccaplmaster",nPROOFWorkers.Data(),ccin2p3UserName.Data(),alienUserName.Data(),"",kFALSE,kTRUE,kTRUE,"OADB:ANALYSIS:ANALYSISalice:AOD:ESD:CORRFW:STEERBase:EMCALUtils:PHOSUtils:PWGCaloTrackCorrBase:PWGGACaloTrackCorrelations:PWGPPEMCAL");
    
    //  TProof* proof = TProof::Open("ccaplmaster",Form("workers=%s",nPROOFWorkers.Data()));
    
    //     //proof->ClearPackages();
    //     proof->UploadPackage("STEERBase");
    //     proof->UploadPackage("ESD");
    //     proof->UploadPackage("AOD");
    //     proof->UploadPackage("ANALYSIS");
    //     proof->UploadPackage("OADB");
    //     proof->UploadPackage("ANALYSISalice");
    //     proof->UploadPackage("CORRFW");
    //     //proof->UploadPackage("JETAN");
    //     proof->UploadPackage("PHOSUtils");
    //     proof->UploadPackage("EMCALUtils");
    //     proof->UploadPackage("PWGCaloTrackCorrBase");
    //     proof->UploadPackage("PWGGACaloTrackCorrelations");
    //     proof->UploadPackage("PWGPPEMCAL");
    
    //     proof->EnablePackage("STEERBase");
    //     proof->EnablePackage("ESD");
    //     proof->EnablePackage("AOD");
    //     proof->EnablePackage("ANALYSIS");
    //     proof->EnablePackage("OADB");
    //     proof->EnablePackage("ANALYSISalice");
    //     proof->EnablePackage("CORRFW");
    //     //proof->EnablePackage("JETAN");
    //     proof->EnablePackage("PHOSUtils");
    //     proof->EnablePackage("EMCALUtils");
    //     proof->EnablePackage("PWGCaloTrackCorrBase");
    //     proof->EnablePackage("PWGGACaloTrackCorrelations");
    //     proof->EnablePackage("PWGPPEMCAL");
    return;
  }  
  
  //--------------------------------------
  // Load the needed libraries via par files if modified
  //--------------------------------------
  //SetupPar("EMCALUtils");
  //SetupPar("EMCALraw");
  //SetupPar("EMCALbase");
  //SetupPar("EMCALsim");
  //SetupPar("EMCALrec");
  
  //SetupPar("PWGPPEMCAL");
  
  //SetupPar("PWGCaloTrackCorrBase");
  //SetupPar("PWGGACaloTrackCorrelations"); 

  //SetupPar("PWGGAGammaConv"); 
  
  // needed for plugin?
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS");
  gSystem->AddIncludePath("-I./");
}

//_________________________________
/// Load par files, create analysis libraries
/// For testing, if par file already decompressed and modified
/// classes then do not decompress.
//_________________________________
void SetupPar(char* pararchivename)
{
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
    
  TString parpar(Form("%s.par", pararchivename)) ; 
  
  if ( gSystem->AccessPathName(pararchivename) )
  {
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
  {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh"))
    {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
  {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

//______________________________________
/// Sets input data and tree strings.
//______________________________________
void CheckInputData(const anaModes mode)
{
  TString ocwd = gSystem->WorkingDirectory();
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if(mode == mLocal){    
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files 
    //to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    if(gSystem->Getenv("INDIR"))  
      kInDir = gSystem->Getenv("INDIR") ; 
    else cout<<"INDIR not set, use default: "<<kInDir<<endl;	
    
    TString sindir(kInDir);
    if     (sindir.Contains("pass1")) kPass = "pass1";
    else if(sindir.Contains("pass2")) kPass = "pass2";
    else if(sindir.Contains("pass3")) kPass = "pass3";
    
    if(gSystem->Getenv("PATTERN"))   
      kPattern = gSystem->Getenv("PATTERN") ; 
    else  cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    cout<<"INDIR   : "<<kInDir<<endl;
    cout<<"NFILES  : "<<kFile<<endl;
    
    char fileE[120] ;   
    char fileA[120] ;   
    char fileG[120] ;
    char fileEm[120] ;   
    for (Int_t event = 0 ; event < kFile ; event++) {
      sprintf(fileE,  "%s/%s%d/AliESDs.root",    kInDir,kPattern,event) ; 
      sprintf(fileA,  "%s/%s%d/AliAOD.root",     kInDir,kPattern,event) ; 
      sprintf(fileG,  "%s/%s%d/galice.root",     kInDir,kPattern,event) ; 
      sprintf(fileEm, "%s/%s%d/embededAOD.root", kInDir,kPattern,event) ; 
      
      TFile * fESD = TFile::Open(fileE) ; 
      TFile * fAOD = TFile::Open(fileA) ; 
      
      //Check if file exists and add it, if not skip it
      if (fESD) 
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        TFile * fG = TFile::Open(fileG);
        if(fG) { kMC = kTRUE; fG->Close();}
        else     kMC = kFALSE;
        
        // Get run number
        TTree* esdTree = (TTree*)fESD->Get("esdTree");
        AliESDEvent* esd = new AliESDEvent();
        esd->ReadFromTree(esdTree);
        esdTree->GetEvent(0);
        kRun = esd->GetRunNumber();
        
        return;
      }
      else if(fAOD)
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        
        // Get run number
        TTree* aodTree = (TTree*)fAOD->Get("aodTree");
        AliAODEvent* aod = new AliAODEvent();
        aod->ReadFromTree(aodTree);
        aodTree->GetEvent(0);
        kRun = aod->GetRunNumber();
        return;
      }
      else if(TFile::Open(fileEm))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC        = kTRUE;
        
        return;
      }
      else if(TFile::Open(fileG))
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC        = kTRUE;
        return;
      }
    }
    
    if(fESD) fESD->Close();
    if(fAOD) fAOD->Close();
    
  }// local files analysis
  
  //------------------------------
  //GRID xml files
  //-----------------------------
  else if(mode == mGRID){
    //Get colection file. It is specified by the environmental
    //variable XML
    
    if(gSystem->Getenv("XML") )
      kXML = gSystem->Getenv("XML");
    else
      sprintf(kXML, "collection.xml") ; 
    
    if (!TFile::Open(kXML)) {
      printf("No collection file with name -- %s -- was found\n",kXML);
      return ;
    }
    else cout<<"XML file "<<kXML<<endl;
    
    //Load necessary libraries and connect to the GRID
    gSystem->Load("libNetx") ; 
    gSystem->Load("libRAliEn"); 
    TGrid::Connect("alien://") ;
    
    //Feed Grid with collection file
    TGridCollection * collection = (TGridCollection*) TAlienCollection::Open(kXML);
    if (! collection) {
      AliError(Form("%s not found", kXML)) ; 
      return kFALSE ; 
    }
    TGridResult* result = collection->GetGridResult("",0 ,0);
    
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      
      if     (alienURL.Contains("pass1")) kPass = "pass1";
      else if(alienURL.Contains("pass2")) kPass = "pass2";
      else if(alienURL.Contains("pass3")) kPass = "pass3";
      
      kRun = AliAnalysisManager::GetRunFromAlienPath(alienURL.Data());
      printf("Run number from alien path = %d\n",kRun);
      
      TFile * fAOD = 0 ; 
      //Check if file exists and add it, if not skip it
      if (alienURL.Contains("AliESDs.root"))  
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        alienURL.ReplaceAll("AliESDs.root","galice.root");
        if(TFile::Open(alienURL)) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if(alienURL.Contains("AliAOD.root"))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        fAOD = TFile::Open(alienURL);
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if(alienURL.Contains("embededAOD.root"))
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC=kTRUE;
        return;
      }
      else if(alienURL.Contains("galice.root"))
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC=kTRUE;
        return;
      } 
    }
  }// xml analysis
  //------------------------------
  //PROOF files
  //-----------------------------
  else if(mode == mPROOF){
    
    TFileCollection* coll  = gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
    
    TIter iter(coll->GetList());
    
    TFileInfo* fileInfo = 0;
    while ((fileInfo = dynamic_cast<TFileInfo*> (iter())))
    {
      if (fileInfo->GetFirstUrl()) {
        TString ProofURL = fileInfo->GetFirstUrl()->GetUrl();
        cout << "================== " << ProofURL << endl ; 
        
        if     (ProofURL.Contains("pass1")) kPass = "pass1";
        else if(ProofURL.Contains("pass2")) kPass = "pass2";
        else if(ProofURL.Contains("pass3")) kPass = "pass3";
        
        kRun = AliAnalysisManager::GetRunFromAlienPath(ProofURL.Data());
        printf("Run number from alien path = %d\n",kRun);
        
        TFile * fAOD = 0 ; 
        //Check if file exists and add it, if not skip it
        if (ProofURL.Contains("AliESDs.root"))  
        {
          kTreeName  = "esdTree";
          kInputData = "ESD";
          alienURL.ReplaceAll("AliESDs.root","galice.root");
          if(TFile::Open(ProofURL)) kMC=kTRUE;
          else kMC = kFALSE;
          
          return;
        }
        else if(ProofURL.Contains("AliAOD.root"))
        {
          kTreeName  = "aodTree";
          kInputData = "AOD";
          fAOD = TFile::Open(ProofURL);
          if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
          else kMC = kFALSE;
          return;
        }
        else if(ProofURL.Contains("embededAOD.root"))
        {
          kTreeName  = "aodTree";
          kInputData = "AOD";
          kMC=kTRUE;
          return;
        }
        else if(ProofURL.Contains("galice.root"))
        {
          kTreeName  = "TE";
          kInputData = "MC";
          kMC=kTRUE;
          return;
        } 
      }
    }
  }// proof analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}

//_____________________________________________________________________
/// Fills chain with data files paths.
//_____________________________________________________________________
void CreateChain(const anaModes mode, TChain * chain, TChain * chainxs)
{
  TString ocwd = gSystem->WorkingDirectory();
  
  if(kInputData == "AOD")
  {
    xsArr = new TArrayF;
    trArr = new TArrayI;
  }
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if(mode == mLocal)
  {
    //If you want to add several ESD files sitting in a common directory INDIR
    //Specify as environmental variables the directory (INDIR), the number of files
    //to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    if(gSystem->Getenv("INDIR"))
      kInDir = gSystem->Getenv("INDIR") ;
    else cout<<"INDIR not set, use default: "<<kInDir<<endl;
    
    if(gSystem->Getenv("PATTERN"))
      kPattern = gSystem->Getenv("PATTERN") ;
    else  cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    if(gSystem->Getenv("NFILES"))
      kFile = atoi(gSystem->Getenv("NFILES")) ;
    else cout<<"NFILES not set, use default: "<<kFile<<endl;
    
    //Check if env variables are set and are correct
    if ( kInDir  && kFile)
    {
      printf("Get %d files from directory %s\n",kFile,kInDir);
      if ( ! gSystem->cd(kInDir) )
      {//check if ESDs directory exist
        printf("%s does not exist\n", kInDir) ;
        return ;
      }
      
      //if(gSystem->Getenv("XSFILE"))
      //kXSFileName = gSystem->Getenv("XSFILE") ;
      //else cout<<" XS file name not set, use default: "<<kXSFileName<<endl;
      char * kGener = gSystem->Getenv("GENER");
      if(kGener)
      {
        cout<<"GENER "<<kGener<<endl;
        if     (!strcmp(kGener,"PYTHIA")) kXSFileName = "pyxsec.root";
        else if(!strcmp(kGener,"HERWIG")) kXSFileName = "hexsec.root";
        else cout<<" UNKNOWN GENER, use default: "<<kXSFileName<<endl;
      }
      else cout<<" GENER not set, use default xs file name: "<<kXSFileName<<endl;
      
      if(kInputData=="AOD")
      {
        kXSFileName = "pyxsec_hists.root";
        xsArr->Set(kFile);
        trArr->Set(kFile);
      }
      
      cout<<"INDIR   : "<<kInDir     <<endl;
      cout<<"NFILES  : "<<kFile      <<endl;
      cout<<"PATTERN : "<<kPattern   <<endl;
      cout<<"XSFILE  : "<<kXSFileName<<endl;
      
      TString datafile="";
      if     (kInputData == "ESD")        datafile = "AliESDs.root" ;
      else if(kInputData.Contains("AOD")) datafile = "AliAOD.root"  ;
      else if(kInputData == "MC")         datafile = "galice.root"  ;
      
      // Loop on ESD/AOD/MC files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ;
      char file[120] ;
      char filexs[120] ;
      
      for (event = 0 ; event < kFile ; event++) {
        sprintf(file,   "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ;
        sprintf(filexs, "%s/%s%d/%s", kInDir,kPattern,event,kXSFileName) ;
        TFile * fData = 0 ;
        // Check if file exists and add it, if not skip it
        if ( fData = TFile::Open(file))
        {
          if ( fData->Get(kTreeName) )
          {
            printf("++++ Adding %s\n", file) ;
            chain->AddFile(file);
            
            if(kInputData != "AOD")
            {
              chainxs->Add(filexs) ;
            }
            else
            {
              TFile*  fxsec = TFile::Open(filexs);
              if(fxsec)
              {
                TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
                if(!key)
                {
                  fxsec->Close();
                  printf("No key!");
                  continue;
                }
                
                TList *list = dynamic_cast<TList*>(key->ReadObj());
                if(!list)
                {
                  fxsec->Close();
                  printf("No list!");
                  continue;
                }
                
                Float_t xsection = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
                Int_t   ntrials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
                fxsec->Close();
                
                xsArr->SetAt(xsection,event);
                trArr->SetAt(ntrials,event);
                
                printf("recovered xs %f, ntrials %d, event %d\n",xsection,ntrials, event);
                //chainxs->Add(tree);
                //fileTMP->Close();
              } // fxsec exists
            } // xs in AODs
            
          }
        }
        else
        {
          printf("---- Skipping %s\n", file) ;
          skipped++ ;
        }
      }
    }
    else {
      TString input = "AliESDs.root" ;
      cout<<">>>>>> No list added, take a single file <<<<<<<<< "<<input<<endl;
      chain->AddFile(input);
    }
    
  }// local files analysis
  
  //------------------------------
  // GRID xml files
  //------------------------------
  else if(mode == mGRID)
  {
    //Get colection file. It is specified by the environmental
    //variable XML
    
    //Feed Grid with collection file
    TGridCollection * collection = (TGridCollection*) TAlienCollection::Open(kXML);
    if (! collection) {
      AliError(Form("%s not found", kXML)) ;
      return kFALSE ;
    }
    
    TGridResult* result = collection->GetGridResult("",0 ,0);
    
    // Makes the ESD chain
    printf("*** Getting the Chain       ***\n");
    for (Int_t index = 0; index < result->GetEntries(); index++) {
      TString alienURL = result->GetKey(index, "turl") ;
      cout << "================== " << alienURL << endl ;
      chain->Add(alienURL) ;
      
      if(kInputData != "AOD")
      {
        alienURL.ReplaceAll("AliESDs.root",kXSFileName);
        alienURL.ReplaceAll("AliAOD.root",kXSFileName);
        chainxs->Add(alienURL) ;
      }
      else
      {
        alienURL.ReplaceAll("AliESDs.root","pyxsec_hists.root");
        alienURL.ReplaceAll("AliAOD.root", "pyxsec_hists.root");
        TFile*  fxsec = TFile::Open(alienURL);
        if(fxsec)
        {
          TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
          if(!key)
          {
            fxsec->Close();
            printf("No key!");
            continue;
          }
          
          TList *list = dynamic_cast<TList*>(key->ReadObj());
          if(!list)
          {
            fxsec->Close();
            printf("No list!");
            continue;
          }
          
          Float_t xsection = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
          Int_t   ntrials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
          fxsec->Close();
          
          xsArr->SetAt(xsection,event);
          trArr->SetAt(ntrials,event);
          
          printf("recovered xs %f, ntrials %d, event %d\n",xsection,ntrials, event);
          
        } // fxsec exists
      } // xs in AODs
    }
  }// xml analysis
  
  //------------------------------
  // PROOF
  //------------------------------
  else if (mode == mPROOF)
  {
    TFileCollection* ds= gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
    
    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/dataset_management/CreateChainFromDataSet.C");
    chain = CreateChainFromDataSet(ds, kTreeName , kDatasetNMaxFiles);
    printf("chain has %d entries\n",chain->GetEntries());
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
}

//______________________________
/// Access one data file and set the year,
/// collision type and run number.
/// It is possible to set them via external parameters.
//______________________________
void CheckEnvironmentVariables()
{
  sprintf(kTrigger,"");
  
  Bool_t bRecalibrate = kFALSE;
  Bool_t bBadChannel = kFALSE;
  
  for (int i=0; i< gApplication->Argc();i++){
    
#ifdef VERBOSEARGS
    
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
    
#endif
    
    TString sRun = "";
    
    if (!(strcmp(gApplication->Argv(i),"--trigger")))
      sprintf(trigger,gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--recalibrate")))
      bRecalibrate = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--badchannel")))
      bBadChannel = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--run")))
    {
      sRun = gApplication->Argv(i+1);
      if(sRun.Contains("LHC10")) {
        kYear = 2010;
      }
      else {
        if(kRun <=0){
          kRun = atoi(gApplication->Argv(i+1));
        }
        else printf("** Run number already set  to %d, do not set to %d\n",kRun,atoi(gApplication->Argv(i+1)));
      }//numeric run
    }//--run available
    
  }// args loop
  
  if     ( kRun < 140000) 
  {
    kYear = 2010;
    if( kRun >= 136851 ) kCollision = "PbPb";
  }
  else if( kRun < 170600)
  {
    kYear = 2011;
    if( kRun >= 166500 ) kCollision = "PbPb";
  }
  else if( kRun < 200000 )
  {
    kYear = 2012;
    if( kRun >= 194000 ) kCollision = "pPb";
  }
  else if( kRun < 247000 )
  {
    kYear = 2015;
    if( kRun >= 244820 ) kCollision = "PbPb";
  }
  else
  {
    kYear = 2016;
  }
  
  if(kMC) sprintf(kTrigger,"");
  
  printf("*********************************************\n");
  //printf("*** Settings trigger %s, recalibrate %d, remove bad channels %d, year %d, collision %s, run %d ***\n",
  //       kTrigger,bRecalibrate,bBadChannel, kYear,kCollision.Data(), kRun);
  printf("*** Settings year %d, collision %s, run %d ***\n",kYear,kCollision.Data(), kRun);
  printf("*********************************************\n");
  
}

//______________________________________________________________________________
/// Read the PYTHIA statistics from the file pyxsec.root created by
/// the function WriteXsection():
/// integrated cross section (xsection) and
/// the  number of Pyevent() calls (ntrials)
/// and calculate the weight per one event xsection/ntrials
/// The spectrum calculated by a user should be
/// multiplied by this weight, something like this:
/// TH1F *userSpectrum ... // book and fill the spectrum
/// userSpectrum->Scale(weight)
///
/// Yuri Kharlov 19 June 2007
/// Gustavo Conesa 15 April 2008
/// Add recovery of xs from pyxsec_hists.root file 15/jan/2015
//______________________________________________________________________________
Bool_t GetAverageXsection(TTree * tree, Double_t & xs, Float_t & ntr, Int_t & n)
{
  Double_t xsection = 0 ;
  UInt_t    ntrials = 0 ;
  Int_t      nfiles = 0 ;
  
  xs  = 0;
  ntr = 0;
  n   = 0;
  if( kInputData != "AOD" &&  tree))
  {
    nfiles =  tree->GetEntries()  ;
    
    tree->SetBranchAddress("xsection",&xsection);
    tree->SetBranchAddress("ntrials" ,&ntrials );
    for(Int_t i = 0; i < nfiles; i++)
    {
      tree->GetEntry(i);
      if(xsection > 0)
      {
        xs  += xsection ;
        ntr += ntrials ;
        n++;
      }
      cout << "xsection " <<xsection<<" ntrials "<<ntrials<<endl;
    } // loop
  }
  else if( kInputData == "AOD" && xsArr))
  {
    nfiles = xsArr->GetSize();
    
    for(Int_t i = 0; i < nfiles; i++)
    {
      if(xsArr->GetAt(i) > 0)
      {
        xs  += xsArr->GetAt(i) ;
        ntr += trArr->GetAt(i) ;
        n++;
      }
      cout << "xsection " <<xsArr->GetAt(i)<<" ntrials "<<trArr->GetAt(i)<<endl;
    } // loop
  }
  else return kFALSE;
  
  xs =   xs /  n;
  ntr =  ntr / n;
  cout << "-----------------------------------------------------------------"<<endl;
  cout << "Average of "<< n <<" files: xsection " <<xs<<" ntrials "<<ntr<<endl;
  cout << "-----------------------------------------------------------------"<<endl;
  
  return kTRUE;
}


