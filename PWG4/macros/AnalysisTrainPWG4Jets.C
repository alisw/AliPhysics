//===================== ANALYSIS TRAIN =========================================
// To use: copy this macro to your work directory, modify the global part to match
// your needs, then run root.
//    root[0] .L AnalysisTrain.C
// Grid full mode as below (other modes: test, offline, submit, terminate)
//    root[1] AnalysisTrainPWG-244Jets("grid", "full")
// CAF mode (requires root v5-23-02 + aliroot v4-16-Rev08)
//    root[2] AnalysisTrainPWG4Jets("proof")
// Local mode requires AliESds.root or AliAOD.root in ./data directory
//    root[3] AnalysisTrainPWG4Jets("local")
// In proof and grid modes, a token is needed and sourcing the produced environment file.
//
// If 'saveTrain' flag is set, the train will generate a directory name and run
// in this directory. A configuration file 'ConfigTrain.C' will be generated. 
// One can replay at any time the train via:
//    root[1] AnalysisTrainPWG4Jets(ana_mode, plugin_mode, "train_default_<date>/ConfigTrain.C")

// For Usage with fastjet run in "offline" first and then "submit"
// jdl and analysismacro are automatically patched after "offline" mode

// =============================================================================
// ### General Steering variables
// =============================================================================
//== general setup variables
TString     kTrainName         = "testAnalysis"; // *CHANGE ME* (no blancs or special characters)
TString     kJobTag            = "PWG4 Jet Tasks analysis train configured"; // *CHANGE ME*

// Usage of par files ONLY in grid mode and ONLY if the code is not available
// in the deployed AliRoot versions. Par file search path: local dir, if not there $ALICE_ROOT.
// To refresh par files, remove the ones in the workdir, then do "make <target.par>" in 
// AliRoot.
Bool_t      kUsePAR             = kFALSE;  // use par files for extra libs
Bool_t      kUseCPAR            = kFALSE;  // use par files for common libs
Bool_t      kFillAOD = kTRUE;  // switch of AOD filling for on the fly analysis
Bool_t      kFilterAOD = kTRUE;
Float_t     kJetTriggerPtCut = 20; // pT for jet trigger in case of iFilter==2
Int_t       kSaveAOD = 8;        // Bit switch 1 = Full AOD 2 = Jet AOD , 4 = PartCorr, 8 = JCORRAN 
//== general input and output variables

Int_t       iAODanalysis       = 1;      // Analysis on input AOD's
Int_t       iFilterAnalysis       = 0;      // Analysis on input AOD's
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iCentralitySelection  = 0;      // Use the centrality
Int_t       iESDfilter         = 0;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iPhysicsSelection  = 1;      // ESD to AOD filter (barrel + muon tracks)
UInt_t      iPhysicsSelectionFlag = 1; // set by pyshics selection and passed to the task, kMB, kUserDefined etc
Bool_t      useTender           = kFALSE; // use tender wagon 
Bool_t      kUseKinefilter     = kFALSE; // use Kinematics filter
Bool_t      kUseMuonfilter     = kFALSE; // use Kinematics filter
TString     kCommonOutputFileName = "PWG4_JetTasksOutput.root";
TString     kCaloQAOutputFileName = "PWG4_CaloQAOutput.root";


//== general process variables

// ### Other flags to steer the analysis
//==============================================================================
Bool_t      kSkipTerminate      = kTRUE; // Do not call Teminate
Bool_t      kUseDate            = kFALSE; // use date in train name
Bool_t      kUseDebug           = kTRUE; // activate debugging
Int_t         kErrorIgnoreLevel = -1; // takes the errror print level from .rootrc 
// From TError.h:
// const Int_t kUnset    =  -1;
// const Int_t kPrint    =   0;
// const Int_t kInfo     =   1000;
// const Int_t kWarning  =   2000;
// const Int_t kError    =   3000;
// const Int_t kBreak    =   4000;
// const Int_t kSysError =   5000;
// const Int_t kFatal    =   6000; 
Int_t         kUseSysInfo         = 0; // activate debugging
Long64_t    kNumberOfEvents     = 1234567890; // Number of events to process from the chain
Bool_t      kUseMC              = kTRUE;  // use MC info
Bool_t      kIsMC               = kTRUE;  // is MC info, if false it overwrites Use(AOD)MC
Bool_t      kUseAODMC           = kTRUE;  // use MC infA
Bool_t      kUseESDTags         = kFALSE; // use ESD tags for selection
Bool_t      kUseTR              = kFALSE;  // use track references
Bool_t      kUseAODTags         = kFALSE;  // use AOD tags
Bool_t      kSaveTrain          = kFALSE;  // save train configuration as: 
Bool_t      kIsPbPb             = kFALSE;  // Pb+Pb


// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iJETAN             = 1;      // Jet analysis (PWG4) // 1 write standard 2 write non-standard jets, 3 wrtie both
Int_t       iJETSUBTRACT        = 1;      // Jet background subtration
TList       kJetListSpectrum;             // list of jets contains TObjString of possible jet finder names
TExMap      kJetMapSpectrum;             // Maps the jet finder pairs to be used in the spectrum task second number negative no pair other wise (j1+1) + (1000 * (j2+1)) +10000 * (j3+1)
TExMap      kJetBackMapSpectrum;             // Maps the jet finder pairs with the background branch used, just for countint of trackrefs
Int_t       kJetMapOffset[3] = {10000,100,1};
TString     kDefaultJetBranch     = "";      // is currently set when filled (iJETAN or clusters) or from config macro 
TString     kDefaultJetBackgroundBranch            = "";      // is currently set when filled (jet clsuters  
TString     kDefaultJetBackgroundBranchCut1        = "";      // is currently set when filled (jet clsuters  
TString     kDefaultJetBackgroundBranchCut2        = "";      // is currently set when filled (jet clsuters  
TString     kDefaultJetBackgroundBranch_extra     = "";      // is currently set when filled (jet clsuters) 
TString     kJetSubtractBranches     = "";      // is currently set when filled (jet clsuters  
TString     kJetSubtractBranchesCut1     = "";      // is currently set when filled (jet clsuters  
TString     kJetSubtractBranchesCut2     = "";      // is currently set when filled (jet clsuters  
TString     kJetSubtractBranches_extra     = "";      // is currently set when filled (jet clsuters  

TString     kDefaultJetBranchMC     = "";      // is currently set when filled (iJETAN or clusters) or from config macro 
TString     kDefaultJetBackgroundBranchMC     = "";      // is currently set when filled (jet clsuters  
TString     kDefaultJetBranchMC2     = "";      // is currently set when filled (iJETAN or clusters) or from config macro 
TString     kDefaultJetBackgroundBranchMC2     = "";      // is currently set when filled (jet clsuters  
TString     kJetSubtractMask1 = "B0";
TString     kJetSubtractMask2 = "B%d";
Int_t       iDIJETAN           = 1;
Int_t       iJETANLib          = 1;
Int_t       iPWG1QASym         = 0;      // Eva's QA task compiled on the fly...
Int_t       iPWG4FastEmbedding = 0;      // Generate non-standard AOD for embedding
Int_t       iPWG4JetTasks      = 0;      // all jet tasks flag for lib laoding
Int_t       iPWG4JetServices   = 0;      // jet spectrum analysis
Int_t       iPWG4JetSpectrum   = 0;      // jet spectrum analysis
Int_t       iPWG4JetResponse   = 0;      // jet response matrix
Int_t       iPWG4JCORRAN       = 0;      // JCORRAN module
Int_t       iPWG4UE            = 0;      // Underlying Event analysis
Int_t       iPWG4LeadingUE     = 0;      // Underlying Event analysis
Int_t       iPWG4CorrectionsUE = 0;      // Underlying Event analysis
Int_t       iPWG4TmpSourceSara = 0;      // Underlying Event analysis not in svn
Int_t       iPWG4Fragmentation = 1;      // Official Fragmentation
Int_t       iPWG4JetChem       = 0;      // Jet chemistry 
Int_t       iPWG4PtQAMC        = 0;      // Marta's QA tasks 
Int_t       iPWG4PtTrackQA     = 0;      // Marta's QA tasks  
Int_t       iPWG4PtSpectra     = 0;      // Marta's QA tasks 
Int_t       iPWG4PtQATPC       = 0;      // Marta's QA tasks 
Int_t       iPWG4Cosmics     = 0;      // Marta's Cosmics Taks 
Int_t       iPWG4ThreeJets     = 0;      // Sona's thrust task
Int_t       iPWG4QGSep     = 0;          // Sona's QG Separation task
Int_t       iPWG4Minijet       = 0;      // Eva's Mini Jet Task cluster task 
Int_t       iPWG4KMeans        = 0;      // Andreas' KMeans task 
Int_t       iPWG4Cluster       = 0;      // CKB cluster task 
Int_t       iEMCUtilLibs       = 0;      // Flag to satisfy dependence on EMC utils
Int_t       iPWG4PartCorrLibs  = 0;      // Gustavo's part corr analysis
Int_t       iPWG4PartCorr      = 0;      // Gustavo's part corr analysis
Int_t       iPWG4CaloQA        = 0;      // Gustavo's part corr analysis
Int_t       iPWG4JetCorr       = 0;     // Paul's jet corr analysis
Int_t       iPWG4Tagged        = 0;      // Gustavo's part corr analysis
Int_t       iPWG4omega3pi      = 0;      // Omega to 3 pi analysis (PWG4) 
Int_t       iPWG4GammaConvLib     = 0;      // Gamma Conversio
Int_t       iPWG4GammaConv     = 0;      // Gamma Conversio
Int_t       iPWG4CaloConv     = 0;      // Gamma Conversio
Int_t       kHighPtFilterMask  = 32;     // change depending on the used AOD Filter
TString     kDeltaAODJetName   = "AliAOD.Jets.root";     
TString     kDeltaAODJetNameInput   = "";     
TString     kDeltaAODJCORRANName   = "AliAOD.JCORRAN.root";     
TString     kDeltaAODPartCorrName   = "AliAOD.PartCorr.root";     
TString     kFastEmbeddingAOD  = "emb/AliAOD.root";
TString     kFastEmbAODList    = "";

//==============================================================================
// ### PROOF Steering varibales
//==============================================================================
//== proof setup variables
TString     kProofCluster      = "alicecaf.cern.ch";
Bool_t      kProofUseAFPAR     = kFALSE;  // use AF special par file
TString     kProofAFversion          = "AF-v4-17";
//== proof input and output variables
TString     kProofDataSet      = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
Bool_t      kProofSaveToAlien   = kFALSE; // save proof outputs in AliEn train_[trainName]_ddMonthyyyy_time.C
TString     kProofOutdir       = "";
Bool_t      kProofClearPackages = kFALSE;
Int_t       kProofEvents = 10000;
Int_t       kProofOffset = 0;
//== proof process variables


//==============================================================================
// ### Grid plugin Steering varibiables
//==============================================================================
//== grid plugin setup variables
Bool_t      kPluginUse         = kTRUE;   // do not change
Bool_t      kPluginUseProductionMode  = kFALSE;   // use the plugin in production mode
TString     kPluginRootVersion       = "v5-30-02";  // *CHANGE ME IF MORE RECENT IN GRID*
TString     kPluginAliRootVersion    = "v4-21-01a-AN";  // *CHANGE ME IF MORE RECENT IN GRID*                                          
Bool_t      kPluginMergeViaJDL       = kTRUE;  // merge via JDL
Bool_t      kPluginFastReadOption   = kFALSE;  // use xrootd tweaks
Bool_t      kPluginOverwriteMode    = kTRUE;  // overwrite existing collections
Int_t       kPluginOutputToRunNumber = 1;     // write the output to subdirs named after run number
// TString kPluginExecutableCommand = "root -b -q";
TString     kPluginExecutableCommand = "source /Users/kleinb/setup_32bit_aliroot_trunk_clean_root_trunk.sh; alienroot -b -q ";

// == grid plugin input and output variables
TString     kGridDatadir      = "/alice/sim/PDC_08b/LHC09a1/AOD/";
Int_t         kGridMaxRunsFromList = 999999999;
Int_t         kGridOffsetRunFromList = 0; // skip the first n runs from the list
TString     kGridLocalRunList = "";
TString     kGridOutdir       = ""; // AliEn output directory. If blank will become output_<kTrainName>
TString     kGridDataSet      = ""; // sub working directory not to confuse different run xmls 
TString     kGridExtraAliendirLevel = ""; // sub working directory not to confuse different run xmls 
Int_t         kGridRunRange[2]       =  {0, -1}; // Set the run range
TString     kGridRunPattern        = "%03d"; // important for leading zeroes!!
TString     kGridPassPattern       = "";
TString     kGridExtraFiles        = ""; // files that will be added to the input list in the JDL...
Int_t         kGridMaxMergeFiles      = 25; // Number of files merged in a chunk grid run range
TString     kGridMergeExclude       = "AliAOD.root"; // Files that should not be merged
TString     kGridOutputStorages      = "disk=2"; // Make replicas on the storages
// == grid process variables
Int_t       kGridRunsPerMaster     = 100; // Number of runs per master job
Int_t       kGridFilesPerJob       = 100; // Maximum number of files per job (gives size of AOD)

//==============================================================================
// ### Local Steering variables
//==============================================================================
//== local setup variables
//== local input and output variables
TString     kLocalXMLDataset   = ""; // Change local xml dataset for local interactive analysis
TString     kLocalDataList   = "local_deltaaod.txt"; // Change local xml dataset for local interactive analysis
// == local process variables

TString kPluginMode = "";
TString kAnalysisMode = "";


// Temporaries.
TString anaPars = "";
TString anaLibs = "";
TString anaLibsExtra = "";
TString anaSources = "";
// Function signatures
class AliAnalysisAlien;
class AliAnalysisManager;

//______________________________________________________________________________
void AnalysisTrainPWG4Jets(const char *analysis_mode="local", 
			   const char *plugin_mode="",
			   const char *config_file="",Int_t iOffset = 0,Int_t iTotal = 0)
{
// Main analysis train macro. If a configuration file is provided, all parameters
// are taken from there but may be altered by CheckModuleFlags.

  // these flag may be needed by the config file
  kPluginMode   = plugin_mode;
  kAnalysisMode = analysis_mode;
  
  if (strlen(config_file) && !LoadConfig(config_file)) return;
  if(iTotal>0)kGridMaxRunsFromList = iTotal; // overwrites the settings from config file
  if(iOffset)kGridOffsetRunFromList = iOffset;
  if(iOffset)kProofOffset = iOffset;


   TString smode(analysis_mode);
   smode.ToUpper();
   if (kSaveTrain)WriteConfig();
   // Check compatibility of selected modules 
   CheckModuleFlags(smode);
   //     gROOT->ProcessLine(".trace");

   printf("==================================================================\n");
   printf("===========    RUNNING ANALYSIS TRAIN %s IN %s MODE   ==========\n", kTrainName.Data(),smode.Data());
   printf("==================================================================\n");
                              printf("=  Configuring analysis train for:                               =\n");
   if (iAODanalysis) printf("=  AOD analysis                                                  =\n");
   else                     printf("=  ESD analysis                                                  =\n");
   if (iPhysicsSelection)   printf("=  Physics selection                                                =\n");
   if(iCentralitySelection)printf("=  Centrality selection                                                =\n");
   if (useTender)   printf("=  Using tender                                                =\n");
   if (iESDfilter)   printf("=  ESD filter                                                    =\n");
   if (iJETAN)       printf("=  Jet analysis                                                  =\n");
   printf("==================================================================\n");

   char *printMask = ":: %20s  %10d\n";
   printf(printMask,"Fill AOD",(UInt_t)kFillAOD);
   printf(printMask,"Save AOD", (UInt_t)kSaveAOD);
   printf(printMask,"MC truth", (UInt_t)kUseMC);
   printf(printMask,"KINE filter", (UInt_t)kUseKinefilter);
   printf(printMask,"track refs", (UInt_t)kUseTR);
   printf(printMask,"tags", (UInt_t)kUseESDTags);
   printf(printMask,"AOD tags", (UInt_t)kUseAODTags);
   printf(printMask,"debugging", (UInt_t)kUseDebug);
   printf(printMask,"PAR files", (UInt_t)kUsePAR);
   printf(printMask,"AliEn plugin", (UInt_t)kPluginUse);
   printf(printMask,"JETAN subtract", (UInt_t)iJETSUBTRACT);
   printf(printMask,"PWG1 QA sym", iPWG1QASym);
   printf(printMask,"PWG4 Source Sara",iPWG4TmpSourceSara);
   printf(printMask,"PWG4 Fragmentation",iPWG4Fragmentation);
   printf(printMask,"PWG4 Jet Chem",iPWG4JetChem);
   printf(printMask,"PWG4 Jet tasks",iPWG4JetTasks);
   printf(printMask,"PWG4 Jet Services",iPWG4JetServices);     
   printf(printMask,"PWG4 Jet Spectrum",iPWG4JetSpectrum);
   printf(printMask,"PWG4 JCORRAN",iPWG4JCORRAN);
   printf(printMask,"PWG4 UE",iPWG4UE); 
   printf(printMask,"PWG4 Leading UE",iPWG4LeadingUE); 
   printf(printMask,"PWG4 Corrections UE",iPWG4CorrectionsUE); 
   printf(printMask,"PWG4 Pt QA MC",iPWG4PtQAMC);
   printf(printMask,"PWG4 Pt QA track",iPWG4PtTrackQA);
   printf(printMask,"PWG4 Pt Spectra",iPWG4PtSpectra);
   printf(printMask,"PWG4 Pt QA TPC",iPWG4PtQATPC);     
   printf(printMask,"PWG4 Cosmics",iPWG4Cosmics);     
   printf(printMask,"PWG4 Three Jets",iPWG4ThreeJets);
   printf(printMask,"PWG4 QGSep",iPWG4QGSep);
   printf(printMask,"PWG4 Minijet",iPWG4Minijet);
   printf(printMask,"PWG4 KMeans",iPWG4KMeans);
   printf(printMask,"PWG4 Cluster",iPWG4Cluster);
   printf(printMask,"PWG4 Part Corr",iPWG4PartCorr);
   printf(printMask,"PWG4 Calo QA",iPWG4CaloQA);
   printf(printMask,"PWG4 Jet Corr",iPWG4JetCorr);
   printf(printMask,"PWG4 JCORRAN",iPWG4JCORRAN);
   printf(printMask,"PWG4 Tagged",iPWG4Tagged);
   printf(printMask,"PWG4 omega to 3 pi ",iPWG4omega3pi);
   printf(printMask,"PWG4 Gamma Conv",iPWG4GammaConv);
   printf(printMask,"PWG4 Calo Conv",iPWG4CaloConv);
   printf(printMask,"HighPt FilterMask",kHighPtFilterMask);    
   
   //==========================================================================
   // Connect to back-end system
   if (!Connect(smode)) {
      ::Error("AnalysisTrain", "Could not connect to %s back-end", analysis_mode);
      return;
   }   

   // Load common libraries and set include path
   if (!LoadCommonLibraries(smode)) {
      ::Error("AnalysisTrain", "Could not load common libraries");
      return;
   }


    
   // Make the analysis manager and connect event handlers
   AliAnalysisManager *mgr  = new AliAnalysisManager("PWG4Train", "pwg4 mini train");
   if (kCommonOutputFileName.Length()>0)mgr->SetCommonFileName(kCommonOutputFileName.Data());
   if (kProofSaveToAlien) mgr->SetSpecialOutputLocation(kProofOutdir);
   mgr->SetNSysInfo(0);
   if (!strcmp(plugin_mode, "test")) mgr->SetNSysInfo(1);
   if (kUseSysInfo)mgr->SetNSysInfo(kUseSysInfo);
   mgr->SetSkipTerminate(kSkipTerminate);

   // Load analysis specific libraries
   if (!LoadAnalysisLibraries(smode)) {
      ::Error("AnalysisTrain", "Could not load analysis libraries");
     return;
   }   

   // Create input handler (input container created automatically)
   if (iAODanalysis) {
   // AOD input handler
      AliAODInputHandler *aodH = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodH);
      if (kDeltaAODJetNameInput.Length()){
	Printf("Adding Friend %s",kDeltaAODJetNameInput.Data());
	aodH->AddFriend(kDeltaAODJetNameInput.Data());
      }
      //      if (iPWG4PartCorr) aodH->AddFriend(Form("deltas/%s"kDeltaAODJetName.Data()));
   } else {   
   // ESD input handler
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      if (kUseESDTags) esdHandler->SetReadTags();
      esdHandler->SetReadFriends(kFALSE);
      mgr->SetInputEventHandler(esdHandler);       
      //      if(iPWG4PtQATPC&& !kTrainName.Contains("pass5"))esdHandler->SetActiveBranches("ESDfriend");

   }

   // Monte Carlo handler
   if (kUseMC && !iAODanalysis) {
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
      mcHandler->SetReadTR(kUseTR); 
   }   
   // AOD output container, created automatically when setting an AOD handler
   if (iAODhandler) {
      // AOD output handler
      AliAODHandler* aodHandler   = new AliAODHandler();
      aodHandler->SetOutputFileName("AliAOD.root");
      aodHandler->SetFillAODforRun(kFillAOD);
      //
      mgr->SetOutputEventHandler(aodHandler);
      AliAnalysisDataContainer * cout_aod = mgr->GetCommonOutputContainer();
      cout_aod->SetSpecialOutput();
   }
   // Debugging if needed

   if (kUseDebug){
     mgr->SetDebugLevel(3);
   }
   if(kUseSysInfo>0){
     mgr->RegisterExtraFile("syswatch.root");
     if(kGridMergeExclude.Length())kGridMergeExclude += " ";
     kGridMergeExclude += "syswatch.root";
   }
   else{
     AliLog::SetGlobalLogLevel(AliLog::kError);
   }


   //==========================================================================
   // Create the chain. In this example it is created only from ALIEN files but
   // can be done to work in batch or grid mode as well.
   TChain *chain = CreateChain(smode, plugin_mode);
        
   //==========================================================================
   // Load the tasks configuration macros for all wagons. These files are supposed now to be
   // in the current workdir, but in AliEn they will be in the file catalog, 
   // mapped from AliRoot and pecified in the jdl input list.
    
   // For now connection to top input container and common AOD output container
   // is done in this macro, but in future these containers will be connected
   // from each task configuration macro.
   
   //                                                                                                                              
   // Tender and supplies. Needs to be called for every event.                                                                     
   //                                                                                                                              

   if (useTender) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C");
      AliAnalysisTaskSE *tender = AddTaskTender(kTRUE);
      //      tender->SelectCollisionCandidates();                                                                                      
      tender->SetDebugLevel(2);
   }

   
   Float_t fTrackEtaWindow = 0.9;
   Float_t fJetEtaWindow   = 0.5;
   Float_t fVertexWindow =  10;
   /*
   if(kIsPbPb){// for pass1
     Float_t fTrackEtaWindow = 0.8;
     Float_t fJetEtaWindow   = 0.4;
   }
   */


   if(iPhysicsSelection && !iAODanalysis){
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
     Int_t iTriggerHIC = 0;
     Bool_t rejectBkg = true;
     if(kIsPbPb){
       iTriggerHIC = 2;
       rejectBkg = false; // for the moment...
     }
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kIsMC,true,rejectBkg);  

     //     mgr->RegisterExtraFile("event_stat.root");
     mgr->AddStatisticsTask();
   }
   else{
     iPhysicsSelectionFlag = AliVEvent::kMB;
   }
   
   if(kIsPbPb&&!iAODanalysis){

     // has to run before AOD filter
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
     AliCentralitySelectionTask *taskC = AddTaskCentrality();
      if (!taskC) ::Warning("AnalysisTrainPWG4Jets", "AliCentralitySelectionTask cannot run for this train conditions - EXCLUDED");


     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
     AliEPSelectionTask *taskEP = AddTaskEventplane();
     if (!taskEP) ::Warning("AnalysisTrainPWG4Jets", "AliEventplan cannot run for this train conditions - EXCLUDED");
   }
   
   if (iESDfilter && !iAODanalysis) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskESDFilterPWG4Train.C");
      // switch on centrality make for PbPb
      AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilterPWG4Train(kUseKinefilter); // carefull, if physics selection is enabled you may get not primary vertex pointer later on...
      taskesdfilter->SetEnableFillAOD(!kFilterAOD);
      taskesdfilter->DisableV0s();
      taskesdfilter->DisableCascades();
      taskesdfilter->DisableKinks();
      taskesdfilter->DisablePmdClusters();
      taskesdfilter->DisableCaloClusters();
      taskesdfilter->DisableCells();

      if(kIsMC){
	mgr->RegisterExtraFile("pyxsec_hists.root");
	if(kGridMergeExclude.Length())kGridMergeExclude += " ";
	kGridMergeExclude += "pyxsec_hists.root";

      }
   }   

   // AOD tags
   if (kUseAODTags) {
      AliAnalysisTaskTagCreator* tagTask = new AliAnalysisTaskTagCreator("AOD Tag Creator");
      mgr->AddTask(tagTask);
      AliAnalysisDataContainer *coutTags = mgr->CreateContainer("cTag",  TTree::Class(), 
                                           AliAnalysisManager::kOutputContainer, "AOD.tag.root");
      mgr->ConnectInput (tagTask, 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(tagTask, 1, coutTags);
   }   

   if (iPWG4FastEmbedding) {
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskFastEmbedding.C");
     AliAnalysisTaskFastEmbedding *taskEmbedding = 0;
     if(kFastEmbAODList.Length()) taskEmbedding = AddTaskFastEmbedding(kFastEmbAODList, 1);
     else                         taskEmbedding = AddTaskFastEmbedding(kFastEmbeddingAOD, 0);
     taskEmbedding->SetJetBranch("jetsAOD_UA104_B0_Filter00256_Cut01000");
     // taskEmbedding->SetEvtSelecMode(AliAnalysisTaskFastEmbedding::kEventsAll);
     // taskEmbedding->SetDebugLevel(10);
   }

    // Jet analysis
   if (iJETAN) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJets.C");
      AliAnalysisTaskJets *taskjets = 0;
      if(iJETAN&1){
	/*
	taskjets = AddTaskJets(kHighPtFilterMask); 
	taskjets->SetName("StandardJets");
	taskjets->SetNonStdBranch("");
	*/
      }
      if(iJETAN&2){
	// Set only few jet finders  backgroudn subtraction w an wo 
	
	taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,1.,0); // no background subtraction
	if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
	TString cTmp("");	
	cTmp = taskjets->GetNonStdBranch();
	if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));

	taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,2.,0); // no background subtraction
	if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
	TString cTmp("");	
	cTmp = taskjets->GetNonStdBranch();
	if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));


	// STANDARD UA jet finders pT cut 1 GeV background mode 2 R = 0.4
	if(kIsPbPb){
	  taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,1.,2); // background subtraction
	  if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
	  cTmp = taskjets->GetNonStdBranch();
	  if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	  taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,2.,2); // background subtraction
	  if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
	  cTmp = taskjets->GetNonStdBranch();
	  if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	}



	// SICONE 
	/*
	taskjets = AddTaskJets("AOD","SISCONE",0.4,kHighPtFilterMask,0.15,0); //no background subtraction to be done later....
	if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
	cTmp = taskjets->GetNonStdBranch();
	if(cTmp.Length()>0)kJetSubtractBranches += Form("%s ",cTmp.Data());
	if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	*/
	if(kUseAODMC){
	  // STANDARD UA jet finders pT cut 1 GeV background mode 2 R = 0.4
	  if(kIsPbPb){
	    taskjets = AddTaskJets("AODMC","UA1",0.4,kHighPtFilterMask,1.,2); // background subtraction
	    cTmp = taskjets->GetNonStdBranch();
	    if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	    taskjets = AddTaskJets("AODMC2","UA1",0.4,kHighPtFilterMask,1.,2); // background subtraction
	    cTmp = taskjets->GetNonStdBranch();
	    if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	  }
	  else{
	    taskjets = AddTaskJets("AODMC","UA1",0.4,kHighPtFilterMask,1.,0); // no background subtraction
	    cTmp = taskjets->GetNonStdBranch();
	    if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	    taskjets = AddTaskJets("AODMC2","UA1",0.4,kHighPtFilterMask,1.,0); // no background subtraction
	    cTmp = taskjets->GetNonStdBranch();
	    if(cTmp.Length())kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	  }
	}
	if(kDeltaAODJetName.Length()>0)mgr->RegisterExtraFile(kDeltaAODJetName.Data()); 
      }
      if (!taskjets) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
   }


   if (iPWG4FastEmbedding && iJETAN) {
     AliAnalysisTaskJets *taskEmbJets = AddTaskJets("AODextra", "FASTJET", 0.4, kHighPtFilterMask);
     taskEmbJets->ReadAODFromOutput();
     kJetSubtractBranches_extra += Form("%s ", taskEmbJets->GetNonStdBranch());

     taskEmbJets = AddTaskJets("AODextraonly", "FASTJET", 0.4, kHighPtFilterMask);
     taskEmbJets->ReadAODFromOutput();
    
     taskEmbJets = AddTaskJets("AODextra", "UA1", 0.4, kHighPtFilterMask,1.,0);
     taskEmbJets->ReadAODFromOutput();
     taskEmbJets = AddTaskJets("AODextraonly", "UA1", 0.4, kHighPtFilterMask,1.,0);
     taskEmbJets->ReadAODFromOutput();
     taskEmbJets = AddTaskJets("AODextra", "UA1", 0.4, kHighPtFilterMask,1.,2);
     taskEmbJets->ReadAODFromOutput();
   }

   if(iPWG4Cluster){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetCluster.C");
     AliAnalysisTaskJetCluster *taskCl = 0;
     Float_t fCenUp = 0;
     Float_t fCenLo = 0;
     if(kIsPbPb&&!kIsMC){
       fCenLo = 0;
       fCenUp = 80;
     }
     if(iPWG4Cluster&1){

       if(kIsPbPb){
	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow,0); // this one is for the background and random jets, random cones with no skip
	 taskCl->SetBackgroundCalc(kTRUE);
	 taskCl->SetNRandomCones(1);
	 //	 taskCl->SetDebugLevel(11);
	 taskCl->SetCentralityCut(fCenLo,fCenUp);

	 taskCl->SetJetTypes(1<<0|1<<2|1<<3);
	 kDefaultJetBackgroundBranch = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),0)));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomCone_random",taskCl->GetJetOutputBranch())));

	 kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomConeSkip00");
	 kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomCone_random");

	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),1.0,fTrackEtaWindow,fVertexWindow,0); // this one is for the background and random jets, random cones with no skip
	 taskCl->SetNRandomCones(1);
	 taskCl->SetBackgroundCalc(kTRUE);
	 taskCl->SetCentralityCut(fCenLo,fCenUp);
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);

	 taskCl->SetJetTypes(1<<0|1<<2|1<<3);
	 kDefaultJetBackgroundBranchCut1 = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	 kJetSubtractBranchesCut1 += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomConeSkip00");
	 kJetSubtractBranchesCut1 += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomCone_random");
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),0)));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomCone_random",taskCl->GetJetOutputBranch())));

	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),2.0,fTrackEtaWindow,fVertexWindow,0); // this one is for the background and random jets, random cones with no skip
	 taskCl->SetNRandomCones(1);
	 taskCl->SetBackgroundCalc(kTRUE);
	 taskCl->SetCentralityCut(fCenLo,fCenUp);
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);

	 taskCl->SetJetTypes(1<<0|1<<2|1<<3);
	 kDefaultJetBackgroundBranchCut2 = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	 kJetSubtractBranchesCut2 += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomConeSkip00");
	 kJetSubtractBranchesCut2 += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomCone_random");
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),0)));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomCone_random",taskCl->GetJetOutputBranch())));


         if (iPWG4FastEmbedding) {
           AliAnalysisTaskJetCluster *taskClEmb = 0;
           taskClEmb = AddTaskJetCluster("AODextra","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background and random jets
           taskClEmb->SetBackgroundCalc(kTRUE);
           taskClEmb->SetCentralityCut(fCenLo,fCenUp);
	   if(iAODanalysis==2)taskClEmb->SetAODTrackInput(kTRUE);
	   kDefaultJetBackgroundBranch_extra = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskClEmb->GetJetOutputBranch());

           taskClEmb = AddTaskJetCluster("AODextraonly","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background and random jets
           taskClEmb->SetBackgroundCalc(kFALSE);
           taskClEmb->SetCentralityCut(fCenLo,fCenUp);
	   if(iAODanalysis==2)taskClEmb->SetAODTrackInput(kTRUE);
	   
           taskClEmb = AddTaskJetCluster("AODextra","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow);
           taskClEmb->SetCentralityCut(fCenLo,fCenUp);
           taskClEmb->SetBackgroundBranch(kDefaultJetBackgroundBranch_extra.Data());
           kJetSubtractBranches_extra += Form("%s ",taskClEmb->GetJetOutputBranch());
	   if(iAODanalysis==2)taskClEmb->SetAODTrackInput(kTRUE);

           taskClEmb = AddTaskJetCluster("AODextraonly","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow);
           taskClEmb->SetCentralityCut(fCenLo,fCenUp);
	   if(iAODanalysis==2)taskClEmb->SetAODTrackInput(kTRUE);
         }

	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.2,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow,0); // this one is for the background and random jets
	 taskCl->SetNRandomCones(1);
	 taskCl->SetBackgroundCalc(kTRUE);
	 taskCl->SetCentralityCut(fCenLo,fCenUp);
	 taskCl->SetJetTypes(1<<0|1<<2|1<<3);
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),0)));
	 kJetListSpectrum.Add(new TObjString(Form("%sRandomCone_random",taskCl->GetJetOutputBranch())));



	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
       }
       else{
	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.6,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background jets
	 taskCl->SetBackgroundCalc(kTRUE);
	 kDefaultJetBackgroundBranch = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));

	 taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); 
	 taskCl->SetBackgroundCalc(kTRUE);
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       } 

       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // thid id sldo got the random cones but avoiding the two leading jets
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());
       taskCl->SetNRandomCones(1);
       if(iFilterAnalysis==2){
	 taskCl->SetJetTriggerPtCut(kJetTriggerPtCut);
       }
       if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);

       taskCl->SetJetTypes(1<<0|1<<2); // only store the RC on full event with 2 removed
       kDefaultJetBranch = taskCl->GetJetOutputBranch();
       kJetSubtractBranches += Form("%s ",taskCl->GetJetOutputBranch());
       //       kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomConeSkip02");
       //       kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomCone_random");
       kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),2)));
       //       kJetListSpectrum.Add(new TObjString(Form("%sRandomCone_random",taskCl->GetJetOutputBranch())));

       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1,kDeltaAODJetName.Data(),1.0,fTrackEtaWindow,fVertexWindow);
       taskCl->SetNRandomCones(1);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranchCut1.Data());

       taskCl->SetJetTypes(1<<0|1<<2); // only store the RC on full event with 2 removed
       kJetSubtractBranchesCut1 += Form("%s ",taskCl->GetJetOutputBranch());
       kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),2)));


       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1,kDeltaAODJetName.Data(),2.0,fTrackEtaWindow,fVertexWindow);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       taskCl->SetNRandomCones(1);
       if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranchCut2.Data());
       taskCl->SetJetTypes(1<<0|1<<2); // only store the RC on full event with 2 removed
       kJetSubtractBranchesCut2 += Form("%s ",taskCl->GetJetOutputBranch());
       kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),2)));

       // tmp track qa...
       /*
       taskCl = AddTaskJetCluster("AOD","",1<<8,iPhysicsSelectionFlag,"ANTIKT",0.4,2,1,kDeltaAODJetName.Data(),2.0);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       taskCl->SetFilterMask(1<<4|1<<8,1);
       */
       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.2,0,1,kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow,2);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       taskCl->SetNRandomCones(1);
       if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());

       taskCl->SetJetTypes(1<<0|1<<2); 
       kJetSubtractBranches += Form("%s ",taskCl->GetJetOutputBranch());
       //       kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomConeSkip00");
       //       kJetSubtractBranches += Form("%s%s ",taskCl->GetJetOutputBranch(),"RandomCone_random");
       kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       kJetListSpectrum.Add(new TObjString(Form("%sRandomConeSkip%02d",taskCl->GetJetOutputBranch(),2)));


       if(kUseAODMC){
	 if(kIsPbPb){
	   taskCl = AddTaskJetCluster("AODMC","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background and random jets
	   taskCl->SetBackgroundCalc(kTRUE);
	   kDefaultJetBackgroundBranchMC = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	   if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);	   
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));

	 taskCl = AddTaskJetCluster("AODMC2","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fVertexWindow,fVertexWindow); // this one is for the background and random jets
	   taskCl->SetBackgroundCalc(kTRUE);
	   kDefaultJetBackgroundBranchMC2 = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch()); 
	   if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	   kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	 }
	 else{
	   taskCl = AddTaskJetCluster("AODMC","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.6,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background and random jets
	   taskCl->SetBackgroundCalc(kTRUE);
	   kDefaultJetBackgroundBranchMC = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));	   

	 taskCl = AddTaskJetCluster("AODMC2","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.6,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); // this one is for the background and random jets
	   taskCl->SetBackgroundCalc(kTRUE);
	   kDefaultJetBackgroundBranchMC2 = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch()); 
	   if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	   kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
	   // pp background calcs...
	 }
	 
	 taskCl = AddTaskJetCluster("AODMC","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow); 
	 if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranchMC.Data());	 
	 kDefaultJetBranchMC = taskCl->GetJetOutputBranch();
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));

	 taskCl = AddTaskJetCluster("AODMC2","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,fVertexWindow);
	 if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranchMC2.Data());	 
	 kDefaultJetBranchMC2 = taskCl->GetJetOutputBranch();
	 if(iAODanalysis==2)taskCl->SetAODTrackInput(kTRUE);
	 kJetListSpectrum.Add(new TObjString(taskCl->GetJetOutputBranch()));
       }
       
       
       if (!taskCl) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskCluster cannot run for this train conditions - EXCLUDED");
     }
   }


   if(iJETSUBTRACT&&kJetSubtractBranches.Length()){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetBackgroundSubtract.C");
     AliAnalysisTaskJetBackgroundSubtract *taskSubtract = 0;
     if(kJetSubtractBranches.Length()){




       for(int iB = 1;iB<3;iB++){
	 if(iB>2)continue;
	 taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranches,iB,kJetSubtractMask1.Data(),kJetSubtractMask2.Data());
	 taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());	 	taskSubtract->SelectCollisionCandidates(iPhysicsSelectionFlag);
	 if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());
	 taskSubtract->SetKeepJets(kTRUE);
	 TString cTmp;       
	 TObjArray *objArr = kJetSubtractBranches.Tokenize(" ");
	 for(int iJB = 0;iJB<objArr->GetEntries();iJB++){
	   TObjString *ostr = (TObjString*)objArr->At(iJB);
	   cTmp = ostr->GetString().Data();
	   cTmp.ReplaceAll(kJetSubtractMask1.Data(),Form(kJetSubtractMask2.Data(),iB));
	   kJetListSpectrum.Add(new TObjString(cTmp.Data()));
	 }
	 
	 //	taskSubtract->SetDebugLevel(3);
	 if(iB==2){
	   if(kJetSubtractBranches.Contains(kDefaultJetBranch.Data())&&kIsPbPb){
	     kDefaultJetBranch.ReplaceAll(taskSubtract->GetToReplace(),Form(taskSubtract->GetReplacementMask(),taskSubtract->GetSubtractionMethod()));
	   }
	 }
	 
       }

       //
       // cut1
       Int_t iB = 2;
       taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranchesCut1,iB,kJetSubtractMask1.Data(),kJetSubtractMask2.Data(),"Cut1000");
       taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranchCut1.Data());	 	
       taskSubtract->SelectCollisionCandidates(iPhysicsSelectionFlag);
       if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());

       objArr = kJetSubtractBranchesCut1.Tokenize(" ");
       for(int iJB = 0;iJB<objArr->GetEntries();iJB++){
	 TObjString *ostr = (TObjString*)objArr->At(iJB);
	 cTmp = ostr->GetString().Data();
	 cTmp.ReplaceAll(kJetSubtractMask1.Data(),Form(kJetSubtractMask2.Data(),iB));
	 kJetListSpectrum.Add(new TObjString(cTmp.Data()));
       }

       taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranchesCut2,iB,kJetSubtractMask1.Data(),kJetSubtractMask2.Data(),"Cut2000");
       taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranchCut2.Data());	 	
       taskSubtract->SelectCollisionCandidates(iPhysicsSelectionFlag);
       if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());

       objArr = kJetSubtractBranchesCut2.Tokenize(" ");
       for(int iJB = 0;iJB<objArr->GetEntries();iJB++){
	 TObjString *ostr = (TObjString*)objArr->At(iJB);
	 cTmp = ostr->GetString().Data();
	 cTmp.ReplaceAll(kJetSubtractMask1.Data(),Form(kJetSubtractMask2.Data(),iB));
	 kJetListSpectrum.Add(new TObjString(cTmp.Data()));
       }

     }
     if(kJetSubtractBranches_extra.Length()){
       taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranches_extra,2,kJetSubtractMask1.Data(),kJetSubtractMask2.Data(),"extra");
       taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranch_extra.Data());
       taskSubtract->SelectCollisionCandidates(iPhysicsSelectionFlag);
       //taskSubtract->SetDebugLevel(3);
       if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());
       if(kJetSubtractBranches_extra.Contains(kDefaultJetBranch.Data())){
	 kDefaultJetBranch.ReplaceAll(taskSubtract->GetToReplace(),Form(taskSubtract->GetReplacementMask(),taskSubtract->GetSubtractionMethod()));
       }
     }
     if (!taskSubtract) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetBackgroundSubtrac cannot run for this train conditions - EXCLUDED");     
   }

   if (iDIJETAN) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskDiJets.C");
      AliAnalysisTaskDiJets *taskdijets = 0;
      if(iDIJETAN&1)taskdijets = AddTaskDiJets(); 
      if (!taskdijets) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
      if(iDIJETAN&2){
	taskdijets = AddTaskDiJets("jetsAOD_CDF07"); 
	taskdijets = AddTaskDiJets("jetsAOD_DA07"); 
	taskdijets = AddTaskDiJets("jetsAOD_FASTJET07"); 
	taskdijets = AddTaskDiJets("jetsAOD_FASTKT07"); 
	taskdijets = AddTaskDiJets("jetsAOD_SISCONE07"); 
	taskdijets = AddTaskDiJets("jetsAOD_UA107");
      }
   }

   if(iPWG1QASym){
     gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskQAsym.C");
     AliAnalysisTaskQASym *taskQASym = AddTaskQAsym(-1);
     if (!taskQASym) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQASym cannot run for this train conditions - EXCLUDED");
   }



   if(iPWG4TmpSourceSara){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskEta.C");
     AliAnalysisTaskEta *taskEta = AddTaskEta();
     if (!taskEta) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskEta cannot run for this train conditions - EXCLUDED");
   }


   if(iPWG4JetServices){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetServices.C");
     AliAnalysisTaskJetServices *taskjetServ = 0;
     taskjetServ = AddTaskJetServices("/Users/kleinb/Dropbox/SharedJets/Christian/Files/PWG4_JetTasksOutput_110818a.root");
     if (!taskjetServ) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetServices cannot run for this train conditions - EXCLUDED");
     if(kGridRunRange[0]>0)taskjetServ->SetRunRange(kGridRunRange[0],kGridRunRange[1]);
     else taskjetServ->SetRunRange(110000,160000);
     taskjetServ->SetMCData(kIsMC);
     taskjetServ->SetUsePhysicsSelection((Bool_t)iPhysicsSelection);
     taskjetServ->SetPhysicsSelectionFlag(iPhysicsSelectionFlag); // 
     taskjetServ->SetNonStdFile(kDeltaAODJetName.Data());
     taskjetServ->SetTrackEtaWindow(fTrackEtaWindow);
       taskjetServ->SetZVertexCut(fVertexWindow);
     taskjetServ->SetFilterMask(kHighPtFilterMask);

     if(kIsPbPb)taskjetServ->SetCollisionType(AliAnalysisTaskJetServices::kPbPb);
     else taskjetServ->SetCollisionType(AliAnalysisTaskJetServices::kPP);
     if(kIsPbPb){
       if(kDeltaAODJetName.Length()>0&&kFilterAOD)taskjetServ->SetFilterAODCollisions(kTRUE);
       //       else if(iAODanalysis)taskjetServ->SetFilterAODCollisions(kTRUE);
       //       taskjetServ->SetDebugLevel(3);
     }
     if(iAODanalysis){
       //  
       taskjetServ->SetAODInput(kTRUE);
     }
   }

   if(iPWG4JetSpectrum){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetSpectrum2.C");
     AliAnalysisTaskJetSpectrum2 *taskjetSpectrum = 0;
     TString bkgClusters = kDefaultJetBackgroundBranch.Data(); 
     bkgClusters.ReplaceAll(Form("%s_",AliAODJetEventBackground::StdBranchName()),"");
     Printf("############# Possible jet branches ###################");
     for(int iJF = 0;iJF < kJetListSpectrum.GetSize();iJF++){
       TObjString *objStr = (TObjString*)kJetListSpectrum.At(iJF);
       Printf("%3d: %s",iJF+1,objStr->GetString().Data());
     }

     //     Printf("Type q to exit");if(getchar()=='q')return;
     if(iPWG4JetSpectrum&1){
       if(kIsPbPb){
	 for(int iJF = 0;iJF < kJetListSpectrum.GetSize();iJF++){
	   Long64_t value = kJetMapSpectrum(iJF+1);
	   Printf("iJF: %d value: %d", iJF+1,(Int_t)value);
	   if(value==0)continue;
	   TObjString *objStr = (TObjString*)kJetListSpectrum.At(iJF);
	   TString bName1 = objStr->GetString().Data();
	   TString bName2[3];
	   
	   TObjString *objStrBkg = (TObjString*)kJetListSpectrum.At(iJF);
	   TString bBkgName1(""); 
	   TString bBkgName2[3];
	   
	   Long64_t valueBkg1 = kJetBackMapSpectrum(iJF+1);
	   if(valueBkg1>0){
	     TObjString *objStrBkg = (TObjString*)kJetListSpectrum.At(valueBkg1-1);
	     bBkgName1 = objStrBkg->GetString().Data();
	   }
	   Int_t iPartner = 0;
	   if(value>0){
	     Int_t iJF2 = -1; 
	     for(int i = 0;i<3;i++){
	       iJF2 = value/kJetMapOffset[i]-1;
	       value = value%kJetMapOffset[i];
	       Printf("%d %d", iJF2+1,(Int_t)value);
	       if(iJF2>=0&&iJF2<kJetListSpectrum.GetSize()){
		 TObjString *objStr2 = (TObjString*)kJetListSpectrum.At(iJF2);
		 bName2[iPartner] = objStr2->GetString().Data();
		 Long64_t valueBkg2 = kJetBackMapSpectrum(iJF2+1);
		 if(valueBkg2>0){
		   TObjString *objStrBkg2 = (TObjString*)kJetListSpectrum.At(valueBkg2-1);
		   bBkgName2[iPartner] = objStrBkg2->GetString().Data();
		 }
		 iPartner++;
	       }
	     }
	   }
	   
	   
	   // loop over all centralities
	   for(int ic = 0;ic<5;ic++){
	     if(ic!=0)continue;
	     Bool_t bDone = kFALSE;
	     for(int i = 0;i<TMath::Max(iPartner,1);i++){
	       if(bName2[i].Length()){
		 taskjetSpectrum = AddTaskJetSpectrum2(bName1.Data(),bName2[i].Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,ic);
		 bDone = kTRUE; 
	       }
	       else{
		 if(!bDone){
		   taskjetSpectrum = AddTaskJetSpectrum2(bName1.Data(),bName2[i].Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,ic);
		   bDone = kTRUE; 
		 }
		 else{
		   continue;
		 }
	       }
	       Printf("%s/%s %s/%s",bName1.Data(),bBkgName1.Data(),bName2[i].Data(),bBkgName2[i].Data());

	       // setting all the other things...
	       taskjetSpectrum->SetTrackEtaWindow(fTrackEtaWindow);
	       taskjetSpectrum->SetJetEtaWindow(fJetEtaWindow);

	       // negative values do not fill the track histos
	       taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetGenFull,-1);
	       taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetGen,-1);	      
	       if(!bName1.Contains("ANTIKT")||bName1.Contains("Cone")){
		 taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetRecFull,-1);
		 taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetRec,-1);	      
	       }


	       if(bName2[i].Length()==0){
		 taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetGenFull,0);
		 taskjetSpectrum->SetFlagJetType(AliAnalysisTaskJetSpectrum2::kJetGen,0);
	       }
	       if(iAODanalysis==1){
		 taskjetSpectrum->SetAODJetInput(kTRUE);
		 taskjetSpectrum->SetAODTrackInput(kTRUE);
	       }
	       else if (iAODanalysis==2){
		 taskjetSpectrum->SetAODTrackInput(kTRUE);
		 taskjetSpectrum->SetAODJetInput(kFALSE);
	       }
	       //	       taskjetSpectrum->SetDebugLevel(11);
	       taskjetSpectrum->SetBranchBkgRec(bBkgName1.Data());
	       taskjetSpectrum->SetBranchBkgGen(bBkgName2[i].Data());
	     } 
	   }
	 }
       }
       else{ // ! PbPb
	 Int_t i = 0;

	 taskjetSpectrum = AddTaskJetSpectrum2(kDefaultJetBranch.Data(),
					       "",kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,i);
	 

	 
	 if(kDefaultJetBranchMC.Length()){
	   taskjetSpectrum = AddTaskJetSpectrum2(kDefaultJetBranch.Data(),
					       kDefaultJetBranchMC.Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,i);
	   //	   taskjetSpectrum->SetMinJetPt(10);
	   taskjetSpectrum->SetTrackEtaWindow(fTrackEtaWindow);
	   taskjetSpectrum->SetJetEtaWindow(fJetEtaWindow);
	   taskjetSpectrum->SetBranchBkgRec(bkgClusters.Data());
	   taskjetSpectrum->SetBranchBkgGen(bkgClusters.Data());
	   if(iAODanalysis)SetAODInput(taskjetSpectrum);
	 }
	 if(kDefaultJetBranchMC2.Length()){
	   taskjetSpectrum = AddTaskJetSpectrum2(kDefaultJetBranch.Data(),
						 kDefaultJetBranchMC2.Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,i);
	   //	   taskjetSpectrum->SetMinJetPt(10);
	   taskjetSpectrum->SetTrackEtaWindow(fTrackEtaWindow);
	   taskjetSpectrum->SetJetEtaWindow(fJetEtaWindow);
	   taskjetSpectrum->SetBranchBkgRec(bkgClusters.Data());
	   taskjetSpectrum->SetBranchBkgGen(bkgClusters.Data());
	   if(iAODanalysis)SetAODInput(taskjetSpectrum);
	 }

	 /*
	 TString tmp2(kDefaultJetBranch.Data());
	 tmp2.ReplaceAll(Form(kJetSubtractMask2.Data(),0),Form(kJetSubtractMask2.Data(),2));
	 taskjetSpectrum = AddTaskJetSpectrum2(tmp2.Data(),kDefaultJetBranch.Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,i);
	 //	 taskjetSpectrum->SetDebugLevel(3);
	 //	   taskjetSpectrum->SetMinJetPt(10);
	 taskjetSpectrum->SetBranchBkgRec(bkgClusters.Data());
	 taskjetSpectrum->SetBranchBkgGen(bkgClusters.Data());
	 taskjetSpectrum->SetTrackEtaWindow(0.8);
	 taskjetSpectrum->SetJetEtaWindow(0.4);

	 // check the old subtracted vs. the new subtracted
	 TString tmp3(kDefaultJetBranch.Data());
	 tmp3.ReplaceAll(Form(kJetSubtractMask2.Data(),0),Form(kJetSubtractMask2.Data(),3));
	 taskjetSpectrum = AddTaskJetSpectrum2(tmp3.Data(),kDefaultJetBranch.Data(),kDeltaAODJetName.Data(),kHighPtFilterMask,AliVEvent::kMB,0,i);
	 //	 taskjetSpectrum->SetDebugLevel(3);
	 //	   taskjetSpectrum->SetMinJetPt(10);
	 taskjetSpectrum->SetBranchBkgRec(bkgClusters.Data());
	 taskjetSpectrum->SetBranchBkgGen(bkgClusters.Data());
	 taskjetSpectrum->SetTrackEtaWindow(0.8);
	 taskjetSpectrum->SetJetEtaWindow(0.4);
	 if(iAODanalysis)SetAODInput(taskjetSpectrum);
	 */
       }
       if (!taskjetSpectrum) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetSpectrum2 cannot run for this train conditions - EXCLUDED");
     }
   }
   AliAnalysisManager::SetCommonFileName("PWG4_Fragmentation.root");
   if(iPWG4Fragmentation){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskFragmentationFunction.C");
       AliAnalysisTaskFragmentationFunction *taskFrag = 0;
       if(kUseAODMC){

       }
       else{


	 for(int ic = 1;ic < 5;ic++){
	   // Anti-kT
	   taskFrag = AddTaskFragmentationFunction(1<<23,kHighPtFilterMask, ic);
	   if(kDeltaAODJetName.Length()>0)taskFrag->SetNonStdFile(kDeltaAODJetName.Data());
	   if(iAODanalysis==2)taskFrag->UseAODInputJets(kFALSE);

	   // UA1
	   taskFrag = AddTaskFragmentationFunction(1<<0,kHighPtFilterMask, ic); 
	   if(kDeltaAODJetName.Length()>0)taskFrag->SetNonStdFile(kDeltaAODJetName.Data());
	   if(iAODanalysis==2)taskFrag->UseAODInputJets(kFALSE);

         // SISCONE 
	   if(ic==1){
	     /*
	       taskFrag = AddTaskFragmentationFunction(1<<28,kHighPtFilterMask, ic);
	       taskFrag = AddTaskFragmentationFunction(1<<29,kHighPtFilterMask, ic);
	       taskFrag = AddTaskFragmentationFunction(1<<30,kHighPtFilterMask, ic);
	     */



	     // Anti-kT B2 - B3
	     if(!iAODanalysis==1){
	       taskFrag = AddTaskFragmentationFunction(1<<26,kHighPtFilterMask, ic);
	       if(kDeltaAODJetName.Length()>0)taskFrag->SetNonStdFile(kDeltaAODJetName.Data());
	       if(iAODanalysis==2)taskFrag->UseAODInputJets(kFALSE);
	       
	       taskFrag = AddTaskFragmentationFunction(1<<27,kHighPtFilterMask, ic);
	       if(kDeltaAODJetName.Length()>0)taskFrag->SetNonStdFile(kDeltaAODJetName.Data());
	       if(iAODanalysis==2)taskFrag->UseAODInputJets(kFALSE);
	     }
	   }

	 } 


       }
       if (!taskFrag) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskFragmentationFunction cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4JetChem){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetChem.C");
     AliAnalysisTask *taskChem = AddTaskJetChem(kHighPtFilterMask);
     if (!taskChem) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetChem cannot run for this train conditions - EXCLUDED");
   }

   if (iPWG4JetResponse) {
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetResponse.C");
     AliAnalysisTask *taskJetResponse = 0;

     if(iJETAN){
        taskJetResponse = AddTaskJetResponse("jets", "FASTJET", 0.4, kHighPtFilterMask, 0.15, 0);
        taskJetResponse = AddTaskJetResponse("jets", "FASTJET", 0.4, kHighPtFilterMask, 0.15, 1);

        taskJetResponse = AddTaskJetResponse("jets", "UA1", 0.4, kHighPtFilterMask, 1., 0);
        taskJetResponse = AddTaskJetResponse("jets", "UA1", 0.4, kHighPtFilterMask, 1., 2);
     }
     if(iPWG4Cluster){
        taskJetResponse = AddTaskJetResponse("clusters", "ANTIKT", 0.4, kHighPtFilterMask, 0.15, 0);
        taskJetResponse = AddTaskJetResponse("clusters", "ANTIKT", 0.4, kHighPtFilterMask, 0.15, 1);
     }
     if (!taskJetResponse) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetResponse cannot run for this train conditions - EXCLUDED");

   }

   if(iPWG4JCORRAN){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJCORRANTask.C");
     AliJCORRANTask* corran = AddTaskJCORRAN(kDeltaAODJCORRANName.Data(),0);
     if(!corran)::Warning("AnalysisTrainPWG4Jets", "AliJCORRANTask cannot run for this train conditions - EXCLUDED");
     else{
       if(kDeltaAODJCORRANName.Length()>0)mgr->RegisterExtraFile(kDeltaAODJCORRANName.Data()); 
     }
   }

   if(iPWG4UE){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskUE.C");
     AliAnalysisTaskUE *taskUE = 0;
     if(iPWG4UE&1)taskUE = AddTaskUE(); 
     if(iPWG4UE&2){
       taskUE  =AddTaskUE("jetsAOD_CDF04","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       //       taskUE  =AddTaskUE("jetsAOD_CDF07","CDF", "LJ", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_SISCONE04","CDF", "LJ", "TRANSV","MSP");
       //       taskUE  =AddTaskUE("jetsAOD_SISCONE07","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       taskUE  =AddTaskUE("jetsAOD_ICDF","CDF","LJ","TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_FASTKT04","CDF", "LJ", "TRANSV","MSP");
       //       taskUE  =AddTaskUE("jetsAOD_FASTKT07","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       taskUE  =AddTaskUE("jetsAOD_NONE","CDF", "MP_eta05", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_NONE","CDF", "MP_eta09", "TRANSV","MSP");
     }

     if (!taskUE) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskUE cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4LeadingUE){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskLeadingTrackUE.C");
     AliAnalysisTaskLeadingTrackUE *taskLeadingUE = AddTaskLeadingTrackUE(kUseMC);
     if (!taskLeadingUE) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTasLeadingTrackkUE cannot run for this train conditions - EXCLUDED");
     //     taskLeadingUE->SetFilterBit(64);
   }


   if(iPWG4CorrectionsUE){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskCorrectionsUE.C");
     AliAnalysisTaskCorrectionsUE *taskCorrectionsUE = 0;
     if(iPWG4CorrectionsUE&1)taskCorrectionsUE = AddTaskCorrectionsUE("jetsAOD_NONE","CDF","MP_eta05","TRANSV","MSP",kFALSE);
     if (!taskCorrectionsUE) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskCorrectionsUE cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4ThreeJets){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskThreeJets.C");
     AliAnalysisTaskThreeJets *taskThree = AddTaskThreeJets();
     if(!taskThree)::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskThreets cannot run for this train conditions - EXCLUDED");
   }
   if(iPWG4QGSep){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskQGSep.C");
     AliAnalysisTaskQGSep *taskQGSep = AddTaskQGSep(kUseMC,iAODanalysis);
     if(!taskQGSep)::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQGSep cannot run for this train conditions - EXCLUDED");
   }
  

   if(iPWG4Minijet){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskMinijet.C");
     AliAnalysisTaskMinijet *taskMini = AddTaskMinijet(-1,"esd",kUseMC,kGridDataSet);
     // if we ha highmult trigger add another task
     if(!taskMini)::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskMinjet cannot run for this train conditions - EXCLUDED");
   }

   AliAnalysisManager::SetCommonFileName("PWG4_HighPtQA.root");
   if(iPWG4PtQAMC){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtQAMC.C");
     AliPWG4HighPtQAMC *taskQAMC = 0;
     if(kUseMC){
       if(iPWG4PtQAMC&1){
	 taskQAMC = AddTaskPWG4HighPtQAMCAll(kGridDataSet.Data());
       }
     }
     if (!taskQAMC) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQAMC cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4PtTrackQA){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtTrackQA.C");
     if(iPWG4PtTrackQA&2)AddTaskPWG4HighPtTrackQAAll(kGridDataSet.Data(),kIsPbPb,iAODanalysis);
     else AddTaskPWG4HighPtTrackQAAllReduced(kGridDataSet.Data(),kIsPbPb,iAODanalysis);
   }

   if(iPWG4PtQATPC){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtQATPConly.C");
     AliPWG4HighPtQATPConly *taskQATPC = 0;
     if(iPWG4PtQATPC&1)taskQATPC = AddTaskPWG4HighPtQATPConly(kGridDataSet.Data(),1);
     if(iPWG4PtQATPC&2)taskQATPC = AddTaskPWG4HighPtQATPConly(kGridDataSet.Data(),2);

 if (!taskQATPC) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQATPC cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4Cosmics){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4CosmicCandidates.C");

     AliPWG4CosmicCandidates *taskPWG4CosmicCandidates = AddTaskPWG4CosmicCandidates(0);
     taskPWG4CosmicCandidates = AddTaskPWG4CosmicCandidates(1);

     if (!taskPWG4CosmicCandidates) ::Warning("AnalysisTrainPWG4Jets", "AddTaskPWG4CosmicCandidates cannot run for this train conditions - EXCLUDED");
   }


   if(iPWG4PtSpectra){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtSpectra.C");
     AddTaskPWG4HighPtSpectraAll(kGridDataSet.Data(),kIsPbPb,iAODanalysis);
   }
   
   if(iPWG4KMeans){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskKMeans.C");
     AliAnalysisTaskKMeans *taskKMeans = AddTaskKMeans();
     if (!taskKMeans) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskKMenans cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4PartCorr){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
     AliAnalysisTaskParticleCorrelation *taskpartcorrPHOS = AddTaskPartCorr("AOD", "PHOS",kFALSE,kIsMC);
     if (!taskpartcorrPHOS) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation PHOS cannot run for this train conditions - EXCLUDED");
     AliAnalysisTaskParticleCorrelation *taskpartcorrEMCAL = AddTaskPartCorr("AOD", "EMCAL",kFALSE,kIsMC);
     if (!taskpartcorrEMCAL) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation EMCAL cannot run for this train conditions - EXCLUDED");
     if(kDeltaAODPartCorrName.Length()>0)mgr->RegisterExtraFile(kDeltaAODPartCorrName.Data());
   } 

   if(iPWG4CaloQA){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
     AliAnalysisTaskParticleCorrelation *taskcaloQA =  AddTaskCalorimeterQA("ESD",kFALSE,kIsMC,kCaloQAOutputFileName.Data());
     if(!taskcaloQA)::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation QA cannot run - EXCLUDED");
     //  if(kCaloQAOutputFileName.Length()>0)mgr->RegisterExtraFile(kCaloQAOutputFileName.Data());
   } 

   if(iPWG4JetCorr){
     //     using namespace JetCorrelHD;
     TString cdir = gSystem->WorkingDirectory();
     gSystem->ChangeDirectory(gSystem->ExpandPathName("$ALICE_ROOT/PWG4/macros/"));
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetCorrel.C");
     AliAnalysisTaskJetCorrel *taskjetcorr = AddTaskJetCorrel();
     gSystem->ChangeDirectory(cdir);
     if (!taskjetcorr) ::Warning("AnalysisTrainNew", "AliAnalysisTaskJetCorrel  cannot run for this train conditions - EXCLUDED");
   } 

   if(iPWG4Tagged){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskTaggedPhotons.C");
     AliAnalysisTaskTaggedPhotons * taskTagged = AddTaskTaggedPhotons(kFALSE); // EMCAL
     taskTagged = AddTaskTaggedPhotons(kTRUE); // PHOS 
     if (!taskTagged) ::Warning("AnalysisTrainNew", "AliAnalysisTaskTaggedPhotons  cannot run for this train conditions - EXCLUDED");     
   }
   if (iPWG4omega3pi) {
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskomega3pi.C");
     AliAnalysisTaskOmegaPi0PiPi *taskomega3pi = AddTaskomega3pi();
     if (!taskomega3pi) ::Warning("AnalysisTrainNew", "AliAnalysisTaskomega3pi cannot run\
 for these train conditions - EXCLUDED");
   }

   // PWG4 gamma conversion analysis
   if (iPWG4GammaConv) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskGammaConversion.C");
      TString cdir = gSystem->WorkingDirectory();
      gSystem->ChangeDirectory(gSystem->ExpandPathName("$ALICE_ROOT/PWG4/macros/"));
      //      TString gcArguments = "-run-on-train -run-jet -run-chic -run-neutralmeson -run-cf";
      //      TString gcArguments = "-run-on-train -run-jet -run-neutralmeson -run-cf -use-own-xyz";
      //      TString gcArguments = "-run-on-train -run-jet -run-neutralmeson -run-cf -use-own-xyz";
      TString gcArguments = "-run-on-train -run-jet -run-omega-meson -use-own-xyz -run-neutralmeson -no-aod";
      //      TString kGCAnalysisCutSelectionId="9003562000100310";
      //      TString kGCAnalysisCutSelectionId="9003562000100312";
      //      gcArguments.Append(Form(" -set-cut-selection %s ",kGCAnalysisCutSelectionId.Data()));
      if(!kIsMC)gcArguments += " -mc-off";
      AliAnalysisTaskGammaConversion * taskGammaConversion = AddTaskGammaConversion(gcArguments,mgr->GetCommonInputContainer());
      gSystem->ChangeDirectory(cdir);
      taskGammaConversion->SelectCollisionCandidates();
      if (!taskGammaConversion) ::Warning("AnalysisTrainNew", "AliAnalysisTaskGammaConversion cannot run for these train conditions - EXCLUDED");
   }   

   if (iPWG4CaloConv) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskCaloConv.C");
      AliAnalysisTaskCaloConv * taskCaloConv = AddTaskCaloConv();
      if (!taskCaloConv) ::Warning("AnalysisTrainNew", "AliAnalysisTaskCaloConv cannot run for these train conditions - EXCLUDED");
   }   



   //==========================================================================
   // FOR THE REST OF THE TASKS THE MACRO AddTaskXXX() is not yet implemented/
   // Run the analysis
   //    
   if (kPluginUse) {
      AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
      AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
   }
   
   if (mgr->InitAnalysis()) {
     mgr->PrintStatus();
     // if (kSaveTrain || strlen(config_file)) gSystem->ChangeDirectory(kTrainName);
     if (!strcmp(plugin_mode,"submit")&&smode=="GRID"){
       TString alien_workdir = gGrid->GetHomeDirectory();
       if (iAODanalysis) alien_workdir += "analysisAOD";
       else              alien_workdir += "analysisESD";
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s%s",kGridDataSet.Data(),kGridExtraAliendirLevel.Data());
       AliAnalysisAlien *gridhandler = (AliAnalysisAlien*)mgr->GetGridHandler();
       printf("=== AnalysisTrainPWG4Jets:: Registering jdl in the work directory alien://%s/%s, should be done by the manager! ===\n",
	      alien_workdir.Data(),gridhandler->GetGridOutputDir());

       TString dest;
       dest = Form("%s/%s/%s.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 Printf("%s exist on grid removing...",dest.Data());
	 gGrid->Rm(dest.Data());
       }
       Printf("%s copy ...",dest.Data());
       TFile::Cp(Form("file:%s.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));


       TString dest;
       dest = Form("%s/%s/%s_merge.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 Printf("%s exist on grid removing...",dest.Data());
	 gGrid->Rm(dest.Data());
       }
       Printf("%s copy ...",dest.Data());
       TFile::Cp(Form("file:%s_merge.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));
       
       dest = Form("%s/%s/%s_merge_final.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 Printf("%s exist on grid removing...",dest.Data());
	 gGrid->Rm(dest.Data());
       }
       Printf("%s copy ...",dest.Data());
       TFile::Cp(Form("file:%s_merge.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));
     }
     AliLog::SetGlobalLogLevel(AliLog::kError);
     if((kUseSysInfo>0&&smode=="LOCAL")||!strcmp(plugin_mode, "test")){
       TFile * fM = TFile::Open("manager_local.root","RECREATE");
       mgr->Write();
       fM->Close();
     }
     
     // grmpf, aliroot error handler overwrites root
     gErrorIgnoreLevel = kErrorIgnoreLevel;
     if(gErrorIgnoreLevel>3000) AliLog::SetGlobalLogLevel(AliLog::kFatal);
     StartAnalysis(smode, chain);
     if((kUseSysInfo>0&&smode=="LOCAL")||!strcmp(plugin_mode, "test")){
       for(int i = 0;i < mgr->GetTopTasks()->GetEntries();i++){
	 mgr->ProfileTask(i);
       }
     }
     if (!strcmp(plugin_mode, "offline")&&smode=="GRID"){
       // Offline mode path files
       //	PatchJDL();
       //       PatchAnalysisMacro();
     }

     if (kSaveTrain && smode=="GRID") {
       AliAnalysisAlien *gridhandler = (AliAnalysisAlien*)mgr->GetGridHandler();
       TString alien_workdir = gGrid->GetHomeDirectory();
       if (iAODanalysis) alien_workdir += "analysisAOD";
       else              alien_workdir += "analysisESD";
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s%s",kGridDataSet.Data(),kGridExtraAliendirLevel.Data());
       //     kGridOutdir = gridhandler->GetGridOutputDir();
       printf("=== Registering ConfigTrain.C in the work directory <%s> ===\n",
                alien_workdir.Data());
       if (AliAnalysisAlien::FileExists(Form("%s/%sConfig.C", alien_workdir.Data(), kTrainName.Data())))
	 gGrid->Rm(Form("%s/%sConfig.C", alien_workdir.Data(), kTrainName.Data()));
       if (strcmp(plugin_mode, "test"))
	 TFile::Cp(Form("file:%sConfig.C",kTrainName.Data()), Form("alien://%s/%sConfig.C", alien_workdir.Data(), kTrainName.Data()));
     }
   }
}

//______________________________________________________________________________
void StartAnalysis(const char *mode, TChain *chain) {
// Start analysis.
   Int_t imode = -1;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   switch (imode) {
      case 0:
         if (!chain) {
            ::Error("AnalysisTrainPWG4Jets.C::StartAnalysis", "Cannot create the chain");
            return;
         }   
         mgr->StartAnalysis(mode, chain,kNumberOfEvents);
         return;
      case 1:
         if (!kProofDataSet.Length()) {
            ::Error("AnalysisTrainPWG4Jets.C::StartAnalysis", "kProofDataSet is empty");
            return;
         }   
         mgr->StartAnalysis(mode, kProofDataSet, kProofEvents,kProofOffset);
         return;
      case 2:
         if (kPluginUse) {
            if (!mgr->GetGridHandler()) {
               ::Error("AnalysisTrainPWG4Jets.C::StartAnalysis", "Grid plugin not initialized");
               return;
            }   
            mgr->StartAnalysis("grid",chain,kNumberOfEvents);
         } else {
            if (!chain) {
               ::Error("AnalysisTrainPWG4Jets.C::StartAnalysis", "Cannot create the chain");
               return;
            }   
	    //            mgr->StartAnalysis(mode, chain);
            mgr->StartAnalysis(mode, chain,kNumberOfEvents);
         }   
         return;
   }      
}          
    
//______________________________________________________________________________
void CheckModuleFlags(const char *mode) {
// Checks selected modules and insure compatibility
   Int_t imode = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;


   if (kUseCPAR) {
     kPluginAliRootVersion    = ""; // NO aliroot if we use CPAR
   }

   if (imode==1) {
      if (!kUsePAR) {
         ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PAR files enabled due to PROOF analysis");
         kUsePAR = kTRUE;
      }   
   }  
   if (imode != 2) {
      ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "AliEn plugin disabled since not in GRID mode");
      kPluginUse = kFALSE; 
   }   


   if(!kIsMC){
     // switch off anthin related to MC
     kUseMC = 0;
     kUseAODMC = 0;
     kUseTR = kFALSE;
   }

   // Decide if we have PbPb
   if(kGridDataSet.CompareTo("LHC10h")==0||kGridDataSet.Contains("LHC10h")) {
     Printf("Using settings for Pb+Pb");
     kIsPbPb = true;
   }

   

   if (iAODanalysis) {
   // AOD analysis
      if (kUseMC)
         ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "MC usage disabled in analysis on AOD's");
      if (kUseAODTags)
         ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "AOD tags usage disabled in analysis on AOD's");
      kUseMC = kFALSE;
      kUseTR = kFALSE;
      kUseAODTags = kFALSE;
      if (iESDfilter)
         ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "ESD filter disabled in analysis on AOD's");
      iESDfilter   = 0;
      if (iPhysicsSelection)
         ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "Physics Selection disabled in analysis on AOD's");
      iPhysicsSelection   = 0;
      if (!iAODhandler) {
         if (iJETAN) 
            ::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "JETAN disabled in analysis on AOD's without AOD handler");
         iJETAN = 0;
         iDIJETAN = 0;
      }
      // Disable tasks that do not work yet on AOD data
      if(iPWG4JCORRAN)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 JCORRAN disabled in analysis on AOD's");
      iPWG4JCORRAN = 0;
      if( iPWG4PtQAMC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtQAMC disabled in analysis on AOD's");
      iPWG4PtQAMC        = 0;
      if( iPWG4PtQATPC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtTPC disabled in analysis on AOD's");
      iPWG4PtQATPC        = 0;
      if( iPWG4Cosmics)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 Comics disabled in analysis on AOD's");
      iPWG4Cosmics        = 0;

      if(iPWG4KMeans)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4KMeans disabled on AOD's");
      iPWG4KMeans       = 0;
      if (iPWG4JetCorr)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4Jetcorr disabled on AOD's");
      iPWG4JetCorr = 0;
      if (iPWG4PartCorr)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4partcorr disabled on AOD's");
      iPWG4PartCorr = 0;
      if (iPWG4CaloQA)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4CaloQA disabled on AOD's");
      iPWG4CaloQA = 0;
      if (iPWG4Tagged)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4tagged disabled on AOD's");
      iPWG4Tagged = 0;
      if (iPWG4omega3pi)
	::Info("AnalysisTrainNew.C::CheckModuleFlags", "PWG4omega3pi disabled on AOD's");
      iPWG4omega3pi = 0;
      if(iPWG1QASym)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG1 QA Sym disabled in analysis on AOD's");
      if (iPWG4GammaConv)::Info("AnalysisPWG4Jets.C::CheckModuleFlags", "PWG4gammaconv disabled on AOD's");
      iPWG4GammaConv = 0;   
      iPWG1QASym     = 0;
      iCentralitySelection = 0;
   } else {   
   // ESD analysis

     if(kIsPbPb){
       iCentralitySelection = 1;
     }

     if (!kUseMC){
       kUseTR = kFALSE;
       
       if(kUseKinefilter)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 Kine Filter disabled in analysis without MC");
       kUseKinefilter = kFALSE;
       if( iPWG4PtQAMC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtQAMC disabled in analysis without MC");
       iPWG4PtQAMC        = 0;
       if( iPWG4CorrectionsUE)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 CorrectionsUE disabled in analysis without MC");
       iPWG4CorrectionsUE = 0;
     }
     if (iJETAN){
       iESDfilter=1;
     }
      if (!iESDfilter){
	kUseKinefilter = kFALSE;
	kUseMuonfilter = kFALSE;
      }
      if(!iJETAN){
	iPWG4JetSpectrum = iPWG4UE =  iPWG4CorrectionsUE = iPWG4ThreeJets = iPWG4QGSep = iDIJETAN = 0;
      }
   }
   iPWG4JetTasks = iPWG4JetServices||iPWG4JetSpectrum||iPWG4UE||iPWG4LeadingUE||iPWG4PtQAMC||iPWG4PtTrackQA||iPWG4PtSpectra||iPWG4PtQATPC||iPWG4Cosmics||iPWG4ThreeJets||iPWG4QGSep||iPWG4JetChem||iPWG4Minijet||iPWG4Fragmentation;
   iPWG4PartCorrLibs = iPWG4PartCorr||iPWG4Tagged||iPWG4CaloQA;
   iPWG4GammaConvLib = iPWG4GammaConv||iPWG4CaloConv;


   iEMCUtilLibs = iPWG4JetTasks||iPWG4PartCorrLibs||iPWG4JCORRAN||iPWG4GammaConvLib||iJETAN;
   iJETANLib = iPWG4JetTasks||iJETAN||iDIJETAN;

   if (iESDfilter) {iAODhandler=1;}
   if (kUseKinefilter && !kUseMC) kUseKinefilter = kFALSE;
   if (kUseAODTags && !iAODhandler) kUseAODTags = kFALSE;


   
}

//______________________________________________________________________________
Bool_t Connect(const char *mode) {
// Connect <username> to the back-end system.
   Int_t imode = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   TString username = gSystem->Getenv("alien_API_USER");
   switch (imode) {
      case 0:
         break;
      case 1:
         if  (!username.Length()) {
            ::Error(Form("AnalysisTrainPWG4Jets.C::Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         ::Info("AnalysisTrainPWG4Jets.C::Connect", "Connecting user <%s> to PROOF cluster <%s>", 
                username.Data(), kProofCluster.Data());
         gEnv->SetValue("XSec.GSI.DelegProxy", "2");
//         TProof::Open(Form("%s@%s:31093", username.Data(), kProofCluster.Data()));       
         TProof::Open(Form("%s@%s", username.Data(), kProofCluster.Data()));       
         if (!gProof) {
            if (strcmp(gSystem->Getenv("XrdSecGSISRVNAMES"), "lxfsrd0506.cern.ch"))
               ::Error(Form("AnalysisTrainPWG4Jets.C::Connect <%s>", mode), "Environment XrdSecGSISRVNAMES different from lxfsrd0506.cern.ch");
            return kFALSE;
         }
	 if(kProofClearPackages)gProof->ClearPackages();

	 if(kProofSaveToAlien){
	   TGrid::Connect("alien://");
	   if (gGrid) {
	     TString homedir = gGrid->GetHomeDirectory();
	     TString workdir = homedir + kTrainName;
	     if (!gGrid->Cd(workdir)) {
               gGrid->Cd(homedir);
               if (gGrid->Mkdir(workdir)) {
		 gGrid->Cd(kTrainName);
		 ::Info("AnalysisTrainPWG4Jets::Connect()", "Directory %s created", gGrid->Pwd());
               }
	     }
	     gGrid->Mkdir("proof_output");
	     gGrid->Cd("proof_output");
	     kProofOutdir = Form("alien://%s", gGrid->Pwd());
	   }   
	 }
         break;
      case 2:      
         if  (!username.Length()) {
            ::Error(Form("AnalysisTrainPWG4Jets.C::Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         if (kPluginUse && !gSystem->Getenv("alien_CLOSE_SE")) {
            ::Error(Form("AnalysisTrainPWG4Jets.C::Connect <%s>", mode), 
                           "When using the AliEn plugin it is preferable to define the \
                           variable alien_CLOSE_SE in your environment.");
            return kFALSE;
         }
         ::Info("AnalysisTrainPWG4Jets.C::Connect", "Connecting user <%s> to AliEn ...", 
                username.Data());
         TGrid::Connect("alien://");
         if (!gGrid || !gGrid->IsConnected()) return kFALSE;
         break;
      default:
         ::Error("AnalysisTrainPWG4Jets.C::Connect", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   ::Info("AnalysisTrainPWG4Jets.C::Connect","Connected in %s mode", mode);
   return kTRUE;
}

//______________________________________________________________________________
Bool_t LoadCommonLibraries(const char *mode)
{
   if (!gSystem->Getenv("ALICE_ROOT")) {
      ::Error("AnalysisTrainPWG4Jets.C", "Analysis train requires that ALICE_ROOT is set to pick up Configurations"); 
      return kFALSE;
   }   
   
   // Load common analysis libraries.
   Int_t imode = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   Bool_t success = kTRUE;
   // ROOT libraries
   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libMinuit.so");
   
   // Load framework classes. Par option ignored here.
   switch (imode) {
      case 0:
      case 2:
         if (kUseCPAR) {
            success &= LoadLibrary("STEERBase", mode, kTRUE);
            success &= LoadLibrary("ESD", mode, kTRUE);
            success &= LoadLibrary("AOD", mode, kTRUE);
            success &= LoadLibrary("ANALYSIS", mode, kTRUE);
            success &= LoadLibrary("OADB", mode, kTRUE);
            success &= LoadLibrary("ANALYSISalice", mode, kTRUE);
            success &= LoadLibrary("ROOTFILES", mode, kTRUE);
	    //            success &= LoadLibrary("EventMixing", mode,kTRUE);
            success &= LoadLibrary("CORRFW", mode, kTRUE);
         } else {   
            success &= LoadLibrary("libSTEERBase.so", mode);
            success &= LoadLibrary("libESD.so", mode);
            success &= LoadLibrary("libAOD.so", mode);
            success &= LoadLibrary("libANALYSIS.so", mode);
            success &= LoadLibrary("libOADB.so", mode, kTRUE);
            success &= LoadLibrary("libANALYSISalice.so", mode);
	    //            success &= LoadLibrary("libEventMixing.so", mode);
            success &= LoadLibrary("libCORRFW.so", mode);
            gROOT->ProcessLine(".include $ALICE_ROOT/include");
         }   
         break;
      case 1:
         Int_t ires = -1;
         if (kProofUseAFPAR && !gSystem->AccessPathName(kProofAFversion)) ires = gProof->UploadPackage(kProofAFversion);
         if (ires < 0) {
            success &= LoadLibrary("STEERBase", mode);
            success &= LoadLibrary("ESD", mode);
            success &= LoadLibrary("AOD", mode);
            success &= LoadLibrary("ANALYSIS", mode);
            success &= LoadLibrary("ANALYSISalice", mode);
            success &= LoadLibrary("EventMixing", mode);
            success &= LoadLibrary("CORRFW", mode);
         } else { 
            ires = gProof->EnablePackage(kProofAFversion);
            if (ires<0) success = kFALSE;
            success &= LoadLibrary("CORRFW", mode);
         }
         break;         
      default:
         ::Error("AnalysisTrainPWG4Jets.C::LoadCommonLibraries", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   if (success) {
      ::Info("AnalysisTrainPWG4Jets.C::LoadCommodLibraries", "Load common libraries:    SUCCESS");
      ::Info("AnalysisTrainPWG4Jets.C::LoadCommodLibraries", "Include path for Aclic compilation:\n%s",
              gSystem->GetIncludePath());
   } else {           
      ::Info("AnalysisTrainPWG4Jets.C::LoadCommodLibraries", "Load common libraries:    FAILED");
   }   
      
   return success;
}

//______________________________________________________________________________
Bool_t LoadAnalysisLibraries(const char *mode)
{
// Load common analysis libraries.
  Bool_t success = kTRUE;
  if (useTender) {
      if (!LoadLibrary("TENDER", mode, kTRUE) ||
          !LoadLibrary("TENDERSupplies", mode, kTRUE)) return kFALSE;
   }
   if (iESDfilter) {
     /*
      if (!LoadLibrary("PWG3base", mode, kTRUE) ||
          !LoadLibrary("PWG3muon", mode, kTRUE)) return kFALSE;
     */
   }   
   // JETAN
   if (iJETANLib) {
     // this part needs some rework in case we do not need the fastjed finders for processing
     if(iEMCUtilLibs){
       if (!LoadLibrary("EMCALUtils", mode, kTRUE) ||
	   !LoadLibrary("PHOSUtils", mode, kTRUE)) return kFALSE;
     }
     if (!LoadLibrary("JETAN", mode, kTRUE)) return kFALSE;
     if (!strcmp(mode, "PROOF")){
       gProof->Exec("gSystem->Load\(\"/afs/cern.ch/user/d/dperrino/public/libCGAL.so\"\)", kTRUE); 
       gProof->Exec("gSystem->Load\(\"/afs/cern.ch/user/d/dperrino/public/libfastjet.so\"\)", kTRUE); 
       // problem when loading siscone copiled with different gcc version??
       // gProof->Exec("gSystem->Load\(\"/afs/cern.ch/user/d/dperrino/public/libsiscone.so\"\)", kTRUE); 
       gProof->Exec("gSystem->Load\(\"/afs/cern.ch/user/d/dperrino/public/libSISConePlugin.so\"\)", kTRUE);      
     }
     if(!kUsePAR){ 
       if (!LoadLibrary("CGAL", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("fastjet", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("siscone", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("SISConePlugin", mode, kTRUE)) return kFALSE;
     }
     else{
       // par files plus FASTJET needs some extra work... need to change
       // the loading sequence in the auto generated .C file
       if (!LoadLibrary("libCGAL.so", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("libfastjet.so", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("libsiscone.so", mode, kTRUE)) return kFALSE;
       if (!LoadLibrary("libSISConePlugin.so", mode, kTRUE)) return kFALSE;
     }
     if (!LoadLibrary("FASTJETAN", mode, kTRUE)) return kFALSE;
   }
   if(iPWG4JetTasks){
     if (!LoadLibrary("PWG4JetTasks", mode, kTRUE)) return kFALSE;
   }

   if(iPWG1QASym){
     if (!LoadSource(Form("%s/PWG1/AliAnalysisTaskQASym.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }
   if(iPWG4TmpSourceSara){
     if(!kUsePAR)gSystem->AddIncludePath("-I$ALICE_ROOT/include/JetTasks"); // ugly hack!!
     if(!LoadSource(Form("%s/PWG4/JetTasks/AliAnalysisTaskEta.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }

   if (iPWG4PartCorrLibs) {   
      if (!LoadLibrary("PWG4PartCorrBase", mode, kTRUE) ||
          !LoadLibrary("PWG4PartCorrDep", mode, kTRUE)) return kFALSE;
   }
   if(iPWG4JCORRAN){
     // PWG4 particle correlations
     if(!LoadLibrary("PWG4JCORRAN",mode,kTRUE))return kFALSE;
   }
   if (iPWG4JetCorr) { 
     if (!LoadLibrary("PWG4JetCorrel", mode, kTRUE)) return kFALSE;
   }  
   if (iPWG4omega3pi) {
     if (!LoadLibrary("PWG4omega3pi", mode, kTRUE)) return kFALSE;
   }
   if (iPWG4GammaConvLib) {
      if (!LoadLibrary("PWG4GammaConv", mode, kTRUE)) return kFALSE;
   }      

   ::Info("AnalysisTrainPWG4Jets.C::LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
   return kTRUE;
}

//______________________________________________________________________________
Bool_t LoadLibrary(const char *module, const char *mode, Bool_t rec=kFALSE)
{
// Load a module library in a given mode. Reports success.
   Int_t imode = -1;
   Int_t result;
   TString smodule(module);
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   TString mod(module);
   if (!mod.Length()) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadLibrary", "Empty module name");
      return kFALSE;
   }   
   // If a library is specified, just load it
   if (smodule.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      result = gSystem->Load(mod);
      if (result < 0) {
         ::Error("AnalysisTrainPWG4Jets.C::LoadLibrary", "Could not load library %s", module);
         return kFALSE;
      }
      if (rec) anaLibs += Form("%s.so ",mod.Data()); 
      if (rec) anaLibsExtra += Form("%s.so ",mod.Data()); 
      return kTRUE;
   } 
   // Check if the library is already loaded
   if (strlen(gSystem->GetLibraries(Form("%s.so", module), "", kFALSE)) > 0)
      return kTRUE;    
   switch (imode) {
      case 0:
      case 2:
         if (kUsePAR) {
            result = SetupPar(module);
            if (rec) anaPars += Form("%s.par ", module);
         } else {
            result = gSystem->Load(Form("lib%s.so", module));
            if (rec) anaLibs += Form("lib%s.so ", module);
         }   
         break;
      case 1:
	if(!gSystem->AccessPathName(module)){
	  ::Info("AnalysisTrainPWG4Jets.C::LoadLibrary", "Removing directory %s",module);
	  gSystem->Exec(Form("rm -rf %s",module));
	}
         result = gProof->UploadPackage(module);
         if (result<0) {
            result = gProof->UploadPackage(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", module)));
            if (result<0) {
               ::Error("AnalysisTrainPWG4Jets.C::LoadLibrary", "Could not find module %s.par in current directory nor in $ALICE_ROOT", module);
               return kFALSE;
            }
         }   
         result = gProof->EnablePackage(module);
         break;
      default:
         return kFALSE;
   }         
   if (result < 0) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadLibrary", "Could not load module %s", module);
      return kFALSE;
   }
   return kTRUE;
}           



//______________________________________________________________________________
Bool_t LoadSource(const char *source, const char *mode, Bool_t rec=kFALSE)
{
// Load a module library in a given mode. Reports success.
   Int_t imode = -1;
   Int_t result = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   TString ssource(source);
   TString basename = gSystem->BaseName(ssource.Data());
   if (!ssource.Length()) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadSource", "Empty task name");
      return kFALSE;
   }   
   // we have a source code so compile it
   if (ssource.EndsWith(".cxx")) {
     // need to copy it here other wise the path is also used on grid...
     ssource.Remove(ssource.Index(".cxx"));
     basename.Remove(basename.Index(".cxx"));
     Printf("LoadSources:: Copying...  path %s{cxx,h}",ssource.Data());
     gSystem->Exec(Form("cp %s.cxx . ",ssource.Data()));
     gSystem->Exec(Form("cp %s.h . ",ssource.Data()));
     // Path to code
     // only needed for local compilation, in grid and proof mode 
     // the task headers are uploaded 
     //     path.Remove(path.Index(gSystem->BaseName(path.Data())));
     // Printf("LoadSources:: Including path %s",path.Data());
     //  if(path.Length()>0)gROOT->ProcessLine(Form(".include %s",path.Data()));
     Printf("LoadSources:: Loading...  path %s",basename.Data());
     switch (imode) {
     case 0:
       result = gROOT->LoadMacro(Form("%s.cxx++g",basename.Data()));
       break;
     case 1:
       result = gProof->LoadMacro(Form("%s.cxx++g",basename.Data()));
       break;
     case 2:
       result = gROOT->LoadMacro(Form("%s.cxx++g",basename.Data()));
       if (rec){
	 // what we want to compile
	 anaSources += Form("%s.cxx ",basename.Data()); 
	 // what we need as input...
	 anaLibs += Form("%s.cxx %s.h ",basename.Data(),basename.Data()); 
       }
       break;
     default:
       return kFALSE;
     }
   } 
   if (result < 0) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadSources", "Could not load source %s", source);
      return kFALSE;
   }
   return kTRUE;
}           


//______________________________________________________________________________
TChain *CreateChain(const char *mode, const char *plugin_mode)
{
// Create the input chain
   Int_t imode = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   TChain *chain = NULL;
   // Local chain
   switch (imode) {
      case 0:
         if (iAODanalysis) {
            if (!kLocalXMLDataset.Length()) {
               // Local AOD
               chain = new TChain("aodTree");
	       TString line;
	       ifstream in;
	       in.open(kLocalDataList.Data());
	       Int_t ic = 0;
	       while (in.good()) {
		 in >> line;
		 if (line.Length() == 0) continue;
		 Printf("%d adding %s",ic,line.Data());
		 chain->Add(line.Data());
		 ic++;
	       }       
            } else {
               // Interactive AOD
               chain = CreateChainSingle(kLocalXMLDataset, "aodTree");
            }
         } else {      
	   if (!kLocalXMLDataset.Length()) {
	     // Local ESD
	     chain = new TChain("esdTree");
	     TString line;
	     ifstream in;
	     in.open(kLocalDataList.Data());
	     while (in.good()) {
	       in >> line;
	       if (line.Length() == 0) continue;
	       cout << " line = " << line << endl;
	       chain->Add(line.Data());
	     }       
	   } else {
	     // Interactive ESD
               chain = CreateChainSingle(kLocalXMLDataset, "esdTree");
	   }   
         }
         break;
      case 1:
         break;
      case 2:
         if (kPluginUse) {
//            AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
//            AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
         } else {
            TString           treeName = "esdTree";
            if (iAODanalysis) treeName = "aodTree";
            chain = CreateChainSingle("wn.xml", treeName);
         }
         break;      
      default:   
   }
   if (chain && chain->GetNtrees()) return chain;
   return NULL;
}   

//______________________________________________________________________________
TChain* CreateChainSingle(const char* xmlfile, const char *treeName)
{
   printf("*******************************\n");
   printf("*** Getting the ESD Chain   ***\n");
   printf("*******************************\n");
   TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

   if (!myCollection) {
      ::Error("AnalysisTrainPWG4Jets.C::CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
      return NULL ;
   }

   TChain* chain = new TChain(treeName);
   myCollection->Reset() ;
   while ( myCollection->Next() ) chain->Add(myCollection->GetTURL("")) ;
   chain->ls();
   return chain;
}

//______________________________________________________________________________
Int_t SetupPar(char* pararchivename)
{
   if (!pararchivename || !strlen(pararchivename)) return -1;
   char processline[1024];
   if (gSystem->AccessPathName(Form("%s.par", pararchivename))) {
      if (!gSystem->AccessPathName(Form("%s/%s.par", gSystem->Getenv("ALICE_ROOT"),pararchivename))) {
         ::Info("AnalysisTrainPWG4Jets.C::SetupPar", "Getting %s.par from $ALICE_ROOT", pararchivename);
	TFile::Cp(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", pararchivename)), 
		  Form("%s.par",pararchivename));
      } else {
         ::Error("AnalysisTrainPWG4Jets.C::SetupPar", "Cannot find %s.par", pararchivename);
         return -1;
      }   
   }
   if (kPluginUse && kSaveTrain) gSystem->Exec(Form("ln -s ../%s.par %s",pararchivename, kTrainName.Data()));
   gSystem->Exec(Form("tar xvzf %s.par", pararchivename));

   TString ocwd = gSystem->WorkingDirectory();
   if (!gSystem->ChangeDirectory(pararchivename)) return -1;
	
   // check for BUILD.sh and execute
   if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");	    
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
         Error("runProcess","Cannot Build the PAR Archive! - Abort!");
         return -1;
      }
   }

	// check for SETUP.C and execute
	if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
	    printf("*******************************\n");
	    printf("*** Setup PAR archive       ***\n");
	    printf("*******************************\n");
	    gROOT->Macro("PROOF-INF/SETUP.C");
	}	
	if (!gSystem->ChangeDirectory(ocwd.Data())){
	  Error("Change back directory",Form("Cannot change to %s",ocwd.Data()));
	  return -1;
	}
   return 0;
}

//______________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(plugin_mode);
   plugin->SetCheckCopy(kFALSE); 
   if (kPluginUseProductionMode) plugin->SetProductionMode();
   plugin->SetJobTag(kJobTag);
   plugin->SetNtestFiles(1);
//   plugin->SetPreferedSE("ALICE::NIHAM::File");
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   //   plugin->SetAPIVersion("V1.0x");
//   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion(kPluginRootVersion);
   plugin->SetAliROOTVersion(kPluginAliRootVersion);
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir(kGridDatadir.Data());
// Set data search pattern
   if (iAODanalysis) plugin->SetDataPattern(Form(" %s/*/*AliAOD.root",kGridPassPattern.Data()));
   else              plugin->SetDataPattern(Form(" %s/*/*AliESDs.root",kGridPassPattern.Data()));
// ...then add run numbers to be considered
//   plugin->SetRunRange(kGridRunRange[0], kGridRunRange[1]);
   for (Int_t i=kGridRunRange[0]; i<=kGridRunRange[1]; i++) {
     Printf("AnalysisTrainPWG4Jets Adding run number %s", Form(kGridRunPattern.Data(),i));
     plugin->AddRunNumber(Form(kGridRunPattern.Data(),i));
   }   

   if(kGridLocalRunList.Length()>0){
     ifstream in1;
     in1.open(kGridLocalRunList.Data());
     int iRun = 0;
     int icount = 0;
     Int_t nRun = 0;
     // just use run numbers, negatives will be excluded
     while(in1>>iRun){
       if(iRun>=0){
	 if(iRun>=0&&nRun>=kGridOffsetRunFromList&&(nRun<kGridMaxRunsFromList)){
	   Printf("AnalysisTrainPWG4Jets Adding run number from File %d: %s",nRun,Form(kGridRunPattern.Data(),iRun));
	   plugin->AddRunNumber(Form(kGridRunPattern.Data(),iRun));

	 }
	 else{
	   Printf("AnalysisTrainPWG4Jets Skipping run number from File %d: %d",nRun, iRun);
	 }
	 nRun++;
       }
       else{
	 Printf("AnalysisTrainPWG4Jets Skipping run number from File %d: %d",nRun, iRun);
       }
     }
   }
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   TString alien_workdir = "";
   
   if (iAODanalysis)  alien_workdir += "analysisAOD";
   else               alien_workdir += "analysisESD";    
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s%s",kGridDataSet.Data(),kGridExtraAliendirLevel.Data());
   plugin->SetGridWorkingDir(alien_workdir.Data());

   // Declare alien output directory. Relative to working directory.
   if (!kGridOutdir.Length()) kGridOutdir = Form("output_%s",kTrainName.Data());
   plugin->SetGridOutputDir(kGridOutdir);

   // Add external packages
   plugin->AddExternalPackage("boost::v1_43_0");
   plugin->AddExternalPackage("cgal::v3.6");
   plugin->AddExternalPackage("fastjet::v2.4.2");


   // set extra libs before par file compilation
   anaLibs += kGridExtraFiles;
   anaLibs     = anaLibs.Strip();   
   Printf("anaLibs %s",anaLibs.Data());
   Printf("anaLibsExtra %s",anaLibsExtra.Data());

   if (anaLibs.Length())          plugin->SetAdditionalLibs(anaLibs.Data());
   if (anaLibsExtra.Length())     plugin->SetAdditionalRootLibs(anaLibsExtra.Data());

   TString ana_sources = "";
   TString ana_add = "";
   if (kUsePAR && anaPars.Length()) {
      printf("%s\n", anaPars.Data());
      TObjArray *arr;
      TObjString *objstr;
      arr = anaPars.Tokenize(" ");
      TIter next(arr);
      while ((objstr=(TObjString*)next())) plugin->EnablePackage(objstr->GetString());
      delete arr;
   } 
   
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   ana_sources = anaSources.Strip();
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.

   if (ana_sources.Length()) plugin->SetAnalysisSource(ana_sources);
   plugin->SetExecutableCommand(kPluginExecutableCommand.Data());  
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetMergeExcludes(kGridMergeExclude);
   plugin->SetMaxMergeFiles(kGridMaxMergeFiles);
   plugin->SetNrunsPerMaster(kGridRunsPerMaster);
   plugin->SetMergeViaJDL(kPluginMergeViaJDL);
   // Use fastread option
   plugin->SetFastReadOption(kPluginFastReadOption);
   // UseOverwrite mode
   plugin->SetOverwriteMode(kPluginOverwriteMode); 
   // Optionally define the files to be archived.
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:AliAOD.root,AOD.tag.root@ALICE::NIHAM::File");
   plugin->SetOutputToRunNo(kPluginOutputToRunNumber);     // write the output to subdirs named after run number
   
   // Put default output files to archive
   TString listhists = "";
   TString listaods  = "";
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TIter next(mgr->GetOutputs());
   AliAnalysisDataContainer *output;
   while ((output=(AliAnalysisDataContainer*)next())) {
      const char *filename = output->GetFileName();
      if (!(strcmp(filename, "default"))) {
	if (!mgr->GetOutputEventHandler()) continue;
         filename = mgr->GetOutputEventHandler()->GetOutputFileName();
         if (listaods.Length()) listaods += " ";
	 listaods += filename;
      } else {
	if(!listhists.Contains(filename)){
	  if (listhists.Length()) listhists += " ";
	  listhists += filename;
	}
      }
   }

   if (mgr->GetExtraFiles().Length()) {
     if (listaods.Length()) listaods += " ";
     listaods += mgr->GetExtraFiles();
   }

   // if we do not fill the aod we do not need to store it
   //   kGridMergeExclude = listaods;
   
   if(kSaveAOD>=0){
     TString outputFiles =  "";
     outputFiles += mgr->GetExtraFiles();
     if (listhists.Length()) outputFiles += " ";
     outputFiles += listhists;
     plugin->SetDefaultOutputs(kFALSE);
     Printf("%s:%d Starting with the files %s",(char*)__FILE__,__LINE__,outputFiles.Data());
     // remove
     // no harm done when we try to remove something that is not there :)
     if(!(kSaveAOD&(1<<0))){
       outputFiles.ReplaceAll("AliAOD.root ","");
       listaods.ReplaceAll("AliAOD.root ","");
     }
     if(!(kSaveAOD&(1<<1))){
       if(kDeltaAODJetName.Length())outputFiles.ReplaceAll(kDeltaAODJetName.Data(),"");
       if(kDeltaAODJetName.Length())listaods.ReplaceAll(kDeltaAODJetName.Data(),"");
       
     }
     if(!(kSaveAOD&(1<<2))){
       if(kDeltaAODPartCorrName.Length())outputFiles.ReplaceAll(kDeltaAODPartCorrName.Data(),"");
       if(kDeltaAODPartCorrName.Length())listaods.ReplaceAll(kDeltaAODPartCorrName.Data(),"");
     }
     if(!(kSaveAOD&(1<<3))){
       if(kDeltaAODJCORRANName.Length())outputFiles.ReplaceAll(kDeltaAODJCORRANName.Data(),"");
       if(kDeltaAODJCORRANName.Length())listaods.ReplaceAll(kDeltaAODJCORRANName.Data(),"");
     }
     
     // 
     plugin->SetDefaultOutputs(kFALSE);
     Printf("%s:%d Saving the files %s",(char*)__FILE__,__LINE__,outputFiles.Data());
     plugin->SetOutputFiles(outputFiles.Data());
   }

   TString outputArchive;
   outputArchive = Form("log_archive.zip:std*@%s","disk=1");
   listaods.ReplaceAll(" ", ",");
   listhists.ReplaceAll(" ", ",");
   if (listhists.Length()){
     outputArchive += " ";
     outputArchive += "root_archive.zip:";
     outputArchive += listhists;
     if (listaods.Length()){
       outputArchive += ",";
       outputArchive += listaods;
     }
     outputArchive += Form("@%s", kGridOutputStorages.Data());
   }
   else{
     
     if (listaods.Length()){
       // we have only aod'ish output
       outputArchive += " ";
       outputArchive += "root_archive.zip:";
       outputArchive += listaods;
       outputArchive += Form("@%s", kGridOutputStorages.Data());
     }
     else{
       // no other outputs than std..
      ::Fatal("AnalysisTrainPWG4Jets", "No task output !");
     }
   }

   plugin->SetDefaultOutputs(kFALSE);
   plugin->SetOutputArchive(outputArchive);


   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", kTrainName.Data()));
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(kGridFilesPerJob);
   // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   //   plugin->SetMaxInitFailed(5);
   // Optionally resubmit threshold.
   // plugin->SetMasterResubmitThreshold(90);
   // Optionally set time to live (default 30000 sec)
   plugin->SetTTL(80000); // 22h...
   // Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName(Form("%s.jdl", kTrainName.Data()));
   // Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable(Form("%s.sh", kTrainName.Data()));
   // Optionally modify job price (default 1)
   plugin->SetPrice(1);      
   // Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}

//______________________________________________________________________________
void WriteConfig()
{
// Write train configuration in a file. The file name has the format:
// train_[trainName]_ddMonthyyyy_time.C
   if (kUseDate) {
      gSystem->Exec("date +%d%b%Y_%Hh%M > date.tmp");
      ifstream fdate("date.tmp");
      if (!fdate.is_open()) {
         ::Error("AnalysisTrainPWG4Jets.C::Export","Could not generate file name");
         return;
      }
      const char date[64];
      fdate.getline(date,64);
      fdate.close();
      gSystem->Exec("rm date.tmp");
      kTrainName = Form("train_%s_%s", kTrainName.Data(), date);
   } else {
      kTrainName = Form("train_%s", kTrainName.Data());
   }   
   TString cdir = gSystem->WorkingDirectory();
   gSystem->MakeDirectory(kTrainName);
   gSystem->ChangeDirectory(kTrainName);
   ofstream out;
   out.open(Form("%sConfig.C",kTrainName.Data()), ios::out); 
   if (out.bad()) {
      ::Error("AnalysisTrainPWG4Jets.C::Export", "Cannot open ConfigTrain.C for writing");
      return;
   }
   out << "{" << endl;
   out << "   kTrainName      = " << "\"" << kTrainName.Data() << "\";" << endl;
   out << "   kProofCluster   = " << "\"" << kProofCluster.Data() << "\";" << endl;
   out << "   kProofUseAFPAR        = " << kProofUseAFPAR << ";" << endl;
   if (kProofUseAFPAR) 
      out << "   kProofAFversion       = " << kProofAFversion.Data() << ";" << endl;
   out << "   kProofDataSet   = " << "\"" << kProofDataSet.Data() << "\";" << endl;
   out << "   kPluginUse       = " << kPluginUse << ";" << endl;
   out << "   kUsePAR          = " << kUsePAR << ";" << endl;
   out << "   kUseCPAR         = " << kUseCPAR << ";" << endl;
   out << "   kPluginRootVersion    = " << "\"" << kPluginRootVersion.Data() << "\";" << endl;
   out << "   kPluginAliRootVersion = " << "\"" << kPluginAliRootVersion.Data() << "\";" << endl;
   out << "   kGridDatadir   = " << "\"" << kGridDatadir.Data() << "\";" << endl;
   if (!kGridOutdir.Length()) kGridOutdir = Form("output_%s",kTrainName.Data());
   out << "   kGridOutdir    = " << "\"" << kGridOutdir.Data() << "\";" << endl;
   out << "   kGridMaxMergeFiles   = " << kGridMaxMergeFiles << ";" << endl;
   out << "   kGridMergeExclude    = " << "\"" << kGridMergeExclude.Data() << "\";" << endl;
   out << "   kGridRunsPerMaster  = " << kGridRunsPerMaster << ";" << endl;
   out << "   kGridFilesPerJob    = " << kGridFilesPerJob << ";" << endl;
   out << "   kGridRunRange[0]    = " << kGridRunRange[0] << ";" << endl;
   out << "   kGridRunRange[1]    = " << kGridRunRange[1] << ";" << endl;
   out << "   kUseDebug          = " << kUseDebug << ";" << endl;
   out << "   kUseMC           = " << kUseMC << ";" << endl;
   out << "   kUseESDTags         = " << kUseESDTags << ";" << endl;
   out << "   kUseKinefilter      = " << kUseKinefilter << ";" << endl;
   out << "   kUseTR           = " << kUseTR << ";" << endl;
   out << "   kUseAODTags      = " << kUseAODTags << ";" << endl;
   out << "   kSaveTrain       = " << "kFALSE;" << endl << endl;
   out << "   // Analysis modules" << endl;
   out << "   iAODanalysis    = " << iAODanalysis << ";" << endl;
   out << "   iAODhandler     = " << iAODhandler << ";" << endl;
   out << "   iESDfilter      = " << iESDfilter << ";" << endl;
   out << "   iJETAN          = " << iJETAN << ";" << endl;
   out << "// Configuration fot the wagons" << endl;
   out << "}" << endl;
   ::Info("AnalysisTrainPWG4Jets.C::WriteConfig", "Train configuration wrote to file %s", Form("config_%s.C", kTrainName.Data()));
   gSystem->ChangeDirectory(cdir);
}   

//______________________________________________________________________________
Bool_t LoadConfig(const char *filename)
{
// Read train configuration from file
   if (gSystem->AccessPathName(filename)) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadConfig", "Config file name not found");
      return kFALSE;
   }   
   gROOT->ProcessLine(Form(".x %s", filename));
   ::Info("AnalysisTrainPWG4Jets.C::LoadConfig", "Train configuration loaded from file %s", filename);
   return kTRUE;
}

Bool_t PatchJDL(){
  Printf(">>> Patching JDL");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisAlien*    gridHandler = (AliAnalysisAlien*)mgr->GetGridHandler();
  TGridJDL *jdl = gridHandler->GetGridJDL();
  if(iJETAN)jdl->AddToPackages("fastjet","v2.4.0");
  gridHandler->WriteJDL(kFALSE);
  Printf("<<<  Patching JDL");
  return kTRUE;
}

Bool_t PatchAnalysisMacro(){
  Printf(">>> Patching AnalysisMacro");
  gSystem->Exec(Form("mv %s.C %s.C_tmp",kTrainName.Data(),kTrainName.Data()));

  ifstream in1; 
  in1.open(Form("%s.C_tmp", kTrainName.Data()));
  char cLine[250];
  TString st;
  while(in1.getline(cLine,250)){
    st += cLine;
    st += "\n";
  }
  Int_t index= -1;
  index = st.Index("gSystem->Load(\"libPhysics\");");
  index += strlen("gSystem->Load(\"libPhysics\");");
  /*
    TObjArray *arr;
    TObjString *objstr;
    arr = anaLibs.Tokenize(" ");
    TIter next(arr);

    add += "\n\n // added by CKB \n";
    while ((objstr=(TObjString*)next())){
      if(objstr->GetString().Contains("PWG3"))continue;
      if(objstr->GetString().EndsWith(".so"))add += Form("gSystem->Load(\"%s\");\n",objstr->GetString().Data());
    }
    delete arr; 
    */
    //    add += Form("AliLog::SetGlobalLogLevel(%d);\n",AliLog::GetGlobalLogLevel());
  TString add = "";

  if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
  add += "\n\n // added by CKB \n";
  if(kErrorIgnoreLevel>0)add += Form("gErrorIgnoreLevel = %d;\n",kErrorIgnoreLevel);
  add += "\n gSystem->AddIncludePath(\"./\"); \n";
  add += "\n gSystem->SetFPEMask(); \n";


  if(gGrid&&kPluginAliRootVersion.Length()==0){
    /*
    add += "\n // Dirty hack for TRD reference data \n";
    add += "\n gSystem->Setenv(\"ALICE_ROOT\",\"";
    add += Form("alien://%s/rootfiles/",gGrid->GetHomeDirectory());
    add += "\"); \n";
    */
  }

  add += "// BKC \n\n";
  st.Insert(index,add.Data());

  if(kUseDebug){
    //    st.Insert(index,"\n gROOT->ProcessLine(\".trace\"); // CKB \n");
  }

  if(kUseCPAR&&kPluginAliRootVersion.Length()==0){
    index = st.Index("gSystem->AddIncludePath(\"-I$"); // uncommen $ALICE_ROOT include for par files
    if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
    st.Insert(index,"// CKB comment out whehn no aliroot is provided \n //");
  }

  if(AliLog::GetGlobalLogLevel()==AliLog::kFatal){
    index = st.Index("AliLog::SetGlobal"); // ncomment setting of log level, do my own
    if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
    st.Insert(index,"AliLog::SetGlobalLogLevel(AliLog::kFatal);// CKB \n  // CKB comment out for own setting \n  //");
  }

  ofstream out;
  out.open(Form("%s.C", kTrainName.Data()));
  if (out.bad()) {
    return kFALSE;
  }
  out << st << endl;
  Printf("<<< Patching AnalysisMacro");





  Printf(">>> Patching Merge Macro");
  gSystem->Exec(Form("mv %s_merge.C %s_merge.C_tmp",kTrainName.Data(),kTrainName.Data()));

  ifstream in2; 
  in2.open(Form("%s_merge.C_tmp", kTrainName.Data()));
  TString st2;
  while(in2.getline(cLine,250)){
    st2 += cLine;
    st2 += "\n";
  }
  index = st2.Index("gSystem->Load(\"libPhysics\");");
  index += strlen("gSystem->Load(\"libPhysics\");");
  TString add2 = "";
  add2 += "\n gSystem->AddIncludePath(\"./\"); \n";
  if(gGrid&&kPluginAliRootVersion.Length()==0){
    /*
    add2 += "\n // Dirty hack for TRD reference data \n";
    add2 += "\n gSystem->Setenv(\"ALICE_ROOT\",\"";
    add2 += Form("alien://%s/rootfiles/",gGrid->GetHomeDirectory());
    add2 += "\"); \n";
    */
  }
  add2 += "// BKC \n\n";
  if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
  st2.Insert(index,add.Data());


  if(kUseDebug){
    //    st2.Insert(index,"\n gROOT->ProcessLine(\".trace\"); // CKB \n");
  }

  if(kUseCPAR&&kPluginAliRootVersion.Length()==0){
    index = st2.Index("gSystem->AddIncludePath(\"-I$"); // uncomment $ALICE_ROOT include for par files
    if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
    st2.Insert(index,"// CKB comment out whehn no aliroot is provided \n //");
  }

  // do not exclude the extra files from merign, this is done explicitly in this train script
  //  index = st2.Index("mergeExcludes +="); // uncommen $ALICE_ROOT include for par files
  //  if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
  // st2.Insert(index,"// CKB comment out, handled explicitly by the train macro \n //");


  ofstream out2;
  out2.open(Form("%s_merge.C", kTrainName.Data()));
  if (out2.bad()) {
    return kFALSE;
  }
  out2 << st2 << endl;
  Printf("<<< Patching Merging Macro");



  return kTRUE;

}
