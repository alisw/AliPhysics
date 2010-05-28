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
Bool_t      kFillAOD = kFALSE;  // switch of AOD filling for on the fly analysis

//== general input and output variables

Int_t       iAODanalysis       = 1;      // Analysis on input AOD's
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iESDfilter         = 0;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iPhysicsSelection  = 1;      // ESD to AOD filter (barrel + muon tracks)
Bool_t      kUseKinefilter     = kFALSE; // use Kinematics filter
Bool_t      kUseMuonfilter     = kFALSE; // use Kinematics filter
TString     kCommonOutputFileName = "PWG4_JetTasksOutput.root";


//== general process variables

// ### Other flags to steer the analysis
//==============================================================================
Bool_t      kSkipTerminate      = kTRUE; // Do not call Teminate
Bool_t      kUseDate            = kFALSE; // use date in train name
Bool_t      kUseDebug           = kTRUE; // activate debugging
Int_t       kUseSysInfo         = 0; // activate debugging
Long_t      kNumberOfEvents     = 1234567890; // number of events to process from the chain
Bool_t      kUseMC              = kTRUE;  // use MC info
Bool_t      kIsMC               = kTRUE;  // is MC info, if false it overwrites Use(AOD)MC
Bool_t      kUseAODMC           = kTRUE;  // use MC infA
Bool_t      kUseESDTags         = kFALSE; // use ESD tags for selection
Bool_t      kUseTR              = kFALSE;  // use track references
Bool_t      kUseAODTags         = kFALSE;  // use AOD tags
Bool_t      kSaveTrain          = kFALSE;  // save train configuration as: 


// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iJETAN             = 1;      // Jet analysis (PWG4) // 1 write standard 2 write non-standard jets, 3 wrtie both
Int_t       iDIJETAN           = 1;
Int_t       iJETANLib          = 1;
Int_t       iPWG1QASym         = 0;      // Eva's QA task compiled on the fly...
Int_t       iPWG4JetTasks      = 0;      // all jet tasks flag for lib laoding
Int_t       iPWG4JetServices   = 0;      // jet spectrum analysis
Int_t       iPWG4JetSpectrum   = 0;      // jet spectrum analysis
Int_t       iPWG4JCORRAN       = 0;      // JCORRAN module
Int_t       iPWG4UE            = 0;      // Underlying Event analysis
Int_t       iPWG4TmpSourceSara = 0;      // Underlying Event analysis not in svn
Int_t       iPWG4TmpSourceFrag = 0;      // Bastian's Fragmentation function not in svn
Int_t       iPWG4JetChem       = 0;      // Jet chemistry 
Int_t       iPWG4PtQAMC        = 0;      // Marta's QA tasks 
Int_t       iPWG4PtSpectra     = 0;      // Marta's QA tasks 
Int_t       iPWG4PtQATPC       = 0;      // Marta's QA tasks 
Int_t       iPWG4ThreeJets     = 0;      // Sona's thrust task
Int_t       iPWG4KMeans        = 0;      // Andreas' KMeans task 
Int_t       iPWG4Cluster       = 0;      // CKB cluster task 
Int_t       iPWG4PartCorrLibs  = 0;      // Gustavo's part corr analysis
Int_t       iPWG4PartCorr      = 0;      // Gustavo's part corr analysis
Int_t       iPWG4CaloQA        = 0;      // Gustavo's part corr analysis
Int_t       iPWG4JetCorr       = 0;     // Paul's jet corr analysis
Int_t       iPWG4Tagged        = 0;      // Gustavo's part corr analysis
Int_t       iPWG4omega3pi      = 0;      // Omega to 3 pi analysis (PWG4) 
Int_t       iPWG4GammaConv     = 0;      // Gamma Conversio
Int_t       kHighPtFilterMask  = 16;     // change depending on the used AOD Filter
TString     kDeltaAODJetName   = "AliAOD.Jets.root";     
TString     kDeltaAODPartCorrName   = "deltaAODPartCorr.root";     


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
TString     kPluginRootVersion       = "v5-26-00b-6";  // *CHANGE ME IF MORE RECENT IN GRID*
TString     kPluginAliRootVersion    = "v4-19-13-AN";  // *CHANGE ME IF MORE RECENT IN GRID*                                          
Bool_t      kPluginMergeViaJDL       = kTRUE;  // merge via JDL
Bool_t      kPluginFastReadOption   = kFALSE;  // use xrootd tweaks
Bool_t      kPluginOverwriteMode    = kTRUE;  // overwrite existing collections
Int_t       kPluginOutputToRunNumber = 1;     // write the output to subdirs named after run number
// TString kPluginExecutableCommand = "root -b -q";
TString     kPluginExecutableCommand = "source /Users/kleinb/setup_32bit_aliroot_trunk_clean_root_trunk.sh; alienroot -b -q ";

// == grid plugin input and output variables
TString     kGridDatadir      = "/alice/sim/PDC_08b/LHC09a1/AOD/";
TString     kGridLocalRunList = "";
TString     kGridOutdir       = ""; // AliEn output directory. If blank will become output_<kTrainName>
TString     kGridDataSet      = ""; // sub working directory not to confuse different run xmls 
Int_t       kGridRunRange[2]       =  {0, -1}; // Set the run range
TString     kGridRunPattern        = "%03d"; // important for leading zeroes!!
TString     kGridPassPattern       = "";
TString     kGridExtraFiles        = ""; // files that will be added to the input list in the JDL...
Int_t       kGridMaxMergeFiles      = 50; // Number of files merged in a chunkgridrunragn
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



// Temporaries.
TString anaPars = "";
TString anaLibs = "";
TString anaLibsExtra = "";
TString anaSources = "";
// Function signatures
class AliAnalysisAlien;

//______________________________________________________________________________
void AnalysisTrainPWG4Jets(const char *analysis_mode="local", 
			   const char *plugin_mode="",
			   const char *config_file="",Int_t iOffset = 0)
{
// Main analysis train macro. If a configuration file is provided, all parameters
// are taken from there but may be altered by CheckModuleFlags.

   if (strlen(config_file) && !LoadConfig(config_file)) return;

   if(iOffset)kProofOffset = iOffset;
   TString smode(analysis_mode);
   smode.ToUpper();
   if (kSaveTrain)WriteConfig();
   // Check compatibility of selected modules
   CheckModuleFlags(smode);
   //   gROOT->ProcessLine(".trace");

   printf("==================================================================\n");
   printf("===========    RUNNING ANALYSIS TRAIN %s IN %s MODE   ==========\n", kTrainName.Data(),smode.Data());
   printf("==================================================================\n");
   printf("=  Configuring analysis train for:                               =\n");
   if (iAODanalysis) printf("=  AOD analysis                                                  =\n");
   else              printf("=  ESD analysis                                                  =\n");
   if (iPhysicsSelection)   printf("=  Physics selection                                                    =\n");
   if (iESDfilter)   printf("=  ESD filter                                                    =\n");
   if (iJETAN)       printf("=  Jet analysis                                                  =\n");
   printf("==================================================================\n");
   printf(":: use Fill AOD      %d\n", (UInt_t)kFillAOD);
   printf(":: use MC truth      %d\n", (UInt_t)kUseMC);
   printf(":: use KINE filter   %d\n", (UInt_t)kUseKinefilter);
   printf(":: use track refs    %d\n", (UInt_t)kUseTR);
   printf(":: use tags          %d\n", (UInt_t)kUseESDTags);
   printf(":: use AOD tags      %d\n", (UInt_t)kUseAODTags);
   printf(":: use debugging     %d\n", (UInt_t)kUseDebug);
   printf(":: use PAR files     %d\n", (UInt_t)kUsePAR);
   printf(":: use AliEn plugin  %d\n", (UInt_t)kPluginUse);
   printf(":: use PWG1 QA sym       %d\n", iPWG1QASym);
   printf(":: use PWG4 Source Sara  %d\n",iPWG4TmpSourceSara);
   printf(":: use PWG4 Source BB    %d\n",iPWG4TmpSourceFrag);
   printf(":: use PWG4 Jet Chem     %d\n",iPWG4JetChem);
   printf(":: use PWG4 Jet tasks    %d\n",iPWG4JetTasks);
   printf(":: use PWG4 Jet Services %d\n",iPWG4JetServices);     
   printf(":: use PWG4 Jet Spectrum %d\n",iPWG4JetSpectrum);
   printf(":: use PWG4 JCORRAN %d\n",iPWG4JCORRAN);
   printf(":: use PWG4 UE           %d \n",iPWG4UE); 
   printf(":: use PWG4 Pt QA MC     %d\n",iPWG4PtQAMC);
   printf(":: use PWG4 Pt Spectra   %d\n",iPWG4PtSpectra);
   printf(":: use PWG4 Pt QA TPC    %d\n",iPWG4PtQATPC);     
   printf(":: use PWG4 Three Jets   %d\n",iPWG4ThreeJets);
   printf(":: use PWG4 KMeans       %d\n",iPWG4KMeans);
   printf(":: use PWG4 Cluster       %d\n",iPWG4Cluster);
   printf(":: use PWG4 Part Corr    %d\n",iPWG4PartCorr);
   printf(":: use PWG4 Calo QA    %d\n",iPWG4CaloQA);
   printf(":: use PWG4 Jet Corr    %d\n",iPWG4JetCorr);
   printf(":: use PWG4 Tagged       %d\n",iPWG4Tagged);
   printf(":: use PWG4 omega to 3 pions %d\n",iPWG4omega3pi);

   printf(":: use PWG4 Gamma Conv   %d\n",iPWG4GammaConv);
   printf(":: use HighPt FilterMask %d\n",kHighPtFilterMask);    
   
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
      if (iPWG4JetTasks) aodH->AddFriend(Form("deltas/%s",kDeltaAODJetName.Data()));
      if (iPWG4PartCorr) aodH->AddFriend(Form("deltas/%s"kDeltaAODJetName.Data()));
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
      
      mgr->SetOutputEventHandler(aodHandler);
      //
      if (iAODanalysis) {

	//         aodHandler->SetCreateNonStandardAOD();
	//	if (iJETAN)aodHandler->SetOutputFileName(kDeltaAODJetName.Data());
      } 
      AliAnalysisDataContainer * cout_aod = mgr->GetCommonOutputContainer();
      cout_aod->SetSpecialOutput();
   }
   // Debugging if needed

   if (kUseDebug){
     mgr->SetDebugLevel(3);
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
   if(iPhysicsSelection && !iAODanalysis){
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kIsMC,kTRUE,kTRUE); // last flag also adds information on  
   }

   if (iESDfilter && !iAODanalysis) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
      AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kUseKinefilter,kUseMuonfilter);
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
    
    // Jet analysis
   if (iJETAN) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJets.C");
      AliAnalysisTaskJets *taskjets = 0;
      if(iJETAN&1)taskjets = AddTaskJets(kHighPtFilterMask); 
      if(iJETAN&2){
	UInt_t selection = 0;
	if(!kFillAOD) selection = 0xffffff; //&~(1<<13);
	else selection = 0xffffff&~(1<<13);// selection = 1<<0|1<<1|1<<2|1<<5|1<<6|1<<7|1<<8|1<<9|1<<10|1<<11|1<<12;
	AddTaskJetsDelta(kDeltaAODJetName.Data(),kHighPtFilterMask,kUseAODMC,selection); 
      }
      if (!taskjets) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
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
     gROOT->LoadMacro("$ALICE_ROOT/PWG1/PWG1macros/AddTaskQAsym.C");
     AliAnalysisTaskQASym *taskQASym = AddTaskQAsym(-1);
     if (!taskQASym) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQASym cannot run for this train conditions - EXCLUDED");
   }



   if(iPWG4TmpSourceSara){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskEta.C");
     AliAnalysisTaskEta *taskEta = AddTaskEta();
     if (!taskEta) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskEta cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4TmpSourceFrag){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskFragFunc.C");
     UInt_t selection = 1<<2;
     AliAnalysisTaskFragFunc *taskFrag = AddTaskFragFunc(kHighPtFilterMask, kUseAODMC,iPhysicsSelection, selection);
     if (!taskFrag) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskFragFunc cannot run for this train conditions - EXCLUDED");
   }


   if(iPWG4JetChem){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetChem.C");
     AliAnalysisTask *taskChem = AddTaskJetChem();
     if (!taskChem) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetChem cannot run for this train conditions - EXCLUDED");
   }


   if(iPWG4JetServices){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetServices.C");
     AliAnalysisTaskJetServices *taskjetServ = 0;
     taskjetServ = AddTaskJetServices();
     if (!taskjetServ) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetServices cannot run for this train conditions - EXCLUDED");
     if(kGridRunRange[0]>0)taskjetServ->SetRunRange(kGridRunRange[0],kGridRunRange[1]);
     else taskjetServ->SetRunRange(104000,125000);
     if(!kIsMC) taskjetServ->SetRealData(kTRUE);
     taskjetServ->SetUsePhysicsSelection((Bool_t)iPhysicsSelection);
     taskjetServ->SetDebugLevel(0);
   }

   if(iPWG4JetSpectrum){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetSpectrum2.C");
     AliAnalysisTaskJetSpectrum2 *taskjetSpectrum = 0;
     if(iPWG4JetSpectrum&1){
       taskjetSpectrum = AddTaskJetSpectrum2("jets","",kHighPtFilterMask,iPhysicsSelection);      
       if(!iAODanalysis){
	 taskjetSpectrum = AddTaskJetSpectrum2("jets","tracks32",32,iPhysicsSelection);       // tmp hack to give it a different name
	 //	 taskjetSpectrum = AddTaskJetSpectrum2("jets","tracks64",64,iPhysicsSelection);      
       }
       if (!taskjetSpectrum) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskJetSpectrum2 cannot run for this train conditions - EXCLUDED");
       taskjetSpectrum->SetDebugLevel(1);
     }

     if(iPWG4JetSpectrum&2){
       UInt_t selection = 0;
       if(!iAODanalysis) selection = 0xffffff;
       else selection = 1<<0|1<<1|1<<2|1<<3|1<<4|1<<5|1<<7|1<<8|1<<9;
       AddTaskJetSpectrum2Delta(kHighPtFilterMask,kUseAODMC,iPhysicsSelection,selection);
     }
   }
   if(iPWG4JCORRAN){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJCORRANTask.C");
     AliJCORRANTask* corran = AddTaskJCORRAN(kTRUE);
     if(!corran)::Warning("AnalysisTrainPWG4Jets", "AliJCORRANTask cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4UE){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskUE.C");
     AliAnalysisTaskUE *taskUE = 0;
     if(iPWG4UE&1)taskUE = AddTaskUE(); 
     if(iPWG4UE&2){
       taskUE  =AddTaskUE("jetsAOD_CDF04","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       taskUE  =AddTaskUE("jetsAOD_CDF07","CDF", "LJ", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_SISCONE04","CDF", "LJ", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_SISCONE07","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       taskUE  =AddTaskUE("jetsAOD_ICDF","CDF","LJ","TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_FASTKT04","CDF", "LJ", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_FASTKT07","CDF", "LJ", "TRANSV","MSP"); //finder not yet in train
       taskUE  =AddTaskUE("jetsAOD_NONE","CDF", "MP_eta05", "TRANSV","MSP");
       taskUE  =AddTaskUE("jetsAOD_NONE","CDF", "MP_eta09", "TRANSV","MSP");
     }

     if (!taskUE) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskUE cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4ThreeJets){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskThreeJets.C");
     AliAnalysisTaskThreeJets *taskThree = AddTaskThreeJets();
     if(!taskThree)::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskThreets cannot run for this train conditions - EXCLUDED");
   }
   if(iPWG4PtQAMC){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtQAMC.C");
     AliPWG4HighPtQAMC *taskQAMC = AddTaskPWG4HighPtQAMC();
     if (!taskQAMC) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQAMC cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4PtQATPC){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtQATPConly.C");
     AliPWG4HighPtQATPConly *taskQATPC = 0;
     if(iPWG4PtQATPC&1)taskQATPC = AddTaskPWG4HighPtQATPConly(1);
     if(iPWG4PtQATPC&2)taskQATPC = AddTaskPWG4HighPtQATPConly(2);

 if (!taskQATPC) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskQATPC cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4PtSpectra){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPWG4HighPtSpectra.C");
     AliPWG4HighPtSpectra *taskPtSpectra = AddTaskPWG4HighPtSpectra();
     if (!taskPtSpectra) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskPtSpectra cannot run for this train conditions - EXCLUDED");
   }
   if(iPWG4KMeans){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskKMeans.C");
     AliAnalysisTaskKMeans *taskKMeans = AddTaskKMeans();
     if (!taskKMeans) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskKMenans cannot run for this train conditions - EXCLUDED");
   }

   if(iPWG4Cluster){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetCluster.C");
     AliAnalysisTaskJetCluster *taskCl = 0;
     if(iPWG4Cluster&1){
       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelection,"KT");
     }
     if(iPWG4Cluster&2){
       UInt_t selection = 0;
       if(!iAODanalysis) selection = 0xffffff;
       else selection = 1<<0|1<<1|1<<2|1<<3|1<<4|1<<5|1<<7|1<<8|1<<9;
       AddTaskJetClusterDelta(kHighPtFilterMask,kUseAODMC,iPhysicsSelection,"KT",selection);
     }
     if(iPWG4Cluster&4){
       UInt_t selection = 0;
       if(!iAODanalysis) selection = 0xffffff;
       else selection = 1<<0|1<<1|1<<2|1<<3|1<<4|1<<5|1<<7|1<<8|1<<9;
       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelection,"ANTIKT");
       AddTaskJetClusterDelta(kHighPtFilterMask,kUseAODMC,iPhysicsSelection,"ANTIKT",selection);
     }


     if (!taskCl) ::Warning("AnalysisTrainPWG4Jets", "AliAnalysisTaskCluster cannot run for this train conditions - EXCLUDED");

   }
   if(iPWG4PartCorr){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
     AliAnalysisTaskParticleCorrelation *taskpartcorrPHOS = AddTaskPartCorr("AOD", "PHOS",kFALSE,kIsMC);
     if (!taskpartcorrPHOS) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation PHOS cannot run for this train conditions - EXCLUDED");
     AliAnalysisTaskParticleCorrelation *taskpartcorrEMCAL = AddTaskPartCorr("AOD", "EMCAL",kFALSE,kIsMC);
     if (!taskpartcorrEMCAL) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation EMCAL cannot run for this train conditions - EXCLUDED");
     if(kDeltaAODPartCorrName.Length()>0)mgr->RegisterExtraFile(kDeltaAODPartCorrName.Data()); // hmm this is written anyway.... but at least we do not register it...
   } 

   if(iPWG4CaloQA){
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
     AliAnalysisTaskParticleCorrelation *taskcaloQA =  AddTaskCalorimeterQA("ESD",kFALSE,kIsMC);
     if(!taskcaloQA)::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation QA cannot run - EXCLUDED");
   } 

   if(iPWG4JetCorr){
     //     using namespace JetCorrelHD;

     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetCorrel.C");
     AliAnalysisTaskJetCorrel *taskjetcorr = AddTaskJetCorrel();
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
      TString gcArguments = "-run-on-train -run-jet -run-omega-meson -run-neutralmeson";
      TString kGCAnalysisCutSelectionId="9002111000";
      gcArguments.Append(Form("-set-cut-selection  %s ",kGCAnalysisCutSelectionId.Data()));
      if(!kIsMC)gcArguments += " -mc-off";
      AliAnalysisTaskGammaConversion * taskGammaConversion = AddTaskGammaConversion(gcArguments,mgr->GetCommonInputContainer());
      gSystem->ChangeDirectory(cdir);
      taskGammaConversion->SelectCollisionCandidates();
      if (!taskGammaConversion) ::Warning("AnalysisTrainNew", "AliAnalysisTaskGammaConversion cannot run for these train conditions - EXCLUDED");
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
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s",kGridDataSet.Data());
       AliAnalysisAlien *gridhandler = (AliAnalysisAlien*)mgr->GetGridHandler();
       printf("=== AnalysisTrainPWG4Jets:: Registering jdl in the work directory alien://%s/%s, should be done by the manager! ===\n",
	      alien_workdir.Data(),gridhandler->GetGridOutputDir());

       TString dest;
       dest = Form("%s/%s/%s.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 //	 Printf("%s exist on grid removing...",dest.Data());
	 //	 gGrid->Rm(dest.Data());
       }
       TFile::Cp(Form("file:%s.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));
       
       dest = Form("%s/rootfiles/STEER/LQ1dRefv1.root",gGrid->GetHomeDirectory());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 Printf("%s exist on grid removing...",dest.Data());
	 gGrid->Rm(dest.Data());
       }
       TFile::Cp(Form("file:%s/STEER/LQ1dRef_v1.root",
		      gSystem->ExpandPathName("$ALICE_ROOT")), 
		 Form("alien://%s",dest.Data()));
     }
     AliLog::SetGlobalLogLevel(AliLog::kError);
     if((kUseSysInfo>0&&smode=="LOCAL")||!strcmp(plugin_mode, "test")){
       TFile * fM = TFile::Open("manager_local.root","RECREATE");
       mgr->Write();
       fM->Close();
     }

     StartAnalysis(smode, chain);
       
     if((kUseSysInfo>0&&smode=="LOCAL")||!strcmp(plugin_mode, "test")){
       for(int i = 0;i < mgr->GetTopTasks()->GetEntries();i++){
	 mgr->ProfileTask(i);
       }
     }
     if (!strcmp(plugin_mode, "offline")&&smode=="GRID"){
       // Offline mode path files
       //	PatchJDL();
       PatchAnalysisMacro();
     }

     if (kSaveTrain && smode=="GRID") {
       AliAnalysisAlien *gridhandler = (AliAnalysisAlien*)mgr->GetGridHandler();
       TString alien_workdir = gGrid->GetHomeDirectory();
       if (iAODanalysis) alien_workdir += "analysisAOD";
       else              alien_workdir += "analysisESD";
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s",kGridDataSet.Data());
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
            mgr->StartAnalysis("grid");
         } else {
            if (!chain) {
               ::Error("AnalysisTrainPWG4Jets.C::StartAnalysis", "Cannot create the chain");
               return;
            }   
            mgr->StartAnalysis(mode, chain);
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
      iPWG$JCORRAN = 0;
      if( iPWG4PtQAMC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtQAMC disabled in analysis on AOD's");
      iPWG4PtQAMC        = 0;
      if( iPWG4PtQATPC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtTPC disabled in analysis on AOD's");
      iPWG4PtQATPC        = 0;
      if( iPWG4PtSpectra)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtQAMC disabled in analysis on AOD's");
      iPWG4PtSpectra     = 0;
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
   } else {   
   // ESD analysis
     if (!kUseMC){
       kUseTR = kFALSE;
       
       if(kUseKinefilter)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 Kine Filter disabled in analysis without MC");
       kUseKinefilter = kFALSE;
       if( iPWG4PtQAMC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "PWG4 PtQAMC disabled in analysis without MC");
       iPWG4PtQAMC        = 0;

     }
     if (!kUseTR) {
       if(iPWG4PtQAMC)::Info("AnalysisTrainPWG4Jets.C::CheckModuleFlags", "iPWG4QATPCMC disabled if not reading track references");
       iPWG4PtQAMC        = 0;
     }   
     if (iJETAN){
       iESDfilter=1;
     }
      if (!iESDfilter){
	kUseKinefilter = kFALSE;
	kUseMuonfilter = kFALSE;
      }
      if(!iJETAN){
	iPWG4JetSpectrum = iPWG4UE = iPWG4ThreeJets = iDIJETAN = 0;
      }
   }
   iPWG4JetTasks = iPWG4JetServices||iPWG4JetSpectrum||iPWG4UE||iPWG4PtQAMC||iPWG4PtSpectra||iPWG4PtQATPC||iPWG4ThreeJets||iPWG4JetChem;
   iPWG4PartCorrLibs = iPWG4PartCorr||iPWG4Tagged||iPWG4CaloQA;
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
// Load common analysis libraries.
   Int_t imode = -1;
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   if (!gSystem->Getenv("ALICE_ROOT")) {
      ::Error("AnalysisTrainPWG4Jets.C::LoadCommonLibraries", "Analysis train requires that analysis libraries are compiled with a local AliRoot"); 
      return kFALSE;
   }   
   Bool_t success = kTRUE;
   // ROOT libraries
   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   
   // Load framework classes. Par option ignored here.
   switch (imode) {
      case 0:
      case 2:
         if (kUseCPAR) {
            success &= LoadLibrary("STEERBase", mode, kTRUE);
            success &= LoadLibrary("ESD", mode, kTRUE);
            success &= LoadLibrary("AOD", mode, kTRUE);
            success &= LoadLibrary("ANALYSIS", mode, kTRUE);
            success &= LoadLibrary("ANALYSISalice", mode, kTRUE);
            success &= LoadLibrary("CORRFW", mode, kTRUE);
         } else {   
            success &= LoadLibrary("libSTEERBase.so", mode);
            success &= LoadLibrary("libESD.so", mode);
            success &= LoadLibrary("libAOD.so", mode);
            success &= LoadLibrary("libANALYSIS.so", mode);
            success &= LoadLibrary("libANALYSISalice.so", mode);
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
   if (iESDfilter) {
      if (!LoadLibrary("PWG3base", mode, kTRUE) ||
          !LoadLibrary("PWG3muon", mode, kTRUE)) return kFALSE;
   }   
   // JETAN
   if (iJETANLib) {
     // this part needs some rework in case we do not need the fastjed finders for processing
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
   if(iPWG4JCORRAN){
     // PWG4 particle correlations
     if(!LoadLibrary("PWG4JCORRAN",mode,kTRUE))return kFALSE;
   }

   if(iPWG1QASym){
     if (!LoadSource(Form("%s/PWG1/AliAnalysisTaskQASym.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }
   if(iPWG4TmpSourceSara){
     if(!kUsePAR)gSystem->AddIncludePath("-I$ALICE_ROOT/include/JetTasks"); // ugly hack!!
     if(!LoadSource(Form("%s/PWG4/JetTasks/AliAnalysisTaskEta.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }
   if(iPWG4TmpSourceFrag){
     if(!kUsePAR)gSystem->AddIncludePath("-I$ALICE_ROOT/include/JetTasks"); // ugly hack!!
     if(!LoadSource(Form("%s/PWG4/JetTasks/AliAnalysisTaskFragFunc.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }


   /*
   if(iPWG4JetChem){
     if(!kUsePAR)gSystem->AddIncludePath("-I$ALICE_ROOT/include/JetTasks"); // ugly hack!!
     if(!LoadSource(Form("%s/PWG4/JetTasks/AliAnalysisTaskJetChem.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }
   */

   if (iPWG4PartCorrLibs) {   
      if (!LoadLibrary("EMCALUtils", mode, kTRUE) ||
          !LoadLibrary("PHOSUtils", mode, kTRUE) ||
          !LoadLibrary("PWG4PartCorrBase", mode, kTRUE) ||
          !LoadLibrary("PWG4PartCorrDep", mode, kTRUE)) return kFALSE;
   }
   if (iPWG4JetCorr) { 
     if (!LoadLibrary("PWG4JetCorrel", mode, kTRUE)) return kFALSE;
   }  
   if (iPWG4omega3pi) {
     if (!LoadLibrary("PWG4omega3pi", mode, kTRUE)) return kFALSE;
   }
   if (iPWG4GammaConv) {
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
	       while (in.good()) {
		 in >> line;
		 if (line.Length() == 0) continue;
		 // cout << " line = " << line << endl;
		 chain->Add(line.Data());
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
	if (!gSystem->ChangeDirectory(ocwd.Data())) return -1;
   return 0;
}

//______________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(plugin_mode);
   if (kPluginUseProductionMode) plugin->SetProductionMode();
   plugin->SetJobTag(kJobTag);
   plugin->SetNtestFiles(2);
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
   if (iAODanalysis) plugin->SetDataPattern(" *AliAOD.root");
   else              plugin->SetDataPattern(Form(" %s/*/*ESD.tag.root",kGridPassPattern.Data()));
// ...then add run numbers to be considered
//   plugin->SetRunRange(kGridRunRange[0], kGridRunRange[1]);
   for (Int_t i=kGridRunRange[0]; i<=kGridRunRange[1]; i++) {
     Printf("AnalysisTrainPWG4Jets Adding run number %s", Form(kGridRunPattern.Data(),i));
     plugin->AddRunNumber(Form(kGridRunPattern.Data(),i));
   }   

   if(kGridLocalRunList.Length()>0){
     ifstream in1;
     in1.open(kGridLocalRunList.Data());
     int iRun;

     /*
     char c;
     char cLine[250];
     while(!in1.eof()){
       c = in1.get();
       if ( (c >= '0') && (c <= '9') )
	 {
	   in1.putback (c);
	   in1>>iRun;
	   Printf("AnalysisTrainPWG4Jets Adding run number from File %s", Form(kGridRunPattern.Data(),iRun));
	   plugin->AddRunNumber(Form(kGridRunPattern.Data(),iRun));
       }
     else
       {
	 in1.putback (c);
	 in1.getline(cLine,250);

       }
     }
     */

     // just use run numbers, negatives will be excluded
     while(in1>>iRun){
       if(iRun>0){
	   Printf("AnalysisTrainPWG4Jets Adding run number from File %s", Form(kGridRunPattern.Data(),iRun));
	   plugin->AddRunNumber(Form(kGridRunPattern.Data(),iRun));
       }
       else{
	 Printf("AnalysisTrainPWG4Jets Skipping run number from File %d", iRun);
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
   if(kGridDataSet.Length()>0)alien_workdir += Form("/%s",kGridDataSet.Data());
   plugin->SetGridWorkingDir(alien_workdir.Data());

   // Declare alien output directory. Relative to working directory.
   if (!kGridOutdir.Length()) kGridOutdir = Form("output_%s",kTrainName.Data());
   plugin->SetGridOutputDir(kGridOutdir);

   // Add external packages
   if (iJETAN||iDIJETAN) {
      plugin->AddExternalPackage("boost::v1_38_0");
      plugin->AddExternalPackage("cgal::v3.3.1");
      plugin->AddExternalPackage("fastjet::v2.4.0");
   }   


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
   plugin->SetDefaultOutputs();
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

   if(kUseSysInfo>0){
     if (listhists.Length()) listhists += " ";
     listhists += "syswatch.root";
   }

   if(kIsMC){
     if (listaods.Length()) listaods += " ";
     listaods += "pyxsec_hists.root";
   }


   if (mgr->GetExtraFiles().Length()) {
     if (listaods.Length()) listaods += " ";
     listaods += mgr->GetExtraFiles();
   }

   // if we do not fill the aod we do not need to store it
   kGridMergeExclude = listaods;
   
   if(!kFillAOD){
     listaods="";
     plugin->SetDefaultOutputs(kFALSE);
     plugin->SetOutputFiles(listhists.Data());
   }


   listaods.ReplaceAll(" ", ",");
   listhists.ReplaceAll(" ", ",");
   if (listhists.Length()) listhists = Form("hist_archive.zip:%s@%s", listhists.Data(), kGridOutputStorages.Data());;
   if (listaods.Length())  listaods  = Form("aod_archive.zip:%s@%s", listaods.Data(), kGridOutputStorages.Data());;

   if (!listhists.Length() && !listaods.Length()) {
      ::Fatal("AnalysisTrainPWG4Jets", "No task output !");
   }



   TString outputArchive = "log_archive.zip:stdout,stderr@disk=2";
   if(kUseSysInfo>0)outputArchive = "log_archive.zip:stdout,stderr,syswatch.log@disk=2";
   if (listaods.Length()) {
      outputArchive += " ";
      outputArchive += listaods;
   }   
   if (listhists.Length()) {
      outputArchive += " ";
      outputArchive += listhists;
   }   
   plugin->SetOutputArchive(outputArchive);
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", kTrainName.Data()));
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(kGridFilesPerJob);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
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
  Int_t index;
  index = st.Index("gSystem->Load(\"libPhysics\");");
  index += strlen("gSystem->Load(\"libPhysics\");");
  if(iJETAN&&kUsePAR){
    TObjArray *arr;
    TObjString *objstr;
    arr = anaLibs.Tokenize(" ");
    TIter next(arr);
    TString add = "";
    add += "\n\n // added by CKB \n";
    while ((objstr=(TObjString*)next())){
      if(objstr->GetString().Contains("PWG3"))continue;
      if(objstr->GetString().EndsWith(".so"))add += Form("gSystem->Load(\"%s\");\n",objstr->GetString().Data());
    }
    delete arr; 
    //    add += Form("AliLog::SetGlobalLogLevel(%d);\n",AliLog::GetGlobalLogLevel());
    add += "gSystem->AddIncludePath(\"./\") \n";
    if(gGrid&&kPluginAliRootVersion.Length()==0){
      add += "// Dirty hack for TRD reference data";
      add += "gSystem->Setenv(\"ALICE_ROOT\",\"";
      add += Form("alien://%s/rootfiles/",gGrid->GetHomeDirectory());
      add += "\"); \n";
    }
    add += "// BKC \n\n";
    st.Insert(index,add.Data());
  }

  if(kUseDebug){
    st.Insert(index,"\n gROOT->ProcessLine(\".trace\"); // CKB \n");
  }

  index = st.Index("gSystem->AddIncludePath");
  st.Insert(index,"// CKB ");

  ofstream out;
  out.open(Form("%s.C", kTrainName.Data()));
  if (out.bad()) {
    return kFALSE;
  }
  out << st << endl;
  Printf("<<< Patching AnalysisMacro");
  return kTRUE;

}
