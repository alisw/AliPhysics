//===================== ANALYSIS TRAIN =========================================
// To use: copy this macro to your work directory, modify the global part to match
// your needs, then run root.
//    root[0] .L AnalysisTrainHMPID.C
// Grid full mode as below (other modes: test, offline, submit, terminate)
//    root[1] AnalysisTrainHMPID("grid", "full")
// CAF mode (requires root v5-23-02 + aliroot v4-16-Rev08)
//    root[2] AnalysisTrainHMPID("proof")
// Local mode requires AliESds.root or AliAOD.root in ./data directory
//    root[3] AnalysisTrainHMPID("local")
// In proof and grid modes, a token is needed and sourcing the produced environment file.
//

// =============================================================================
// ### General Steering variables
// =============================================================================
//== general setup variables
TString     kTrainName         = "hmpidAnalysis"; // (no blancs or special characters)
TString     kJobTag            = "HMPID Tasks analysis train configured"; //
Bool_t      kUsePAR            = kFALSE;  // use par files for extra libs
Bool_t      kUseCPAR           = kFALSE;  // use par files for common libs
Bool_t      kFillAOD           = kFALSE;  // switch of AOD filling for on the fly analysis

Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iESDfilter         = 0;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iPhysicsSelection  = 1;      // Physics selection task
Bool_t      kUseKinefilter     = kFALSE; // use Kinematics filter
Bool_t      kUseMuonfilter     = kFALSE; // use muon filter
TString     kCommonOutputFileName = "HmpidOutput.root";


//== general process variables

// ### Other flags to steer the analysis
//==============================================================================
Bool_t      kSkipTerminate      = kFALSE; // Do not call Teminate
Bool_t      kDebugLevel         = kTRUE; // activate debugging
Int_t       kUseSysInfo         = 0; // activate debugging
Bool_t      kUseMC              = kTRUE;  // use MC info
Bool_t      kIsMC               = kTRUE;  // is MC info, if false it overwrites Use(AOD)MC
Bool_t      kUseESDTags         = kTRUE; // use ESD tags for selection
Bool_t      kUseTR              = kFALSE;  // use track references

// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iHMPID             = 1;      // Basic HMPID analysis task
Int_t       iJETAN             = 0;      // Jet analysis (PWG4) // 1 write standard 2 write non-standard jets
Int_t       iJETANLib          = 0;
Int_t       kHighPtFilterMask  = 16;     // change depending on the used AOD Filter


//==============================================================================
// ### PROOF Steering varibales
//==============================================================================
//== proof setup variables
TString     kProofCluster      = "alice-caf.cern.ch";
Bool_t      kProofUseAFPAR     = kTRUE;  // use AF special par file
TString     kProofAFversion    = "VO_ALICE@AliRoot::v4-20-08-AN";
//== proof input and output variables
TString     kProofDataSet      = "/alice/sim/LHC10d2_117220";
//== proof process variables
Bool_t      kProofClearPackages = kFALSE;
Int_t       kProofEvents = 10000;
Int_t       kProofOffset = 0;

//==============================================================================
// ### Grid plugin Steering varibiables
//==============================================================================
//== grid plugin setup variables
Bool_t      kPluginUse         = kTRUE;   // do not change
Bool_t      kPluginUseProductionMode  = kFALSE;   // use the plugin in production mode
TString     kPluginRootVersion       = "v5-27-06b";  // *CHANGE ME IF MORE RECENT IN GRID*
TString     kPluginAliRootVersion    = "v4-21-05-AN";  // *CHANGE ME IF MORE RECENT IN GRID*                                          
Bool_t      kPluginMergeViaJDL       = kTRUE;  // merge via JDL
Bool_t      kPluginFastReadOption    = kFALSE;  // use xrootd flags to reduce timeouts
Bool_t      kPluginOverwriteMode     = kTRUE;  // overwrite existing collections
Int_t       kPluginOutputToRunNumber = 1;     // write the output to subdirs named after run number
TString kPluginExecutableCommand = "root -b -q";

// == grid plugin input and output variables
TString     kGridDatadir      = "/alice/sim/LHC10d4a";
TString     kGridLocalRunList = "";
TString     kGridWorkDir      = "HmpidAnalysis/LHC10d4a";   // Alien working directory
TString     kGridOutdir       = ""; // AliEn output directory. If blank will become output_<kTrainName>
TString     kGridDataSet      = ""; // sub working directory not to confuse different run xmls 
Int_t       kGridRunRange[2]       = {120820, 120820}; // Set the run range
TString     kGridRunPattern        = "%03d"; // important for leading zeroes!!
TString     kGridPassPattern       = "";
TString     kGridExtraFiles        = ""; // files that will be added to the input list in the JDL...
Int_t       kGridMaxMergeFiles      = 12; // Number of files merged in a chunk grid run range
TString     kGridMergeExclude       = "AliAOD.root"; // Files that should not be merged
TString     kGridOutputStorages      = "disk=2"; // Make replicas on the storages
// == grid process variables
Int_t       kGridRunsPerMaster     = 1; // Number of runs per master job
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
void AnalysisTrainHMPID(const char *analysis_mode="local", const char *plugin_mode="",
                        const char *config_file="",Int_t iOffset = 0)
{
// Main analysis train macro. If a configuration file is provided, all parameters
// are taken from there but may be altered by CheckModuleFlags.

   if (strlen(config_file) && !LoadConfig(config_file)) return;

   if(iOffset)kProofOffset = iOffset;
   TString smode(analysis_mode);
   smode.ToUpper();
   // Check compatibility of selected modules
   CheckModuleFlags(smode);

   printf("==================================================================\n");
   printf("===========    RUNNING ANALYSIS TRAIN %s IN %s MODE   ==========\n", kTrainName.Data(),smode.Data());
   printf("==================================================================\n");
   printf("=  Configuring analysis train for:                               =\n");
   printf("=  ESD analysis                                                  =\n");
   if (iPhysicsSelection)   printf("=  Physics selection                                                    =\n");
   if (iESDfilter)   printf("=  ESD filter                                                    =\n");
   if (iJETAN)       printf("=  Jet analysis                                                  =\n");
   printf("==================================================================\n");
   printf(":: use MC truth      %d\n", (UInt_t)kUseMC);
   printf(":: use KINE filter   %d\n", (UInt_t)kUseKinefilter);
   printf(":: use track refs    %d\n", (UInt_t)kUseTR);
   printf(":: use tags          %d\n", (UInt_t)kUseESDTags);
   printf(":: use debugging     %d\n", (UInt_t)kDebugLevel);
   printf(":: use PAR files     %d\n", (UInt_t)kUsePAR);
   printf(":: use AliEn plugin  %d\n", (UInt_t)kPluginUse);
   
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
   AliAnalysisManager *mgr  = new AliAnalysisManager("HMPIDTrain", "HMPID train");
   if (kCommonOutputFileName.Length()>0)mgr->SetCommonFileName(kCommonOutputFileName.Data());
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
   // ESD input handler
   AliESDInputHandler *esdHandler = new AliESDInputHandler();
   if (kUseESDTags) esdHandler->SetReadTags();
   esdHandler->SetReadFriends(kFALSE);
   mgr->SetInputEventHandler(esdHandler);       

   // Monte Carlo handler
   if (kUseMC) {
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
      AliAnalysisDataContainer *cout_aod = mgr->GetCommonOutputContainer();
      cout_aod->SetSpecialOutput();
   }

   // Debugging if needed
   if (kDebugLevel){
      mgr->SetDebugLevel(3);
   }
   if(kUseSysInfo>0){
      mgr->RegisterExtraFile("syswatch.root");
      if(kGridMergeExclude.Length())kGridMergeExclude += " ";
      kGridMergeExclude += "syswatch.root";
   } else {
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
   if(iPhysicsSelection){
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kIsMC,kTRUE,kTRUE); // last flag also adds information on  
   }

   if (iESDfilter) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
      AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kUseKinefilter,kUseMuonfilter);
      if(kIsMC){
         mgr->RegisterExtraFile("pyxsec_hists.root");
         if(kGridMergeExclude.Length())kGridMergeExclude += " ";
         kGridMergeExclude += "pyxsec_hists.root";
      }
   }   

    // Jet analysis
   if (iJETAN) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJets.C");
      AliAnalysisTaskJets *taskjets = 0;
      if (iJETAN&1) taskjets = AddTaskJets(kHighPtFilterMask); 
      if (!taskjets) ::Warning("AnalysisTrainHMPID", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
   }

   if(iHMPID){
     gROOT->LoadMacro("$ALICE_ROOT/HMPID/AddTaskHMPID.C");
     AliHMPIDAnalysisTask *taskHmpid = AddTaskHMPID(kUseMC);
     if (!taskHmpid) ::Warning("AnalysisTrainHMPID", "AliAnalysisTaskHMPID cannot run for this train conditions - EXCLUDED");
   }

   if (kPluginUse) {
      AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
      AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
   }

   if (mgr->InitAnalysis()) {
     mgr->PrintStatus();
     if (!strcmp(plugin_mode,"submit") && smode=="GRID"){
       TString alien_workdir = gGrid->GetHomeDirectory();
       alien_workdir += kGridWorkDir.Data();
       if(kGridDataSet.Length()>0)alien_workdir += Form("/%s",kGridDataSet.Data());
       AliAnalysisAlien *gridhandler = (AliAnalysisAlien*)mgr->GetGridHandler();
       printf("=== AnalysisTrainHMPID:: Registering jdl in the work directory alien://%s/%s, should be done by the manager! ===\n",
	      alien_workdir.Data(),gridhandler->GetGridOutputDir());

       TString dest;
       dest = Form("%s/%s/%s.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
     //  Printf("%s exist on grid removing...",dest.Data());
     //  gGrid->Rm(dest.Data());
       }
       TFile::Cp(Form("file:%s.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));


       TString dest;
       dest = Form("%s/%s/%s_merge.jdl",alien_workdir.Data(),gridhandler->GetGridOutputDir(),kTrainName.Data());
       if(AliAnalysisAlien::FileExists(dest.Data())){
	 //	 Printf("%s exist on grid removing...",dest.Data());
	 //	 gGrid->Rm(dest.Data());
       }
       TFile::Cp(Form("file:%s_merge.jdl",kTrainName.Data()),Form("alien://%s",dest.Data()));
     }

     AliLog::SetGlobalLogLevel(AliLog::kError);
     if((kUseSysInfo>0 && smode=="LOCAL") || !strcmp(plugin_mode, "test")){
       TFile *fM = TFile::Open("manager_local.root","RECREATE");
       mgr->Write();
       fM->Close();
     }

     StartAnalysis(smode, chain);
       
     if((kUseSysInfo>0 && smode=="LOCAL") || !strcmp(plugin_mode, "test")){
       for(int i = 0;i < mgr->GetTopTasks()->GetEntries();i++){
         mgr->ProfileTask(i);
       }
     }
     if (!strcmp(plugin_mode, "offline") && smode=="GRID"){
       // Offline mode path files
       //	PatchJDL();
       PatchAnalysisMacro();
     }
   }
}


//______________________________________________________________________________
void StartAnalysis(const char *mode, TChain *chain) {

   Int_t imode = -1;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!strcmp(mode, "LOCAL")) imode = 0;
   if (!strcmp(mode, "PROOF")) imode = 1;
   if (!strcmp(mode, "GRID"))  imode = 2;
   switch (imode) {
      case 0:
         if (!chain) {
            ::Error("AnalysisTrainHMPID.C::StartAnalysis", "Cannot create the chain");
            return;
         }
         mgr->StartAnalysis(mode, chain,kNumberOfEvents);
         return;
      case 1:
         if (!kProofDataSet.Length()) {
            ::Error("AnalysisTrainHMPID.C::StartAnalysis", "kProofDataSet is empty");
            return;
         }
         mgr->StartAnalysis(mode, kProofDataSet, kProofEvents,kProofOffset);
         return;
      case 2:
         if (kPluginUse) {
            if (!mgr->GetGridHandler()) {
               ::Error("AnalysisTrainHMPID.C::StartAnalysis", "Grid plugin not initialized");
               return;
            }
            mgr->StartAnalysis("grid");
         } else {
            if (!chain) {
               ::Error("AnalysisTrainHMPID.C::StartAnalysis", "Cannot create the chain");
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
         ::Info("AnalysisTrainHMPID.C::CheckModuleFlags", "PAR files enabled due to PROOF analysis");
         kUsePAR = kTRUE;
      }   
   }  
   if (imode != 2) {
      ::Info("AnalysisTrainHMPID.C::CheckModuleFlags", "AliEn plugin disabled since not in GRID mode");
      kPluginUse = kFALSE; 
   }

   if(!kIsMC){
     // switch off anthin related to MC
     kUseMC = 0;
     kUseTR = kFALSE;
   }

 // ESD analysis
   if (!kUseMC){
     kUseTR = kFALSE;
     if(kUseKinefilter)::Info("AnalysisTrainHMPID.C::CheckModuleFlags", "Kine Filter disabled in analysis without MC");
     kUseKinefilter = kFALSE;
   }
   if (iJETAN){
     iESDfilter=1;
   }
   if (!iESDfilter){
     kUseKinefilter = kFALSE;
     kUseMuonfilter = kFALSE;
   }

   iJETANLib = iJETAN && 1;
   if (iESDfilter) {iAODhandler=1;}
   if (kUseKinefilter && !kUseMC) kUseKinefilter = kFALSE;
   
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
         ::Error(Form("AnalysisTrainHMPID.C::Connect <%s>", mode), "Make sure you:\n \
                        1. Have called: alien-token-init <username>\n \
                        2. Have called: >source /tmp/gclient_env_$UID");
         return kFALSE;
       }
       ::Info("AnalysisTrainHMPID.C::Connect", "Connecting user <%s> to PROOF cluster <%s>", 
                username.Data(), kProofCluster.Data());
       gEnv->SetValue("XSec.GSI.DelegProxy", "2");
       TProof::Open(Form("%s@%s", username.Data(), kProofCluster.Data()));       
       if (!gProof) {
         if (strcmp(gSystem->Getenv("XrdSecGSISRVNAMES"), "lxfsrd0506.cern.ch"))
           ::Error(Form("AnalysisTrainHMPID.C::Connect <%s>", mode), "Environment XrdSecGSISRVNAMES different from lxfsrd0506.cern.ch");
           return kFALSE;
         }
       if(kProofClearPackages)gProof->ClearPackages();
       break;
     case 2:      
       if  (!username.Length()) {
         ::Error(Form("AnalysisTrainHMPID.C::Connect <%s>", mode), "Make sure you:\n \
                        1. Have called: alien-token-init <username>\n \
                        2. Have called: >source /tmp/gclient_env_$UID");
         return kFALSE;
       }
       if (kPluginUse && !gSystem->Getenv("alien_CLOSE_SE")) {
         ::Error(Form("AnalysisTrainHMPID.C::Connect <%s>", mode), 
                        "When using the AliEn plugin it is preferable to define the \
                        variable alien_CLOSE_SE in your environment.");
         return kFALSE;
       }
       ::Info("AnalysisTrainHMPID.C::Connect", "Connecting user <%s> to AliEn ...", 
                username.Data());
       TGrid::Connect("alien://");
       if (!gGrid || !gGrid->IsConnected()) return kFALSE;
       break;
     default:
       ::Error("AnalysisTrainHMPID.C::Connect", "Unknown run mode: %s", mode);
       return kFALSE;
   }
   ::Info("AnalysisTrainHMPID.C::Connect","Connected in %s mode", mode);
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
      ::Error("AnalysisTrainHMPID.C::LoadCommonLibraries", "Analysis train requires that analysis libraries are compiled with a local AliRoot"); 
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
         if (!kProofUseAFPAR) {
            success &= LoadLibrary("STEERBase", mode);
            success &= LoadLibrary("ESD", mode);
            success &= LoadLibrary("AOD", mode);
            success &= LoadLibrary("ANALYSIS", mode);
            success &= LoadLibrary("ANALYSISalice", mode);
            success &= LoadLibrary("CORRFW", mode);
         } else { 
            success &= !gProof->EnablePackage(kProofAFversion);
            success &= LoadLibrary("CORRFW", mode);
         }
         break;         
      default:
         ::Error("AnalysisTrainHMPID.C::LoadCommonLibraries", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   if (success) {
      ::Info("AnalysisTrainHMPID.C::LoadCommodLibraries", "Load common libraries:    SUCCESS");
      ::Info("AnalysisTrainHMPID.C::LoadCommodLibraries", "Include path for Aclic compilation:\n%s",
              gSystem->GetIncludePath());
   } else {           
      ::Info("AnalysisTrainHMPID.C::LoadCommodLibraries", "Load common libraries:    FAILED");
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

   if(iHMPID){
     if (!LoadSource(Form("%s/HMPID/AliHMPIDAnalysisTask.cxx",gSystem->ExpandPathName("$ALICE_ROOT")), mode, kTRUE))return kFALSE;
   }

   if (iJETANLib) {
     if (!LoadLibrary("JETAN", mode, kTRUE)) return kFALSE;
   }

   ::Info("AnalysisTrainHMPID.C::LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
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
      ::Error("AnalysisTrainHMPID.C::LoadLibrary", "Empty module name");
      return kFALSE;
   }   
   // If a library is specified, just load it
   if (smodule.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      result = gSystem->Load(mod);
      if (result < 0) {
         ::Error("AnalysisTrainHMPID.C::LoadLibrary", "Could not load library %s", module);
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
	       ::Info("AnalysisTrainHMPID.C::LoadLibrary", "Removing directory %s",module);
	       gSystem->Exec(Form("rm -rf %s",module));
         }
         result = gProof->UploadPackage(module);
         if (result<0) {
            result = gProof->UploadPackage(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", module)));
            if (result<0) {
               ::Error("AnalysisTrainHMPID.C::LoadLibrary", "Could not find module %s.par in current directory nor in $ALICE_ROOT", module);
               return kFALSE;
            }
         }   
         result = gProof->EnablePackage(module);
         break;
      default:
         return kFALSE;
   }         
   if (result < 0) {
      ::Error("AnalysisTrainHMPID.C::LoadLibrary", "Could not load module %s", module);
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
      ::Error("AnalysisTrainHMPID.C::LoadSource", "Empty task name");
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
       result = gProof->Load(Form("%s.cxx++g",basename.Data()));
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
      ::Error("AnalysisTrainHMPID.C::LoadSources", "Could not load source %s", source);
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
         break;
      case 1:
         break;
      case 2:
         if (kPluginUse) {
            AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
            AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
         } else {
            TString treeName = "esdTree";
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
      ::Error("AnalysisTrainHMPID.C::CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
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
         ::Info("AnalysisTrainHMPID.C::SetupPar", "Getting %s.par from $ALICE_ROOT", pararchivename);
         TFile::Cp(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", pararchivename)), 
                   Form("%s.par",pararchivename));
      } else {
         ::Error("AnalysisTrainHMPID.C::SetupPar", "Cannot find %s.par", pararchivename);
         return -1;
      }   
   }
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
   plugin->SetNtestFiles(1);
//   plugin->SetPreferedSE("ALICE::NIHAM::File");
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion(kPluginRootVersion);
   plugin->SetAliROOTVersion(kPluginAliRootVersion);

// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir(kGridDatadir.Data());
// Set data search pattern
   plugin->SetDataPattern(Form(" %s/*/*ESDs.root",kGridPassPattern.Data()));
// ...then add run numbers to be considered
   plugin->SetRunRange(kGridRunRange[0], kGridRunRange[1]);

   if(kGridLocalRunList.Length()>0){
     ifstream in1;
     in1.open(kGridLocalRunList.Data());
     int iRun;
     // just use run numbers, negatives will be excluded
     while(in1>>iRun){
       if(iRun>0){
       Printf("AnalysisTrainHMPID Adding run number from File %s", Form(kGridRunPattern.Data(),iRun));
       plugin->AddRunNumber(Form(kGridRunPattern.Data(),iRun));
       } else{
         Printf("AnalysisTrainHMPID Skipping run number from File %d", iRun);
       }
     }
   }

// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("Hijing.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   TString alien_workdir = "";

   alien_workdir += kGridWorkDir.Data();
   if(kGridDataSet.Length()>0)alien_workdir += Form("/%s",kGridDataSet.Data());
   plugin->SetGridWorkingDir(alien_workdir.Data());

   // Declare alien output directory. Relative to working directory.
   if (!kGridOutdir.Length()) kGridOutdir = Form("output_%s",kTrainName.Data());
   plugin->SetGridOutputDir(kGridOutdir);

   // set extra libs before par file compilation
   anaLibs += kGridExtraFiles;
   anaLibs  = anaLibs.Strip();   
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
   plugin->SetUseSubmitPolicy(kFALSE);
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
//   plugin->SetDefaultOutputs(kFALSE);

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

   TString outputArchive;
   outputArchive = Form("log_archive.zip:std*r@%s",kGridOutputStorages.Data());
   listaods.ReplaceAll(" ", ",");
   listhists.ReplaceAll(" ", ",");
   if (listhists.Length()) listhists = Form("hist_archive.zip:%s@%s", listhists.Data(), kGridOutputStorages.Data());
   if (listaods.Length())  listaods  = Form("aod_archive.zip:%s@%s", listaods.Data(), kGridOutputStorages.Data());

   if (!listhists.Length() && !listaods.Length()) {
      ::Fatal("AnalysisTrainHMPID", "No task output !");
   }

   if (listaods.Length()) {
      outputArchive += " ";
      outputArchive += listaods;
   }
   if (listhists.Length()) {
      outputArchive += " ";
      outputArchive += listhists;
   }
//   plugin->SetOutputArchive(outputArchive);

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
   plugin->SetCheckCopy(kFALSE);
   return plugin;
}

//______________________________________________________________________________
Bool_t LoadConfig(const char *filename)
{
// Read train configuration from file
   if (gSystem->AccessPathName(filename)) {
      ::Error("AnalysisTrainHMPID.C::LoadConfig", "Config file name not found");
      return kFALSE;
   }   
   gROOT->ProcessLine(Form(".x %s", filename));
   ::Info("AnalysisTrainHMPID.C::LoadConfig", "Train configuration loaded from file %s", filename);
   return kTRUE;
}

//______________________________________________________________________________
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

//______________________________________________________________________________
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
  add += "\n gSystem->AddIncludePath(\"./\"); \n";
  if(gGrid && kPluginAliRootVersion.Length()==0){
    add += "\n // Dirty hack for TRD reference data \n";
    add += "\n gSystem->Setenv(\"ALICE_ROOT\",\"";
    add += Form("alien://%s/rootfiles/",gGrid->GetHomeDirectory());
    add += "\"); \n";
  }
  add += "// BKC \n\n";
  st.Insert(index,add.Data());

  if(kUseCPAR && kPluginAliRootVersion.Length()==0){
    index = st.Index("gSystem->AddIncludePath(\"-I$"); // uncommen $ALICE_ROOT include for par files
    if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
    st.Insert(index,"// CKB comment out whehn no aliroot is provided \n //");
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
    add2 += "\n // Dirty hack for TRD reference data \n";
    add2 += "\n gSystem->Setenv(\"ALICE_ROOT\",\"";
    add2 += Form("alien://%s/rootfiles/",gGrid->GetHomeDirectory());
    add2 += "\"); \n";
  }
  add2 += "// BKC \n\n";
  if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
  st2.Insert(index,add.Data());

  if(kUseCPAR&&kPluginAliRootVersion.Length()==0){
    index = st2.Index("gSystem->AddIncludePath(\"-I$"); // uncommen $ALICE_ROOT include for par files
    if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
    st2.Insert(index,"// CKB comment out whehn no aliroot is provided \n //");
  }

  // do not exclude the extra files from merign, this is done explicitly in this train script
  index = st2.Index("mergeExcludes +="); // uncommen $ALICE_ROOT include for par files
  if(index<0)Printf("%s:%d index out of bounds",(char*)__FILE__,__LINE__);
  st2.Insert(index,"// CKB comment out, handled explicitly by the train macro \n //");


  ofstream out2;
  out2.open(Form("%s_merge.C", kTrainName.Data()));
  if (out2.bad()) {
    return kFALSE;
  }
  out2 << st2 << endl;
  Printf("<<< Patching Merging Macro");


  return kTRUE;

}
