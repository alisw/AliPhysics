//===================== ANALYSIS TRAIN =========================================
// To use: copy this macro to your work directory, modify the global part to match
// your needs, then run root.
//    root[0] .L AnalysisTrain.C
// Grid full mode as below (other modes: test, offline, submit, terminate)
//    root[1] AnalysisTrainNew("grid", "full")
// CAF mode (requires root v5-23-02 + aliroot v4-16-Rev08)
//    root[2] AnalysisTrainNew("proof")
// Local mode requires AliESds.root or AliAOD.root in ./data directory
//    root[3] AnalysisTrainNew("local")
// In proof and grid modes, a token is needed and sourcing the produced environment file.
//
// If 'saveTrain' flag is set, the train will generate a directory name and run
// in this directory. A configuration file 'ConfigTrain.C' will be generated. 
// One can replay at any time the train via:
//    root[1] AnalysisTrainNew(ana_mode, plugin_mode, "train_default_<date>/ConfigTrain.C")

//==================   TRAIN NAME   ============================================
TString     train_name         = "default"; // enters file names, so no blancs or special characters
//==============================================================================

// ### Settings that make sense in PROOF only
//==============================================================================
TString     proof_cluster      = "alicecaf.cern.ch";
Bool_t      useAFPAR           = kFALSE;  // use AF special par file
TString     AFversion          = "AF-v4-16";
// Change CAF dataset here
TString     proof_dataset      = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";
TString     proof_outdir       = "";

// ### Settings that make sense when using the Alien plugin
//==============================================================================
Bool_t      usePLUGIN          = kTRUE;   // do not change
// Usage of par files ONLY in grid mode and ONLY if the code is not available
// in the deployed AliRoot versions. Par file search path: local dir, if not there $ALICE_ROOT.
// To refresh par files, remove the ones in the workdir, then do "make <target.par>" in 
// AliRoot.
Bool_t      usePAR             = kFALSE;  // use par files for extra libs
Bool_t      useCPAR            = kFALSE;  // use par files for common libs
TString     root_version       = "v5-23-04";
TString     aliroot_version    = "v4-17-01";
// Change production base directory here
TString     alien_datadir      = "/alice/sim/PDC_09/LHC09a4/";
// Use up to 10 non-zero run numbers
Int_t       run_numbers[10]    = {81272,    81273 ,     81274,     0,     0,
                                      0,     0,     0,     0,     0};
// ### Settings that make sense only for local analysis
//==============================================================================
// Change local xml dataset for local interactive analysis
TString     local_xmldataset   = "";

// ### Other flags to steer the analysis
//==============================================================================
Bool_t      useDBG             = kTRUE;  // activate debugging
Bool_t      useMC              = kTRUE;  // use MC info
Bool_t      useTAGS            = kFALSE; // use ESD tags for selection
Bool_t      useKFILTER         = kTRUE;  // use Kinematics filter
Bool_t      useTR              = kFALSE; // use track references
Bool_t      useCORRFW          = kFALSE; // do not change
Bool_t      useAODTAGS         = kTRUE;  // use AOD tags
Bool_t      saveTrain          = kTRUE;  // save train configuration as: 
Bool_t      saveProofToAlien   = kFALSE; // save proof outputs in AliEn
                                         // train_[trainName]_ddMonthyyyy_time.C
// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iAODanalysis       = 0;      // Analysis on input AOD's
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iESDfilter         = 1;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iMUONcopyAOD       = 0;      // Task that copies only muon events in a separate AOD (PWG3)
Int_t       iJETAN             = 0;      // Jet analysis (PWG4) - needs ESD filter
Int_t       iPWG4partcorr      = 1;      // Gamma-hadron correlations task (PWG4)
Int_t       iPWG3vertexing     = 0;      // Vertexing HF task (PWG2)
Int_t       iPWG2femto         = 0;      // Femtoscopy task (PWG2)
Int_t       iPWG2spectra       = 0;      // Spectra PWG2 tasks (protons, cascades, V0 check, strange)
Int_t       iPWG2flow          = 0;      // Flow analysis task (PWG2)
Int_t       iPWG2res           = 0;      // Resonances task (PWG2)
Int_t       iPWG2kink          = 0;      // Kink analysis task (PWG2)

// Temporaries.
TString anaPars = "";
TString anaLibs = "";
// Function signatures
class AliAnalysisGrid;

//______________________________________________________________________________
void AnalysisTrainNew(const char *analysis_mode="grid", 
                      const char *plugin_mode="full",
                      const char *config_file="")
{
// Main analysis train macro. If a configuration file is provided, all parameters
// are taken from there but may be altered by CheckModuleFlags.
   if (strlen(config_file) && !LoadConfig(config_file)) return;
   TString smode(analysis_mode);
   smode.ToUpper();
   if (saveTrain)              WriteConfig();
   // Check compatibility of selected modules
   CheckModuleFlags(smode);

   printf("==================================================================\n");
   printf("===========    RUNNING ANALYSIS TRAIN %s IN %s MODE   ==========\n", train_name.Data(),smode.Data());
   printf("==================================================================\n");
   printf("=  Configuring analysis train for:                               =\n");
   if (iAODanalysis) printf("=  AOD analysis                                                  =\n");
   else              printf("=  ESD analysis                                                  =\n");
   if (iESDfilter)   printf("=  ESD filter                                                    =\n");
   if (iMUONcopyAOD) printf("=  MUON copy AOD                                                 =\n");
   if (iJETAN)       printf("=  Jet analysis                                                  =\n");
   if (iPWG2spectra) printf("=  PWG2 proton, checkCascade, checkV0, strange                   =\n");
   if (iPWG2femto)   printf("=  PWG2 femtoscopy                                               =\n");
   if (iPWG2flow)    printf("=  PWG2 flow                                                     =\n");
   if (iPWG2res)     printf("=  PWG2 resonances                                               =\n");
   if (iPWG2kink)    printf("=  PWG2 kink analysis                                            =\n");
   if (iPWG3vertexing) printf("=  PWG3 vertexing                                            =\n");
   if (iPWG4partcorr)  printf("=  PWG4 gamma-hadron, pi0 and gamma-jet correlations         =\n");
   printf("==================================================================\n");
   printf(":: use MC truth      %d\n", (UInt_t)useMC);
   printf(":: use KINE filter   %d\n", (UInt_t)useKFILTER);
   printf(":: use track refs    %d\n", (UInt_t)useTR);
   printf(":: use tags          %d\n", (UInt_t)useTAGS);
   printf(":: use AOD tags      %d\n", (UInt_t)useAODTAGS);
   printf(":: use debugging     %d\n", (UInt_t)useDBG);
   printf(":: use PAR files     %d\n", (UInt_t)usePAR);
   printf(":: use AliEn plugin  %d\n", (UInt_t)usePLUGIN);

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
   AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");
   if (saveProofToAlien) mgr->SetSpecialOutputLocation(proof_outdir);

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
   } else {   
   // ESD input handler
      AliESDInputHandler *esdHandler = new AliESDInputHandler();
      if (useTAGS) esdHandler->SetReadTags();
      mgr->SetInputEventHandler(esdHandler);       
   }
   // Monte Carlo handler
   if (useMC && !iAODanalysis) {
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
      mcHandler->SetReadTR(useTR); 
   }   
   // AOD output container, created automatically when setting an AOD handler
   if (iAODhandler) {
      // AOD output handler
      AliAODHandler* aodHandler   = new AliAODHandler();
      aodHandler->SetOutputFileName("AliAOD.root");
      mgr->SetOutputEventHandler(aodHandler);
      if (iAODanalysis) {
         aodHandler->SetCreateNonStandardAOD();
         if (iPWG3vertexing) aodHandler->SetOutputFileName("AliAOD.VertexingHF.root");
      } 
   }
   // Debugging if needed
   if (useDBG) mgr->SetDebugLevel(3);

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

   if (iESDfilter && !iAODanalysis) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
      AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(useKFILTER);
   }   

   // AOD tags
   if (useAODTAGS) {
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
      AliAnalysisTaskJets *taskjets = AddTaskJets("AOD", "UA1");
      if (!taskjets) ::Warning("AnalysisTrainNew", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
   }
       
   // Proton analysis
   if (iPWG2spectra) {
      // protons
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/AddTaskProtons.C");
      AliAnalysisTaskProtons *taskprotons = AddTaskProtons();
      if (!taskprotons) ::Warning("AnalysisTrainNew", "AliAnalysisTaskProtons cannot run for this train conditions - EXCLUDED");
      // cascades
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/AddTaskCheckCascade.C");
      AliAnalysisTaskCheckCascade *taskcheckcascade = AddTaskCheckCascade();      
      if (!taskcheckcascade) ::Warning("AnalysisTrainNew", "AliAnalysisTaskCheckCascade cannot run for this train conditions - EXCLUDED");
      // v0's
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/AddTaskCheckV0.C");
      AliAnalysisTaskCheckV0 *taskcheckV0 = AddTaskCheckV0();
      if (!taskcheckV0) ::Warning("AnalysisTrainNew", "AliAnalysisTaskCheckV0 cannot run for this train conditions - EXCLUDED");
      // strangeness
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/AddTaskStrange.C");
      AliAnalysisTaskStrange *taskstrange = AddTaskStrange();
      if (!taskstrange) ::Warning("AnalysisTrainNew", "AliAnalysisTaskStrange cannot run for this train conditions - EXCLUDED");
   }   
   
   // Femtoscopy analysis modules
   if (iPWG2femto) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/FEMTOSCOPY/macros/AddTaskFemto.C");
      AliAnalysisTaskFemto *taskfemto = AddTaskFemto();
      if (!taskfemto) ::Warning("AnalysisTrainNew", "AliAnalysisTaskFemto cannot run for this train conditions - EXCLUDED");
   }   

   // Kink analysis
   if (iPWG2kink) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/KINK/macros/AddTaskKink.C");
      AliAnalysisKinkESDMC *taskkink = AddTaskKink();
      if (!taskkink) ::Warning("AnalysisTrainNew", "AliAnalysisKinkESDMC cannot run for this train conditions - EXCLUDED");
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/KINK/macros/AddTaskKinkResonance.C");
      AliResonanceKinkPID *taskkinkres = AddTaskKinkResonance();
      if (!taskkinkres) ::Warning("AnalysisTrainNew", "AliResonanceKinkPID cannot run for this train conditions - EXCLUDED");
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/KINK/macros/AddTaskKinkResonanceLikeSign.C");
      AliResonanceKinkLikeSign *taskkinklikesign = AddTaskKinkResonanceLikeSign();
      if (!taskkinklikesign) ::Warning("AnalysisTrainNew", "AliResonanceKinkLikeSign cannot run for this train conditions - EXCLUDED");
   }   
      
  // PWG3 vertexing
   if (iPWG3vertexing) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskVertexingHF.C");
      AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
      if (!taskvertexingHF) ::Warning("AnalysisTrainNew", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
   }   
      
   // PWG4 hadron correlations
   if (iPWG4partcorr) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
      AliAnalysisTaskParticleCorrelation *taskpartcorrPHOS = AddTaskPartCorr("AOD", "PHOS");
      if (!taskpartcorrPHOS) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation PHOS cannot run for this train conditions - EXCLUDED");
      AliAnalysisTaskParticleCorrelation *taskpartcorrEMCAL = AddTaskPartCorr("AOD", "EMCAL");
      if (!taskpartcorrEMCAL) ::Warning("AnalysisTrainNew", "AliAnalysisTaskParticleCorrelation EMCAL cannot run for this train conditions - EXCLUDED");
   }   
   //==========================================================================
   // FOR THE REST OF THE TASKS THE MACRO AddTaskXXX() is not yet implemented/
   // Run the analysis
   //    
   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      if (saveTrain) gSystem->ChangeDirectory(train_name);
      StartAnalysis(smode, chain);
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
            ::Error("AnalysisTrainNew.C::StartAnalysis", "Cannot create the chain");
            return;
         }   
         mgr->StartAnalysis(mode, chain);
         return;
      case 1:
         if (!proof_dataset.Length()) {
            ::Error("AnalysisTrainNew.C::StartAnalysis", "proof_dataset is empty");
            return;
         }   
         mgr->StartAnalysis(mode, proof_dataset, 1000);
         return;
      case 2:
         if (usePLUGIN) {
            if (!mgr->GetGridHandler()) {
               ::Error("AnalysisTrainNew.C::StartAnalysis", "Grid plugin not initialized");
               return;
            }   
            mgr->StartAnalysis("grid");
         } else {
            if (!chain) {
               ::Error("AnalysisTrainNew.C::StartAnalysis", "Cannot create the chain");
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
   if (imode==1) {
      if (!usePAR) {
         ::Info("AnalysisTrainNew.C::CheckModuleFlags", "PAR files enabled due to PROOF analysis");
         usePAR = kTRUE;
      }   
   }  
   if (imode != 2) {
      ::Info("AnalysisTrainNew.C::CheckModuleFlags", "AliEn plugin disabled since not in GRID mode");
      usePLUGIN = kFALSE; 
   }   
   if (iAODanalysis) {
   // AOD analysis
      if (useMC)
         ::Info("AnalysisTrainNew.C::CheckModuleFlags", "MC usage disabled in analysis on AOD's");
      useMC = kFALSE;
      useTR = kFALSE;
      if (iESDfilter)
         ::Info("AnalysisTrainNew.C::CheckModuleFlags", "ESD filter disabled in analysis on AOD's");
      iESDfilter   = 0;
      if (!iAODhandler) {
         if (iJETAN) 
            ::Info("AnalysisTrainNew.C::CheckModuleFlags", "JETAN disabled in analysis on AOD's without AOD handler");
         iJETAN = 0;
      }
      // Disable tasks that do not work yet on AOD data
      if (iPWG2kink)         
         ::Info("AnalysisTrainNew.C::CheckModuleFlags", "PWG2kink disabled in analysis on AOD's");
         iPWG2kink = 0;
   } else {   
   // ESD analysis
      iMUONcopyAOD = 0;
   }       
   if (iJETAN) iESDfilter=1;
   if (iESDfilter) iAODhandler=1;
   if (iPWG2spectra || iPWG2flow || iPWG3vertexing) useCORRFW = kTRUE;
   if (useKFILTER && !useMC) useKFILTER = kFALSE;
   if (useAODTAGS && !iAODhandler) useAODTAGS = kFALSE;
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
            ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         ::Info("AnalysisTrainNew.C::Connect", "Connecting user <%s> to PROOF cluster <%s>", 
                username.Data(), proof_cluster.Data());
         gEnv->SetValue("XSec.GSI.DelegProxy", "2");
//         TProof::Open(Form("%s@%s:31093", username.Data(), proof_cluster.Data()));       
         TProof::Open(Form("%s@%s", username.Data(), proof_cluster.Data()));       
         if (!gProof) {
            if (strcmp(gSystem->Getenv("XrdSecGSISRVNAMES"), "lxfsrd0506.cern.ch"))
               ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), "Environment XrdSecGSISRVNAMES different from lxfsrd0506.cern.ch");
            return kFALSE;
         }
         TGrid::Connect("alien://");
         if (gGrid) {
            TString homedir = gGrid->GetHomeDirectory();
            TString workdir = homedir + train_name;
            if (!gGrid->Cd(workdir)) {
               gGrid->Cd(homedir);
               if (gGrid->Mkdir(workdir)) {
                  gGrid->Cd(train_name);
                  ::Info("AnalysisTrainNew::Connect()", "Directory %s created", gGrid->Pwd());
               }
            }
            gGrid->Mkdir("proof_output");
            gGrid->Cd("proof_output");
            proof_outdir = Form("alien://%s", gGrid->Pwd());
         }   
         break;
      case 2:      
         if  (!username.Length()) {
            ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         if (usePLUGIN && !gSystem->Getenv("alien_CLOSE_SE")) {
            ::Error(Form("AnalysisTrainNew.C::Connect <%s>", mode), 
                           "When using the AliEn plugin it is preferable to define the \
                           variable alien_CLOSE_SE in your environment.");
            return kFALSE;
         }
         ::Info("AnalysisTrainNew.C::Connect", "Connecting user <%s> to AliEn ...", 
                username.Data());
         TGrid::Connect("alien://");
         if (!gGrid || !gGrid->IsConnected()) return kFALSE;
         break;
      default:
         ::Error("AnalysisTrainNew.C::Connect", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   ::Info("AnalysisTrainNew.C::Connect","Connected in %s mode", mode);
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
      ::Error("AnalysisTrainNew.C::LoadCommonLibraries", "Analysis train requires that analysis libraries are compiled with a local AliRoot"); 
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
         if (useCPAR) {
            success &= LoadLibrary("STEERBase", mode, kTRUE);
            success &= LoadLibrary("ESD", mode, kTRUE);
            success &= LoadLibrary("AOD", mode, kTRUE);
            success &= LoadLibrary("ANALYSIS", mode, kTRUE);
            success &= LoadLibrary("ANALYSISalice", mode, kTRUE);
            if (useCORRFW) success &= LoadLibrary("CORRFW", mode, kTRUE);
         } else {   
            success &= LoadLibrary("libSTEERBase.so", mode);
            success &= LoadLibrary("libESD.so", mode);
            success &= LoadLibrary("libAOD.so", mode);
            success &= LoadLibrary("libANALYSIS.so", mode);
            success &= LoadLibrary("libANALYSISalice.so", mode);
            if (useCORRFW) success &= LoadLibrary("libCORRFW.so", mode);
            gROOT->ProcessLine(".include $ALICE_ROOT/include");
         }   
         break;
      case 1:
         Int_t ires = -1;
         if (useAFPAR && !gSystem->AccessPathName(AFversion)) ires = gProof->UploadPackage(AFversion);
         if (ires < 0) {
            success &= LoadLibrary("STEERBase", mode);
            success &= LoadLibrary("ESD", mode);
            success &= LoadLibrary("AOD", mode);
            success &= LoadLibrary("ANALYSIS", mode);
            success &= LoadLibrary("ANALYSISalice", mode);
            if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
         } else { 
            ires = gProof->EnablePackage(AFversion);
            if (ires<0) success = kFALSE;
            if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
         }
         break;         
      default:
         ::Error("AnalysisTrainNew.C::LoadCommonLibraries", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   if (success) {
      ::Info("AnalysisTrainNew.C::LoadCommodLibraries", "Load common libraries:    SUCCESS");
      ::Info("AnalysisTrainNew.C::LoadCommodLibraries", "Include path for Aclic compilation:\n%s",
              gSystem->GetIncludePath());
   } else {           
      ::Info("AnalysisTrainNew.C::LoadCommodLibraries", "Load common libraries:    FAILED");
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
   if (iJETAN) {
      if (!LoadLibrary("JETAN", mode, kTRUE)) return kFALSE;
   }   
            
   // PWG4 particle correlations
   if (iPWG4partcorr) {   
      if (!LoadLibrary("PWG4PartCorrBase", mode, kTRUE) ||
          !LoadLibrary("PWG4PartCorrDep", mode, kTRUE)) return kFALSE;
   }
   // PWG2 task protons 
   if (iPWG2spectra) {
      if (!LoadLibrary("PWG2spectra", mode, kTRUE)) return kFALSE;
   }
   // PWG2 flow
   if (iPWG2flow) {
      if (!LoadLibrary("PWG2AOD", mode, kTRUE) ||
          !LoadLibrary("PWG2flow", mode, kTRUE)) return kFALSE;
   }
   // PWG2 resonances
   if (iPWG2res) {
      if (!LoadLibrary("PWG2resonances", mode, kTRUE)) return kFALSE;
   }   
   // PWG2 kink
   if (iPWG2kink) {
      if (!LoadLibrary("PWG2kink", mode, kTRUE)) return kFALSE;
   }   
   // PWG2 femtoscopy
   if (iPWG2femto) {
      if (!LoadLibrary("PWG2AOD", mode, kTRUE) ||
          !LoadLibrary("PWG2femtoscopy", mode, kTRUE) ||
          !LoadLibrary("PWG2femtoscopyUser", mode, kTRUE)) return kFALSE;
      TFile::Cp(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG2/FEMTOSCOPY/macros/ConfigFemtoAnalysis.C"), Form("%s/ConfigFemtoAnalysis.C", train_name.Data()));
      anaLibs += "ConfigFemtoAnalysis.C ";
   }   
   // Vertexing HF
   if (iPWG3vertexing) {
      if (!LoadLibrary("PWG3base", mode, kTRUE) ||
          !LoadLibrary("PWG3vertexingHF", mode, kTRUE)) return kFALSE;
   }   
   ::Info("AnalysisTrainNew.C::LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
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
      ::Error("AnalysisTrainNew.C::LoadLibrary", "Empty module name");
      return kFALSE;
   }   
   // If a library is specified, just load it
   if (smodule.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      result = gSystem->Load(mod);
      if (result < 0) {
         ::Error("AnalysisTrainNew.C::LoadLibrary", "Could not load library %s", module);
         return kFALSE;
      }
      if (rec) anaLibs += Form("%s.so ",mod.Data()); 
      return kTRUE;
   } 
   // Check if the library is already loaded
   if (strlen(gSystem->GetLibraries(Form("%s.so", module), "", kFALSE)) > 0)
      return kTRUE;    
   switch (imode) {
      case 0:
      case 2:
         if (usePAR) {
            result = SetupPar(module);
            if (rec) anaPars += Form("%s.par ", module);
         } else {
            result = gSystem->Load(Form("lib%s.so", module));
            if (rec) anaLibs += Form("lib%s.so ", module);
         }   
         break;
      case 1:
         result = gProof->UploadPackage(module);
         if (result<0) {
            result = gProof->UploadPackage(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", module)));
            if (result<0) {
               ::Error("AnalysisTrainNew.C::LoadLibrary", "Could not find module %s.par in current directory nor in $ALICE_ROOT", module);
               return kFALSE;
            }
         }   
         result = gProof->EnablePackage(module);
         break;
      default:
         return kFALSE;
   }         
   if (result < 0) {
      ::Error("AnalysisTrainNew.C::LoadLibrary", "Could not load module %s", module);
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
            if (!local_xmldataset.Length()) {
               // Local AOD
               chain = new TChain("aodTree");
               if (gSystem->AccessPathName("data/AliAOD.root")) 
                  ::Error("AnalysisTrainNew.C::CreateChain", "File: AliAOD.root not in ./data dir");
               else {
                  if (!saveTrain) chain->Add("data/AliAOD.root");
                  else            chain->Add("../data/AliAOD.root");
               }   
            } else {
               // Interactive AOD
               chain = CreateChainSingle(local_xmldataset, "aodTree");
            }
         } else {      
            if (!local_xmldataset.Length()) {
               // Local ESD
               chain = new TChain("esdTree");
               if (gSystem->AccessPathName("data/AliESDs.root")) 
                  ::Error("AnalysisTrainNew.C::CreateChain", "File: AliESDs.root not in ./data dir");
               else {
                  if (!saveTrain) chain->Add("data/AliESDs.root");
                  else            chain->Add("../data/AliESDs.root");
               }   
            } else {
               // Interactive ESD
               chain = CreateChainSingle(local_xmldataset, "esdTree");
            }   
         }
         break;
      case 1:
         break;
      case 2:
         if (usePLUGIN) {
            AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
            AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
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
      ::Error("AnalysisTrainNew.C::CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
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
         ::Info("AnalysisTrainNew.C::SetupPar", "Getting %s.par from $ALICE_ROOT", pararchivename);
         TFile::Cp(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", pararchivename)), 
                   Form("%s.par",pararchivename));
      } else {
         ::Error("AnalysisTrainNew.C::SetupPar", "Cannot find %s.par", pararchivename);
         return -1;
      }   
   }
   if (usePLUGIN && saveTrain) gSystem->Exec(Form("ln -s ../%s.par %s",pararchivename, train_name.Data()));
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
AliAnalysisGrid* CreateAlienHandler(const char *plugin_mode)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(plugin_mode);
   plugin->SetNtestFiles(1);
   plugin->SetPreferedSE("ALICE::NIHAM::FILE");
// Set versions of used packages
   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion(root_version);
   plugin->SetAliROOTVersion(aliroot_version);
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir(alien_datadir);
// Set data search pattern
   if (iAODanalysis) plugin->SetDataPattern("*AliAOD.root");
   else              plugin->SetDataPattern("*AliESDs.root");
// ...then add run numbers to be considered
   for (Int_t i=0; i<10; i++) {
      if (run_numbers[i]==0) break;
      plugin->AddRunNumber(run_numbers[i]);
   }   
// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   if (iAODanalysis) plugin->SetGridWorkingDir("analysisAOD");
   else              plugin->SetGridWorkingDir("analysisESD");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir(Form("output_%s",train_name.Data())); // In this case will be $HOME/work/output

   TString ana_sources = "";
   TString ana_add = "";
   if (usePAR && anaPars.Length()) {
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
   ana_sources = ana_sources.Strip();
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   anaLibs     = anaLibs.Strip();   
   if (ana_sources.Length()) plugin->SetAnalysisSource(ana_sources);
   if (anaLibs.Length())     plugin->SetAdditionalLibs(anaLibs);
     
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs();
   plugin->SetMergeExcludes("AliAOD.root");
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(100);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName(Form("%s.jdl", train_name.Data()));
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable(Form("%s.sh", train_name.Data()));
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
   gSystem->Exec("date +%d%b%Y_%Hh%M > date.tmp");
   ifstream fdate("date.tmp");
   if (!fdate.is_open()) {
      ::Error("AnalysisTrainNew.C::Export","Could not generate file name");
      return;
   }
   const char date[64];
   fdate.getline(date,64);
   fdate.close();
   gSystem->Exec("rm date.tmp");
   train_name = Form("train_%s_%s", train_name.Data(), date);
   TString cdir = gSystem->WorkingDirectory();
   gSystem->MakeDirectory(train_name);
   gSystem->ChangeDirectory(train_name);
   ofstream out;
   out.open("ConfigTrain.C", ios::out); 
   if (out.bad()) {
      ::Error("AnalysisTrainNew.C::Export", "Cannot open ConfigTrain.C for writing");
      return;
   }
   out << "{" << endl;
   out << "   train_name      = " << "\"" << train_name.Data() << "\";" << endl;
   out << "   proof_cluster   = " << "\"" << proof_cluster.Data() << "\";" << endl;
   out << "   useAFPAR        = " << useAFPAR << ";" << endl;
   if (useAFPAR) 
      out << "   AFversion       = " << AFversion.Data() << ";" << endl;
   out << "   proof_dataset   = " << "\"" << proof_dataset.Data() << "\";" << endl;
   out << "   usePLUGIN       = " << usePLUGIN << ";" << endl;
   out << "   usePAR          = " << usePAR << ";" << endl;
   out << "   useCPAR         = " << useCPAR << ";" << endl;
   out << "   root_version    = " << "\"" << root_version.Data() << "\";" << endl;
   out << "   aliroot_version = " << "\"" << aliroot_version.Data() << "\";" << endl;
   out << "   alien_datadir   = " << "\"" << alien_datadir.Data() << "\";" << endl;
   for (Int_t i=0; i<10; i++) {
      if (run_numbers[i]) 
         out << "   run_numbers[" << i << "]  = " << run_numbers[i] << ";" << endl;
   }
   out << "   useDBG          = " << useDBG << ";" << endl;
   out << "   useMC           = " << useMC << ";" << endl;
   out << "   useTAGS         = " << useTAGS << ";" << endl;
   out << "   useKFILTER      = " << useKFILTER << ";" << endl;
   out << "   useTR           = " << useTR << ";" << endl;
   out << "   useCORRFW       = " << useCORRFW << ";" << endl;
   out << "   useAODTAGS      = " << useAODTAGS << ";" << endl;
   out << "   saveTrain       = " << "kFALSE;" << endl << endl;
   out << "   // Analysis modules" << endl;
   out << "   iAODanalysis    = " << iAODanalysis << ";" << endl;
   out << "   iAODhandler     = " << iAODhandler << ";" << endl;
   out << "   iESDfilter      = " << iESDfilter << ";" << endl;
   out << "   iMUONcopyAOD    = " << iMUONcopyAOD << ";" << endl;
   out << "   iJETAN          = " << iJETAN << ";" << endl;
   out << "   iPWG4partcorr   = " << iPWG4partcorr << ";" << endl;
   out << "   iPWG2femto      = " << iPWG2femto << ";" << endl;
   out << "   iPWG2spectra    = " << iPWG2spectra << ";" << endl;
   out << "   iPWG2flow       = " << iPWG2flow << ";" << endl;
   out << "   iPWG2res        = " << iPWG2res << ";" << endl;
   out << "   iPWG2kink       = " << iPWG2kink << ";" << endl;
   out << "}" << endl;
   ::Info("AnalysisTrainNew.C::WriteConfig", "Train configuration wrote to file %s", Form("config_%s.C", train_name.Data()));
   gSystem->ChangeDirectory(cdir);
}   

//______________________________________________________________________________
Bool_t LoadConfig(const char *filename)
{
// Read train configuration from file
   if (gSystem->AccessPathName(filename)) {
      ::Error("AnalysisTrainNew.C::LoadConfig", "Config file name not found");
      return kFALSE;
   }   
   gROOT->ProcessLine(Form(".x %s", filename));
   ::Info("AnalysisTrainNew.C::LoadConfig", "Train configuration loaded from file %s", filename);
   return kTRUE;
}
