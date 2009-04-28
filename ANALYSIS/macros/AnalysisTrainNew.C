//===================== ANALYSIS TRAIN ==========================================
// To use: copy this macro to your work directory, modify the global part to match
// your needs, then run root.
//    root[0] .L AnalysisTrain.C
// Grid full mode as below (other modes: test, offline, submit, terminate)
//    root[1] AnalysisTrainNew("grid", "full")
// CAF mode (requires root v5-23-02 + aliroot v4-16-Rev08)
//    root[2] AnalysisTrainNew("proof")
// Local mode requires AliESds.root or AliAOD.root in ./data directory
//    root[3] AnalysisTrainNew("local")

const char *root_version    = "v5-23-02";
const char *aliroot_version = "v4-16-Rev-08";
const char *cluster         = "alicecaf.cern.ch";
//const char *AFversion       = "AF-v4-16";
const char *AFversion       = "";
// Dataset name. If empty assumes local ESD/AOD. In CAF it is a dataset name.
//               In grid has to point to an xml file or be empty if using the plugin
// === PROOF
const char *proof_dataset = "/COMMON/COMMON/LHC09a4_run8100X#/esdTree";

// === ALIEN
const char *dataset   = "";
const char *alien_datadir = "/alice/sim/PDC_09/LHC09a4/";
Int_t run_numbers[10] = {81072,     0,     0,     0,     0,
                             0,     0,     0,     0,     0};

Bool_t useDBG        = kTRUE;
Bool_t useMC         = kTRUE;
Bool_t useTAGS       = kFALSE;
Bool_t useKFILTER    = kTRUE;
Bool_t usePAR        = kFALSE;
Bool_t useTR         = kFALSE;
Bool_t usePLUGIN     = kTRUE;
Bool_t useCORRFW     = kFALSE; // do not change
Bool_t useAODTAGS    = kTRUE;
    
Int_t iAODanalysis   = 0;
Int_t iAODhandler    = 1;
Int_t iESDfilter     = 1;
Int_t iMUONcopyAOD   = 0;
Int_t iJETAN         = 1;
Int_t iPWG4partcorr  = 1;
Int_t iPWG4pi0       = 1;
Int_t iPWG4gammajet  = 1;
Int_t iPWG3vertexing = 1;
Int_t iPWG2femto     = 0;
Int_t iPWG2spectra   = 1;
Int_t iPWG2flow      = 0;
Int_t iPWG2res       = 1;

TString anaPars = "";
TString anaLibs = "";
// Function signatures
class AliAnalysisGrid;

//______________________________________________________________________________
void AnalysisTrainNew(const char *analysis_mode="grid", const char *plugin_mode="test")
{
// Example of running analysis train
   TString smode(analysis_mode);
   smode.ToUpper();    
   // Check compatibility of selected modules
   CheckModuleFlags(smode);

   printf("==================================================================\n");
   printf("===========    RUNNING ANALYSIS TRAIN IN %s MODE   ==========\n", smode.Data());
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
   if (iPWG3vertexing) printf("=  PWG3 vertexing                                            =\n");
   if (iPWG4partcorr)  printf("=  PWG4 gamma-hadron correlations                            =\n");
   if (iPWG4pi0)     printf("=  PWG4 pi0 analysis                                             =\n");
   if (iPWG4gammajet)  printf("=  PWG4 gamma jet analysis                                       =\n");
   printf("==================================================================\n");
   printf(":: use MC truth      %d\n", (UInt_t)useMC);
   printf(":: use KINE filter   %d\n", (UInt_t)useKFILTER);
   printf(":: use track refs    %d\n", (UInt_t)useTR);
   printf(":: use tags          %d\n", (UInt_t)useTAGS);
   printf(":: use AOD tags      %d\n", (UInt_t)useAODTAGS);
   printf(":: use debugging     %d\n", (UInt_t)useDBG);
   printf(":: use PAR files     %d\n", (UInt_t)usePAR);
   printf(":: use AliEn plugin  %d\n", (UInt_t)usePLUGIN);

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
    
   // Load analysis specific libraries
   if (!LoadAnalysisLibraries(smode)) {
      ::Error("AnalysisTrain", "Could not load analysis libraries");
      return;
   }   

    
   //==========================================================================
   // Make the analysis manager and connect event handlers
   AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");

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
      mgr->SetOutputEventHandler(aodHandler);
      aodHandler->SetOutputFileName("AliAOD.root");
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
      AliAnalysisTaskJets *taskjets = AddTaskJets("ESD", "UA1");
   }
       
   // Proton analysis
   if (iPWG2spectra) {
      // protons
      gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/AddTaskProtons.C");
      AliAnalysisTaskProtons *taskprotons = AddTaskProtons();
      // cascades
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCheckCascade.C");
      AliAnalysisTaskCheckCascade *taskcheckcascade = AddTaskCheckCascade();      
      // v0's
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCheckV0.C");
      AliAnalysisTaskCheckV0 *taskcheckV0 = AddTaskCheckV0();
      // strangeness
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskStrange.C");
      AliAnalysisTaskStrange *taskstrange = AddTaskStrange();
   }   
   
   // PWG3 vertexing
   if (iPWG3vertexing) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/AddTaskVertexingHF.C");
      AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
   }   
      
   // PWG4 hadron correlations
   if (iPWG4partcorr) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
      AliAnalysisTaskParticleCorrelation *taskgammahadron = AddTaskPartCorr("GammaHadron", "AOD", "PHOS");
   }   
       
   // PWG4 pi0
   if (iPWG4pi0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
      AliAnalysisTaskParticleCorrelation *taskpi0 = AddTaskPartCorr("Pi0", "AOD", "PHOS");
   }   
    
   // PWG4 gamma jet finder
   if (iPWG4gammajet) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskPartCorr.C");
      AliAnalysisTaskParticleCorrelation *taskgammajet = AddTaskPartCorr("GammaJetFinder", "AOD", "PHOS");
   }   
   //==========================================================================
   // FOR THE REST OF THE TASKS THE MACRO AddTaskXXX() is not yet implemented/
   // Run the analysis
   //    
   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
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
            ::Error("StartAnalysis", "Cannot create the chain");
            return;
         }   
         mgr->StartAnalysis(mode, chain);
         return;
      case 1:
         if (!strlen(proof_dataset)) {
            ::Error("StartAnalysis", "proof_dataset is empty");
            return;
         }   
         mgr->StartAnalysis(mode, proof_dataset, 1000);
         return;
      case 2:
         if (usePLUGIN) {
            if (!mgr->GetGridHandler()) {
               ::Error("StartAnalysis", "Grid plugin not initialized");
               return;
            }   
            mgr->StartAnalysis("grid");
         } else {
            if (!chain) {
               ::Error("StartAnalysis", "Cannot create the chain");
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
      usePAR = kTRUE;
   }   
   if (iAODanalysis) {
   // AOD analysis
      useMC = kFALSE;
      useTR = kFALSE;
      iESDfilter   = 0;
      // Disable tasks that do not work yet on AOD data
   } else {   
   // ESD analysis
      iMUONcopyAOD = 0;
//      iPWG3vertexing = 0;
   }       
   if (iJETAN) iESDfilter=1;
   if (iESDfilter) iAODhandler=1;
   if (iAODanalysis || !iAODhandler) {
      iPWG4partcorr=0;
      iPWG4pi0=0;
   }   
   if (iPWG4gammajet && !iJETAN) {
      ::Error("CheckModuleFlags", "Gamma jet finder needs JETAN. Disabling.");
      iPWG4gammajet = 0;
   }
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
            ::Error(Form("Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         ::Info("Connect", "Connecting user <%s> to PROOF cluster <%s>", 
                username.Data(), cluster);
         TProof::Open(Form("%s@%s", username.Data(), cluster));       
         if (!gProof) {
            if (strcmp(gSystem->Getenv("XrdSecGSISRVNAMES"), "lxfsrd0506.cern.ch"))
               ::Error(Form("Connect <%s>", mode), "Environment XrdSecGSISRVNAMES different from lxfsrd0506.cern.ch");
            return kFALSE;
         }
         break;
      case 2:      
         if  (!username.Length()) {
            ::Error(Form("Connect <%s>", mode), "Make sure you:\n \
                           1. Have called: alien-token-init <username>\n \
                           2. Have called: >source /tmp/gclient_env_$UID");
            return kFALSE;
         }
         if (usePLUGIN && !gSystem->Getenv("alien_CLOSE_SE")) {
            ::Error(Form("Connect <%s>", mode), 
                           "When using the AliEn plugin it is preferable to define the \
                           variable alien_CLOSE_SE in your environment.");
            return kFALSE;
         }
         ::Info("Connect", "Connecting user <%s> to AliEn ...", 
                username.Data());
         TGrid::Connect("alien://");
         if (!gGrid || !gGrid->IsConnected()) return kFALSE;
         break;
      default:
         ::Error("Connect", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   ::Info("Connect","Connected in %s mode", mode);
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
      ::Error("LoadCommonLibraries", "Analysis train requires that analysis libraries are compiled with a local AliRoot"); 
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
         success &= LoadLibrary("libSTEERBase.so", mode);
         success &= LoadLibrary("libESD.so", mode);
         success &= LoadLibrary("libAOD.so", mode);
         success &= LoadLibrary("libANALYSIS.so", mode);
         success &= LoadLibrary("libANALYSISalice.so", mode);
         if (useCORRFW) success &= LoadLibrary("libCORRFW.so", mode);
         gROOT->ProcessLine(".include $ALICE_ROOT/include");
         break;
      case 1:
         Int_t ires = -1;
         if (!gSystem->AccessPathName(AFversion)) ires = gProof->UploadPackage(AFversion);
         if (ires < 0) {
            success &= LoadLibrary("STEERBase", mode);
            success &= LoadLibrary("ESD", mode);
            success &= LoadLibrary("AOD", mode);
            success &= LoadLibrary("ANALYSIS", mode);
            success &= LoadLibrary("ANALYSISalice", mode);
            if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
         } else { 
            ires = gProof->EnablePackage(AFversion);
            if (useCORRFW) success &= LoadLibrary("CORRFW", mode);
         }
         if (ires<0) success = kFALSE;
         break;         
      default:
         ::Error("LoadCommonLibraries", "Unknown run mode: %s", mode);
         return kFALSE;
   }
   if (success) {
      ::Info("LoadCommodLibraries", "Load common libraries:    SUCCESS");
      ::Info("LoadCommodLibraries", "Include path for Aclic compilation:\n%s",
              gSystem->GetIncludePath());
   } else {           
      ::Info("LoadCommodLibraries", "Load common libraries:    FAILED");
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
   if (iPWG4partcorr || iPWG4pi0) {   
      if (!LoadLibrary("PWG4PartCorrBase", mode, kTRUE) ||
          !LoadLibrary("PWG4PartCorrDep", mode, kTRUE)) return kFALSE;
      TFile::Cp(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG4/macros/ConfigAnalysisGammaHadronCorrelation.C"), "ConfigAnalysisGammaHadronCorrelation.C");
      TFile::Cp(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG4/macros/ConfigAnalysisPi0.C"), "ConfigAnalysisPi0.C");
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
   // PWG2 femtoscopy
   if (iPWG2femto) {
      if (!LoadLibrary("PWG2AOD", mode, kTRUE) ||
          !LoadLibrary("PWG2femtoscopy", mode, kTRUE) ||
          !LoadLibrary("PWG2femtoscopyUser", mode, kTRUE)) return kFALSE;
   }   
   // Vertexing HF
   if (iPWG3vertexing) {
      if (!LoadLibrary("PWG3base", mode, kTRUE) ||
          !LoadLibrary("PWG3vertexingHF", mode, kTRUE)) return kFALSE;
   }   
   ::Info("LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
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
      ::Error("LoadLibrary", "Empty module name");
      return kFALSE;
   }   
   // If a library is specified, just load it
   if (smodule.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      result = gSystem->Load(mod);
      if (result < 0) {
         ::Error("LoadLibrary", "Could not load library %s", module);
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
               ::Error("LoadLibrary", "Could not find module %s.par in current directory nor in $ALICE_ROOT", module);
               return kFALSE;
            }
         }   
         result = gProof->EnablePackage(module);
         break;
      default:
         return kFALSE;
   }         
   if (result < 0) {
      ::Error("LoadLibrary", "Could not load module %s", module);
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
            if (!strlen(dataset)) {
               // Local AOD
               chain = new TChain("aodTree");
               if (gSystem->AccessPathName("AliAOD.root")) 
                  ::Error("CreateChain", "File: AliAOD.root not in current dir");
               else chain->Add("AliAOD.root");
            } else {
               // Interactive AOD
               chain = CreateChainSingle(dataset, "aodTree");
            }
         } else {      
            if (!strlen(dataset)) {
               // Local ESD
               chain = new TChain("esdTree");
               if (gSystem->AccessPathName("data/AliESDs.root")) 
                  ::Error("CreateChain", "File: AliESDs.root not in current dir");
               else chain->Add("data/AliESDs.root");
            } else {
               // Interactive ESD
               chain = CreateChainSingle(dataset, "esdTree");
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
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
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
         ::Info("SetupPar", "Getting %s.par from $ALICE_ROOT", pararchivename);
         TFile::Cp(gSystem->ExpandPathName(Form("$ALICE_ROOT/%s.par", pararchivename)), 
                   Form("%s.par",pararchivename));
      } else {
         ::Error("SetupPar", "Cannot find %s.par", pararchivename);
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
//   plugin->SetPreferedSE("ALICE::NIHAM::FILE");
// Set versions of used packages
   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion(root_version);
   plugin->SetAliROOTVersion(aliroot_version);
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir(alien_datadir);
// Set data search pattern
   plugin->SetDataPattern("*ESDs.root");
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
   plugin->SetGridWorkingDir("analysisESD");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

   TString ana_sources = "";
   TString ana_add = "";
   if (usePAR && anaPars.Length()) {
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
//   plugin->SetOutputFiles("AliAOD.root");
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisTrainGrid.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(50);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(30000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("AnalysisTrain.jdl");
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("AnalysisTrain.sh");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   return plugin;
}
