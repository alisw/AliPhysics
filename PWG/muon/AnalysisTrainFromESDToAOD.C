class AliAnalysisGrid; 
const char *dataset   = "";
TString anaLibs = "";
Int_t iESDfilter       = 1;
Int_t iAODTagCreation  = 1;
Int_t iAODAddMCBranch  = 1;

void AnalysisTrainFromESDToAOD(const char *analysisMode = "GRID", Bool_t usePLUGIN = kFALSE, Int_t nev=12345678)

//========================================================================
// The macro produces a standard AOD starting from ESD files. 
// (Simplified version of ANALYSIS/macros/AnalysisTrainNew.C)
//
// Two wagons are attached to the train: 
//      AddTaskESDFilter.C and AddTaskTagCreation.C 
//
// If the iESDfilter flag is activated, in AddTaskESDFilter.C 
// two tasks are executed :
// 1- with the first one (AliAnalysisTaskESDfilter), 
//    all the branches of the AOD are filled apart from the muons 
// 2- with the second task (AliAnalysisTaskESDMuonFilter) 
//    muons tracks are added to the tracks branch 
//
// In AddTaskESDFilter.C there is the possibility to apply cuts 
// on the tracks and muon tracks in order to reject them before 
// filling the AOD. 
//
// - if the flag iAODAddMCBranch is activated the MC branch 
// (containing Kinematics info) is added to the AOD 

// - if the iAODTagCreation flag is activated, in AddTaskTagCreation.C the 
//  AliAnalysisTaskTagCreator task is executed in order to create aod tags.  
//
//
// Options tested: (case sensitive)
//    GRID (with/without AliEn plugin)
//    LOCAL (you have to specify in TChain *CreateChain(...) 
//          the directory where your data are)
//========================================================================

{
    // Global configuration flags 
    //=====================================================================
    Bool_t debug         = kTRUE;
    Bool_t readTR        = kFALSE;      
    Bool_t useKFILTER    = kFALSE;  // add MC Branch 
    if(iAODAddMCBranch) useKFILTER=kTRUE;
    if(strcmp(analysisMode,"LOCAL")==0) usePLUGIN = kFALSE; 
       
    // Load common libraries (STEERBase, ESD, AOD, ANALYSIS. ANALYSISalice)
    //=====================================================================
    LoadCommonLibraries(analysisMode);
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
   
    // Load analysis specific libraries
    //=====================================================================
    if (iESDfilter) {
       if(!strcmp(analysisMode, "LOCAL")){
         gSystem->Load("libPWGHFbase");
         gSystem->Load("libPWGmuon");
       } 
       else if(!strcmp(analysisMode, "GRID"))LoadAnalysisLibraries(analysisMode);
     }
     
    // If Plugin is used, load macro with JDL parameters
    //=====================================================================
    if(usePLUGIN){
      gROOT->LoadMacro("CreateAlienHandler_FromESDToAOD.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandler_FromESDToAOD();  
      if (!alienHandler) return;     
    }

     // Create the chain. This is dependent on the analysis mode.
     //=====================================================================
     if(!usePLUGIN){
       if (!strcmp(analysisMode, "GRID")) TGrid::Connect("alien://");
       TChain* chain = CreateChain(analysisMode,""); 
     }
     
    // Create the train and set-up the handlers
    //=====================================================================
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis Train for standard AOD production");

    // GRID handler
    if(usePLUGIN) mgr->SetGridHandler(alienHandler); 
    
    // ESD input handler
       AliESDInputHandler *esdHandler = new AliESDInputHandler();
       mgr->SetInputEventHandler(esdHandler);      

    // Monte Carlo handler
    if (iAODAddMCBranch) {
       AliMCEventHandler* mcHandler = new AliMCEventHandler();
       mgr->SetMCtruthEventHandler(mcHandler);
       mcHandler->SetReadTR(readTR);
    } 
     
    // AOD output handler
    AliAODHandler* aodHandler   = new AliAODHandler();
    mgr->SetOutputEventHandler(aodHandler);
    aodHandler->SetOutputFileName("AliAODs.root");
       
    // Debugging if requested
    if (debug) mgr->SetDebugLevel(3);
       
    // Load the tasks configuration macros for all included wagons
    //=====================================================================    
    if (iESDfilter) {
       gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
       AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter(useKFILTER);
       
       if (iAODTagCreation) {
// use this line if AddTaskTagCreation.C is available in the grid aliroot version
       if(!strcmp(analysisMode, "LOCAL")){
   	  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muon/AddTaskTagCreation.C");
       } else  {
// uncomment  this line if AddTaskTagCreation.C is available in the grid aliroot version
//   	  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muon/AddTaskTagCreation.C");
// otherwise temporary: (and AddTaskTagCreation.C must be also added in the jdl)
         gROOT->LoadMacro("AddTaskTagCreation.C");
       }	 
	AliAnalysisTaskTagCreator *tagcreator = AddTaskTagCreation();
      }       
    }  

    // Run the analysis
    //=====================================================================    
    if (mgr->InitAnalysis()) {
       mgr->PrintStatus();
       if(usePLUGIN) mgr->StartAnalysis("GRID");
       else mgr->StartAnalysis("local", chain,nev);
    }  
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
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");
   
   // Load framework classes. Par option ignored here.
   switch (imode) {
      case 0:
      case 2:
            success &= LoadLibrary("libSTEERBase.so", mode);
            success &= LoadLibrary("libESD.so", mode);
            success &= LoadLibrary("libAOD.so", mode);
            success &= LoadLibrary("libANALYSIS.so", mode);
            success &= LoadLibrary("libANALYSISalice.so", mode);
            success &= LoadLibrary("libCORRFW.so", mode);
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
            success &= LoadLibrary("CORRFW", mode);
         } else { 
            ires = gProof->EnablePackage(AFversion);
            if (ires<0) success = kFALSE;
         }
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
            if (!strlen(dataset)) {
               // Local ESD
               chain = new TChain("esdTree");
               if (gSystem->AccessPathName("/n60raid3/alice/roberta/MCBranch/AliESDs.root")) 
                  ::Error("CreateChain", "File: AliESDs.root not in ./data dir");
               else chain->Add("/n60raid3/alice/roberta/MCBranch/AliESDs.root");
            } else {
               // Interactive ESD
               chain = CreateChainSingle(dataset, "esdTree");
            }   
         break;
      case 1:
         break;
      case 2:
            TString  treeName = "esdTree";
            chain = CreateChainSingle("wn.xml", treeName);
         break;      
      default:   
   }
   if (chain && chain->GetNtrees()) return chain;
   return NULL;
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
   if (strlen(gSystem->GetLibraries(module, "", kFALSE)) > 0)
      return kTRUE;    
   switch (imode) {
      case 0:
      case 2:
            result = gSystem->Load(Form("lib%s", module));
            if (rec) anaLibs += Form("lib%s.so ", module);
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
TChain* CreateChainSingle(const char* xmlfile, const char *treeName)
{
   printf("*******************************\n");
   printf("*** Getting the ESD Chain   ***\n");
   printf("*******************************\n");
   TGridCollection * myCollection  = gGrid->OpenCollection(xmlfile);

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
Bool_t LoadAnalysisLibraries(const char *mode)
{
// Load common analysis libraries.
   Bool_t success = kTRUE;
      if (!LoadLibrary("PWG3base", mode, kTRUE) ||
          !LoadLibrary("PWG3muon", mode, kTRUE)) return kFALSE;
   ::Info("LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
   return kTRUE;
}

