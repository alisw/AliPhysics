// ### Settings that make sense when using the Alien plugin
//==============================================================================
Int_t       runOnData          = 0;       // Set to 1 if processing real data
Int_t       iCollision         = 0;       // 0=pp, 1=Pb-Pb
//==============================================================================
Bool_t      doCDBconnect        = kFALSE;
Bool_t      usePhysicsSelection = kTRUE; // use physics selection
Bool_t      useTender           = kFALSE; // use tender wagon
Bool_t      useCentrality       = kTRUE; // centrality
Bool_t      useV0tender         = kFALSE;  // use V0 correction in tender
Bool_t      useDBG              = kTRUE;  // activate debugging
Bool_t      useMC               = kTRUE;  // use MC info
Bool_t      useKFILTER          = kTRUE;  // use Kinematics filter
Bool_t      useTR               = kTRUE;  // use track references
Bool_t      useCORRFW           = kFALSE; // do not change
Bool_t      useAODTAGS          = kFALSE; // use AOD tags
Bool_t      useSysInfo          = kFALSE; // use sys info

// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iESDMCLabelAddition= 0;
Int_t       iESDfilter         = 1;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iMUONcopyAOD       = 1;      // Task that copies only muon events in a separate AOD (PWG3)
Int_t       iJETAN             = 0;      // Jet analysis (PWG4)
Int_t       iJETANdelta        = 0;      // Jet delta AODs
Int_t       iPWGHFvertexing     = 0;      // Vertexing HF task (PWG3)
Int_t       iPWGDQJPSIfilter    = 0;      // JPSI filtering (PWG3)
Int_t       iPWGHFd2h           = 0;      // D0->2 hadrons (PWG3)
Bool_t doPIDResponse  = 1;
Bool_t doPIDqa        = 1; //new

// ### Configuration macros used for each module
//==============================================================================
 TString configPWGHFd2h = (iCollision==0)?"$ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF.C"
                          :"$ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF_Pb_AllCent.C";
                                                 
// Temporaries.
class AliOADBPhysicsSelection;                                                                                                                  
AliOADBPhysicsSelection *CreateOADBphysicsSelection();
void AODmerge();
void AddAnalysisTasks(Int_t);
Bool_t LoadCommonLibraries();
Bool_t LoadAnalysisLibraries();
Bool_t LoadLibrary(const char *);
TChain *CreateChain();
const char *cdbPath = "raw://";
Int_t run_number = 0;

//______________________________________________________________________________
void AODtrain(Int_t merge=0)
{
// Main analysis train macro.
   AliLog::SetGlobalDebugLevel(3);
  if (merge) {
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) {
      ::Error("QAtrain", "No grid connection");
      return;
    }
  }
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
   printf("==================================================================\n");
   printf("===========    RUNNING FILTERING TRAIN   ==========\n");
   printf("==================================================================\n");
   printf("=  Configuring analysis train for:                               =\n");
   if (usePhysicsSelection)   printf("=  Physics selection                                                =\n");
   if (useTender)    printf("=  TENDER                                                        =\n");
   if (iESDfilter)   printf("=  ESD filter                                                    =\n");
   if (iMUONcopyAOD) printf("=  MUON copy AOD                                                 =\n");
   if (iJETAN)       printf("=  Jet analysis                                                  =\n");
   if (iJETANdelta)  printf("=     Jet delta AODs                                             =\n");
   if (iPWGHFvertexing) printf("=  PWGHF vertexing                                                =\n");
   if (iPWGDQJPSIfilter) printf("=  PWGDQ j/psi filter                                             =\n");
   if (iPWGHFd2h) printf("=  PWGHF D0->2 hadrons QA                                     =\n");

   // Load common libraries and set include path
   if (!LoadCommonLibraries()) {
      ::Error("AnalysisTrain", "Could not load common libraries");
      return;
   }
    
   // Make the analysis manager and connect event handlers
   AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");
   if (useSysInfo) mgr->SetNSysInfo(100);
   // Load analysis specific libraries
   if (!LoadAnalysisLibraries()) {
      ::Error("AnalysisTrain", "Could not load analysis libraries");
      return;
   }   

   // Create input handler (input container created automatically)
   // ESD input handler
   AliESDInputHandler *esdHandler = new AliESDInputHandler();
   esdHandler->SetNeedField();
   mgr->SetInputEventHandler(esdHandler);       
   // Monte Carlo handler
   if (useMC) {
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
   }
   // Debugging if needed
   if (useDBG) mgr->SetDebugLevel(3);

   AddAnalysisTasks(merge, cdbPath);
   if (merge) {
      AODmerge();
      mgr->InitAnalysis();
      mgr->SetGridHandler(new AliAnalysisAlien);
      mgr->StartAnalysis("gridterminate",0);
      return;
   }   
   // Run the analysis                                                                                                                     
   //
   TChain *chain = CreateChain();
   if (!chain) return;
                                                                                                                                                   
   TStopwatch timer;
   timer.Start();
   mgr->SetSkipTerminate(kTRUE);
   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      mgr->StartAnalysis("local", chain);
   }
   timer.Print();
}                                                                                                                                          
                                                                                                                                            
//______________________________________________________________________________                                                           
void AddAnalysisTasks(Int_t merge, const char *cdb_location){                                                                                                                                          
  // Add all analysis task wagons to the train                                                                                               
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();                                                                     

  //
  // Tender and supplies. Needs to be called for every event.
  //
   AliAnalysisManager::SetCommonFileName("AODQA.root");
   if (useTender) {
      gROOT->LoadMacro("$ALICE_ROOT/TENDER/TenderSupplies/AddTaskTender.C");
      // IF V0 tender needed, put kTRUE below
      AliAnalysisTaskSE *tender = AddTaskTender(useV0tender);
//      tender->SetDebugLevel(2);
   }

  // CDB connection
  //
  if (doCDBconnect && !useTender) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect(cdb_location, run_number);
    if (!taskCDB) return;
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(cdb_location);
//    taskCDB->SetRunNumber(run_number);
  }    

   if (usePhysicsSelection) {
   // Physics selection task
      gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
      mgr->RegisterExtraFile("event_stat.root");
      AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kTRUE);
//      AliOADBPhysicsSelection * oadbDefaultPbPb = CreateOADBphysicsSelection();      
//      physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadbDefaultPbPb,0,0);
//      if (!merge) mgr->AddStatisticsTask(AliVEvent::kAny);
   }
   // Centrality (only Pb-Pb)
   if (iCollision && useCentrality) {
      gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
      taskCentrality->SelectCollisionCandidates(AliVEvent::kAny);
   }

	if(iESDMCLabelAddition) {
    		gROOT->LoadMacro("$ALICE_ROOT/PWG/muondep/AddTaskESDMCLabelAddition.C");
    		AliAnalysisTaskESDMCLabelAddition *esdmclabel = AddTaskESDMCLabelAddition(kFALSE);
	}
	
   if (iESDfilter) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");	
      if (iMUONcopyAOD) {
         printf("Registering delta AOD file\n");
         mgr->RegisterExtraFile("AliAOD.Muons.root");
         mgr->RegisterExtraFile("AliAOD.Dimuons.root");
	 AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(useKFILTER, kTRUE, kFALSE, kFALSE /*usePhysicsSelection*/,kFALSE,kTRUE,kTRUE,kTRUE,1100,1);
      } else {
    AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(useKFILTER, kFALSE, kFALSE, kFALSE /*usePhysicsSelection*/,kFALSE,kTRUE,kTRUE,kTRUE,1100,1);
      }   
   }   

// ********** PWG3 wagons ******************************************************           
   // PWGHF vertexing
   if (iPWGHFvertexing) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
      if (!iPWGHFd2h) TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), "file:ConfigVertexingHF.C");
      AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
      if (!taskvertexingHF) ::Warning("AnalysisTrainNew", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
      else mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
      taskvertexingHF->SelectCollisionCandidates(0);
   }   
      
   // PWGDQ JPSI filtering (only pp)
   if (iPWGDQJPSIfilter && (iCollision==0)) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macros/AddTaskJPSIFilter.C");
      AliAnalysisTaskSE *taskJPSIfilter = AddTaskJPSIFilter();
      if (!taskJPSIfilter) ::Warning("AnalysisTrainNew", "AliAnalysisTaskDielectronFilter cannot run for this train conditions - EXCLUDED");
      else mgr->RegisterExtraFile("AliAOD.Dielectron.root");
      taskJPSIfilter->SelectCollisionCandidates(0);
   }   

   // PWGHF D2h
   if (iPWGHFd2h) {   
     gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/AddD2HTrain.C");
     TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), "file:ConfigVertexingHF.C");
     AddD2HTrain(kFALSE, 1,0,0,0,0,0,0,0,0,0,0);                                 
   }
   
   // ********** PWG4 wagons ******************************************************
   // Jet analysis

   // Configurations flags, move up?
   TString kDeltaAODJetName = "AliAOD.Jets.root"; //
   Bool_t  kIsPbPb = (iCollision==0)?false:true; // can be more intlligent checking the name of the data set
   TString kDefaultJetBackgroundBranch = "";
   TString kJetSubtractBranches = "";
   UInt_t kHighPtFilterMask = 128;// from esd filter
   UInt_t iPhysicsSelectionFlag = AliVEvent::kMB;
   if (iJETAN) {
     gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJets.C");
     // Default jet reconstructor running on ESD's
     AliAnalysisTaskJets *taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,1.,0); // no background subtraction     
     if (!taskjets) ::Fatal("AnalysisTrainNew", "AliAnalysisTaskJets cannot run for this train conditions - EXCLUDED");
     if(kDeltaAODJetName.Length()>0) taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
     if (iJETANdelta) {
        //            AddTaskJetsDelta("AliAOD.Jets.root"); // need to modify this accordingly in the add task jets
        mgr->RegisterExtraFile(kDeltaAODJetName.Data());
        TString cTmp("");
        if(kIsPbPb){
          // UA1 intrinsic background subtraction
          taskjets = AddTaskJets("AOD","UA1",0.4,kHighPtFilterMask,1.,2); // background subtraction
          if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
       }
       // SICONE 
       taskjets = AddTaskJets("AOD","SISCONE",0.4,kHighPtFilterMask,0.15,0); //no background subtration to be done later....                                                                                  
       if(kDeltaAODJetName.Length()>0)taskjets->SetNonStdOutputFile(kDeltaAODJetName.Data());
       cTmp = taskjets->GetNonStdBranch();
       if(cTmp.Length()>0)kJetSubtractBranches += Form("%s ",cTmp.Data());
	 
       // Add the clusters..
       gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetCluster.C");
       AliAnalysisTaskJetCluster *taskCl = 0;
       Float_t fCenUp = 0;
       Float_t fCenLo = 0;
       Float_t fTrackEtaWindow = 0.9;
       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"KT",0.4,0,1, kDeltaAODJetName.Data(),0.15,fTrackEtaWindow,0); // this one is for the background and random jets, random cones with no skip                                                                                 
       taskCl->SetBackgroundCalc(kTRUE);
       taskCl->SetNRandomCones(10);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       taskCl->SetGhostEtamax(fTrackEtaWindow);
       kDefaultJetBackgroundBranch = Form("%s_%s",AliAODJetEventBackground::StdBranchName(),taskCl->GetJetOutputBranch());

       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.4,2,1,kDeltaAODJetName.Data(),0.15);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());
       taskCl->SetNRandomCones(10);
       kJetSubtractBranches += Form("%s ",taskCl->GetJetOutputBranch());

       taskCl = AddTaskJetCluster("AOD","",kHighPtFilterMask,iPhysicsSelectionFlag,"ANTIKT",0.2,0,1,kDeltaAODJetName.Data(),0.15);
       taskCl->SetCentralityCut(fCenLo,fCenUp);
       if(kIsPbPb)taskCl->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());
       kJetSubtractBranches += Form("%s ",taskCl->GetJetOutputBranch());
	 
       // DO THE BACKGROUND SUBTRACTION
       if(kIsPbPb&&kJetSubtractBranches.Length()){
         gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskJetBackgroundSubtract.C");
         AliAnalysisTaskJetBackgroundSubtract *taskSubtract = 0;
         taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranches,1,"B0","B%d");
         taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());
         if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());
       }
     } 
   }
}

//______________________________________________________________________________
Bool_t LoadCommonLibraries()
{

  // Load common analysis libraries.
   if (!gSystem->Getenv("ALICE_ROOT")) {
      ::Error("AnalysisTrainNew.C::LoadCommonLibraries", "Analysis train requires that analysis libraries are compiled with a local AliRoot"); 
      return kFALSE;
   }   
   Bool_t success = kTRUE;
   // Load framework classes. Par option ignored here.
   success &= LoadLibrary("libSTEERBase.so");
   success &= LoadLibrary("libESD.so");
   success &= LoadLibrary("libAOD.so");
   success &= LoadLibrary("libANALYSIS.so");
   success &= LoadLibrary("libOADB.so");
   success &= LoadLibrary("libANALYSISalice.so");
   success &= LoadLibrary("libESDfilter.so");
   success &= LoadLibrary("libCORRFW.so");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
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
Bool_t LoadAnalysisLibraries()
{

  // Load common analysis libraries.
  if (useTender || doCDBconnect) {
    if (!LoadLibrary("Tender") ||
	!LoadLibrary("TenderSupplies")) return kFALSE;
  }       
  // CDBconnect
  if (doCDBconnect && !useTender) {
    if (!LoadLibrary("PWGPP")) return kFALSE;
  }
  
  if (iESDfilter || iPWGMuonTrain) {
    if (!LoadLibrary("PWGmuon")) return kFALSE;
  }   
  if (iESDMCLabelAddition) {
    if (!LoadLibrary("PWGmuondep")) return kFALSE;
  }
  // JETAN
  if (iJETAN) {
    if (!LoadLibrary("JETAN")) return kFALSE;
  }
  if (iJETANdelta) {
    if (!LoadLibrary("JETAN") ||
	!LoadLibrary("CGAL") ||
	!LoadLibrary("fastjet") ||
	!LoadLibrary("siscone") ||
	!LoadLibrary("SISConePlugin") ||
	!LoadLibrary("FASTJETAN")) return kFALSE;
  }     
  // PWG3 Vertexing HF
  if (iPWGHFvertexing || iPWGHFd2h) {
    if (!LoadLibrary("PWGflowBase") ||
	!LoadLibrary("PWGflowTasks") ||
	!LoadLibrary("PWGHFvertexingHF")) return kFALSE;
  }    
  //    if (iPWGHFvertexing || iPWG3d2h) {
  //     if (!LoadLibrary("PWG3base") ||
  //         !LoadLibrary("PWGHFvertexingHF")) return kFALSE;
  //  }   
  // PWG3 dielectron
  if (iPWGDQJPSIfilter) {
    if (!LoadLibrary("PWGDQdielectron")) return kFALSE;
  }   
  
  ::Info("AnalysisTrainNew.C::LoadAnalysisLibraries", "Load other libraries:   SUCCESS");
  return kTRUE;
}

//______________________________________________________________________________
Bool_t LoadLibrary(const char *module)
{
// Load a module library in a given mode. Reports success.
   Int_t result;
   TString mod(module);
   if (!mod.Length()) {
      ::Error("AnalysisTrainNew.C::LoadLibrary", "Empty module name");
      return kFALSE;
   }   
   // If a library is specified, just load it
   if (mod.EndsWith(".so")) {
      mod.Remove(mod.Index(".so"));
      result = gSystem->Load(mod);
      if (result < 0) {
         ::Error("AnalysisTrainNew.C::LoadLibrary", "Could not load library %s", module);
         return kFALSE;
      }
      return kTRUE;
   } 
   // Check if the library is already loaded
   if (strlen(gSystem->GetLibraries(module, "", kFALSE)) > 0) return kTRUE;    
   result = gSystem->Load(Form("lib%s", module));
   if (result < 0) {
      ::Error("AnalysisTrainNew.C::LoadLibrary", "Could not load module %s", module);
      return kFALSE;
   }
   return kTRUE;
}           


//______________________________________________________________________________
TChain *CreateChain()
{
// Create the input chain
   chain = new TChain("esdTree");
   if (gSystem->AccessPathName("AliESDs.root")) 
      ::Error("AnalysisTrainNew.C::CreateChain", "File: AliESDs.root not in ./data dir");
   else 
       chain->Add("AliESDs.root");
   if (chain->GetNtrees()) return chain;
   return NULL;
}   

//______________________________________________________________________________
void AODmerge()
{
// Merging method. No staging and no terminate phase.
  TStopwatch timer;
  timer.Start();
  TString outputDir = "wn.xml";
  TString outputFiles = "EventStat_temp.root,AliAOD.root,AliAOD.Muons.root";
  TString mergeExcludes = "";
  TObjArray *list = outputFiles.Tokenize(",");
  TIter *iter = new TIter(list);
  TObjString *str;
  TString outputFile;
  Bool_t merged = kTRUE;
  while((str=(TObjString*)iter->Next())) {
    outputFile = str->GetString();
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
       printf("Output file <%s> found. Not merging again.",outputFile.Data());
       continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, 0);
    if (!merged) {
       printf("ERROR: Cannot merge %s\n", outputFile.Data());
       return;
    }
  }
  // all outputs merged, validate
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  timer.Print();
}
