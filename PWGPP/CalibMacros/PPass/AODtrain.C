// ### Settings that make sense when using the Alien plugin
//==============================================================================
Int_t       runOnData          = 1;       // Set to 1 if processing real data
Int_t       iCollision         = 0;       // 0=pp, 1=Pb-Pb
Int_t       run_flag           = 1100;    // year (2010 pp)
//==============================================================================
Bool_t      doCDBconnect        =1;
Bool_t      usePhysicsSelection = kTRUE; // use physics selection
Bool_t      useTender           = kFALSE; // use tender wagon
Bool_t      useCentrality       = kFALSE; // centrality
Bool_t      useV0tender         = kFALSE;  // use V0 correction in tender
Bool_t      useDBG              = kTRUE;  // activate debugging
Bool_t      useMC               = kFALSE;  // use MC info
Bool_t      useKFILTER          = kFALSE;  // use Kinematics filter
Bool_t      useTR               = kFALSE;  // use track references
Bool_t      useCORRFW           = kFALSE; // do not change
Bool_t      useAODTAGS          = kFALSE; // use AOD tags
Bool_t      useSysInfo          = kTRUE; // use sys info

// ### Analysis modules to be included. Some may not be yet fully implemented.
//==============================================================================
Int_t       iAODhandler        = 1;      // Analysis produces an AOD or dAOD's
Int_t       iESDfilter         = 1;      // ESD to AOD filter (barrel + muon tracks)
Int_t       iMUONcopyAOD       = 1;      // Task that copies only muon events in a separate AOD (PWG3)
Int_t       iJETAN             = 0;      // Jet analysis (PWG4)
Int_t       iJETANdelta        = 0;      // Jet delta AODs
Int_t       iPWGHFvertexing     = 1;      // Vertexing HF task (PWG3)
Int_t       iPWGDQJPSIfilter    = 0;      // JPSI filtering (PWG3)
Int_t       iPWGHFd2h           = 1;      // D0->2 hadrons (PWG3)
Int_t       iPWGPP               =1;      // high pt filter task
Int_t       iPWGLFForward       = 1;      // Forward mult task (PWGLF)
Int_t       iPWGGAgammaconv    = 1;      // Gamma conversion analysis (PWG4)
Bool_t   doPIDResponse  = 1;
Bool_t   doPIDqa        = 1; //new

// ### Configuration macros used for each module
//==============================================================================
 TString configPWGHFd2h = (iCollision==0)?"$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF.C"
                          :"$ALICE_PHYSICS/PWGHF/vertexingHF/ConfigVertexingHF_highmult.C";
                                                 
// Temporaries.
void AODmerge();
void AddAnalysisTasks();
Bool_t LoadCommonLibraries();
Bool_t LoadAnalysisLibraries();
Bool_t LoadLibrary(const char *);
TChain *CreateChain();
Int_t run_number = 0;

//______________________________________________________________________________
void AODtrain(Int_t merge=0, const char* cdbPath="cvmfs://", int collisionSystem=0)
{
// Main analysis train macro.
// collision system: 0-pp, 1-pbpb
  
  TString cdbPathTString=cdbPath;
  Bool_t initAlienConnection=kFALSE;
  if (cdbPathTString.BeginsWith("alien://")) initAlienConnection=kTRUE;

  if (initAlienConnection || merge) {
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

   // Make the analysis manager and connect event handlers
   AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Production train");
    mgr->SetCacheSize(0);
   if (useSysInfo) {
 //     mgr->SetNSysInfo(100);
      AliSysInfo::SetVerbose(kTRUE);
   }   
   // Create input handler (input container created automatically)
   // ESD input handler
   AliESDInputHandler *esdHandler = new AliESDInputHandler();
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

   AddAnalysisTasks(cdbPath,collisionSystem);
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
void AddAnalysisTasks(const char *cdb_location, int collisionSystem=0){                                                                                                                                          
  // Add all analysis task wagons to the train                                                                                               
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();                                                                     

  //
  // Tender and supplies. Needs to be called for every event.
  //
   AliAnalysisManager::SetCommonFileName("AODQA.root");
   if (useTender) {
      gROOT->LoadMacro("$ALICE_PHYSICS/TENDER/TenderSupplies/AddTaskTender.C");
      // IF V0 tender needed, put kTRUE below
      AliAnalysisTaskSE *tender = AddTaskTender(useV0tender);
//      tender->SetDebugLevel(2);
   }
   
   // Clean Geometry: Ruben
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/CalibMacros/commonMacros/CleanGeom.C++");
//  CleanGeom* clgmTask = new CleanGeom("cleanGeom");
//  mgr->AddTask(clgmTask);
//  AliAnalysisDataContainer *dummyInp = mgr->GetCommonInputContainer();
//  if (dummyInp) mgr->ConnectInput(clgmTask,0,dummyInp);  
 
   
   //
  // PIDResponse(JENS)
  //
  if (doPIDResponse) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); 
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse();
 //    PIDResponse->SetUserDataRecoPass(1);
//    PIDResponse->SelectCollisionCandidates(AliVEvent::kAny);
  }  
 
  //
  // PIDqa(JENS)
  //
  if (doPIDqa) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa *PIDQA = AddTaskPIDqa();
    PIDQA->SelectCollisionCandidates(AliVEvent::kAny);
  }  
  // CDB connection
  //
  if (doCDBconnect && !useTender) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect(cdb_location, run_number);
    if (!taskCDB) return;
    //    AliCDBManager *cdb = AliCDBManager::Instance(); // RS Never do this: task locks the CDB
    //    cdb->SetDefaultStorage(cdb_location);
//    taskCDB->SetRunNumber(run_number);
  }    
 
   if (usePhysicsSelection) {
   // Physics selection task
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      mgr->RegisterExtraFile("event_stat.root");
      AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(useMC);
//      AliOADBPhysicsSelection * oadbDefaultPbPb = CreateOADBphysicsSelection();      
//      physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadbDefaultPbPb,0,0);
      mgr->AddStatisticsTask(AliVEvent::kAny);
   }
   

//Jacek
   if (iPWGPP) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskFilteredTree.C");
      AddTaskFilteredTree("FilterEvents_Trees.root");
   }   
   
   // Centrality (only Pb-Pb)
   if (useCentrality) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
      //taskCentrality->SelectCollisionCandidates(AliVEvent::kAny);
   }
   
// --- PWGLF - Forward (cholm@nbi.dk) -----------------------------
   if (iPWGLFForward && usePhysicsSelection) { 
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");
     UShort_t pwglfForwardSys = 0; // iCollision+1; // pp:1, PbPb:2, pPb:3
     UShort_t pwglfSNN        = 0;            // GeV, 0==unknown
     Short_t  pwglfField      = 0;
     AddTaskForwardMult(useMC && useTR,        // Need track-refs 
			pwglfForwardSys,       // Collision system
			pwglfSNN, 
			pwglfField);
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskCentralMult.C");
        AddTaskCentralMult(useMC, pwglfForwardSys, pwglfSNN, pwglfField);
   }

   // PWGAgammaconv
   if (iPWGGAgammaconv) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_ConversionAODProduction.C");
      AliAnalysisTask *taskconv = AddTask_ConversionAODProduction(iCollision, kFALSE);
      mgr->RegisterExtraFile("AliAODGammaConversion.root");
   }

   if (iESDfilter) {
      //  ESD filter task configuration.
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
      if (iMUONcopyAOD) {
         printf("Registering delta AOD file\n");
         mgr->RegisterExtraFile("AliAOD.Muons.root");
          mgr->RegisterExtraFile("AliAOD.Dimuons.root");
      }
      AliAnalysisTaskESDfilter *taskesdfilter = 
                 AddTaskESDFilter(useKFILTER, 
                                  iMUONcopyAOD,         // write Muon AOD
                                  kFALSE,               // write dimuon AOD 
                                  kFALSE,               // usePhysicsSelection 
                                  kFALSE,               // centrality OBSOLETE
                                  kTRUE,                // enable TPS only tracks
                                  kFALSE,               // disable cascades
                                  kFALSE,               // disable kinks
                                  run_flag);            // run flag (YY00)
        mgr->RegisterExtraFile("AliAOD.Dimuons.root");
         AliEMCALGeometry::GetInstance("","");
   }   

// ********** PWG3 wagons ******************************************************           
   // PWGHF vertexing
   if (iPWGHFvertexing) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskVertexingHF.C");
      if (!iPWGHFd2h) TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), "file:ConfigVertexingHF.C");
      AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF(collisionSystem);
      if (!taskvertexingHF) ::Warning("AnalysisTrainNew", "AliAnalysisTaskSEVertexingHF cannot run for this train conditions - EXCLUDED");
      else mgr->RegisterExtraFile("AliAOD.VertexingHF.root");
      taskvertexingHF->SelectCollisionCandidates(0);
   }   
      
   // PWGDQ JPSI filtering (only pp)
   if (iPWGDQJPSIfilter && (iCollision==0)) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/dielectron/macros/AddTaskJPSIFilter.C");
      AliAnalysisTaskSE *taskJPSIfilter = AddTaskJPSIFilter();
      if (!taskJPSIfilter) ::Warning("AnalysisTrainNew", "AliAnalysisTaskDielectronFilter cannot run for this train conditions - EXCLUDED");
      else mgr->RegisterExtraFile("AliAOD.Dielectron.root");
      taskJPSIfilter->SelectCollisionCandidates(0);
   }   

   // PWGHF D2h
   if (iPWGHFd2h) {   
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/AddD2HTrain.C");
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
   UInt_t kHighPtFilterMask = 272;// from esd filter
   UInt_t iPhysicsSelectionFlag = 0;
   if (iJETAN) {
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskJets.C");
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
       gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskJetCluster.C");
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
         gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/AddTaskJetBackgroundSubtract.C");
         AliAnalysisTaskJetBackgroundSubtract *taskSubtract = 0;
         taskSubtract = AddTaskJetBackgroundSubtract(kJetSubtractBranches,1,"B0","B%d");
         taskSubtract->SetBackgroundBranch(kDefaultJetBackgroundBranch.Data());
         if(kDeltaAODJetName.Length()>0)taskSubtract->SetNonStdOutputFile(kDeltaAODJetName.Data());
       }
     } 
   }
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
  TString outputFiles = "EventStat_temp.root,AODQA.root,AliAOD.root,AliAOD.VertexingHF.root,AliAODGammaConversion.root,FilterEvents_Trees.root,AliAOD.Muons.root";
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
       continue;
    }
  }
  // all outputs merged, validate
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->InitAnalysis();
  mgr->SetGridHandler(new AliAnalysisAlien);
  mgr->StartAnalysis("gridterminate",0);
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  timer.Print();
}
