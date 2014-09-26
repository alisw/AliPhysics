/**
 * @file   AOD.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Sep 24 15:02:00 2014
 * 
 * @brief  Master script for AOD production.
 * 
 * @note Do not modify this script. 
 *
 * This script reads in two other scripts 
 *
 * - GRP.C to load the global run parameters for the selected run,
 *   such as collision system, energy, etc.
 * 
 * - AODConfig.C which defines a number of functions that return
 *   either true or false.  The tasks added depends on these functions
 */
/** Path to CDB */
const char *cdbPath = "raw://";

/**
 * Interface (pure virtual) that all configuration classes must
 * implement.
 */
struct VirtualAODCfg
{
  /** 
   * @{ 
   * @name Plug-in settings 
   * Settings that make sense when using the Alien plugin
   */
  /** @return Connect to CDB */
  virtual Bool_t UseCDBconnect() const = 0;
  /** @return use physics selection */
  virtual Bool_t UsePhysicsSelection() const = 0;
  /** @return use tender wagon */
  virtual Bool_t UseTender() const = 0;
  /** @return centrality */
  virtual Bool_t UseCentrality() const = 0;
  /** @return use V0 correction in tender */
  virtual Bool_t UseV0tender() const = 0;
  /** @return activate debugging */
  virtual Bool_t UseDBG() const = 0;
  /** @return use MC info */
  virtual Bool_t UseMC() const = 0;
  /** @return use Kinematics filter */
  virtual Bool_t UseKFILTER() const = 0;
  /** @return use track references */
  virtual Bool_t UseTR() const = 0;
  /** @return do not change */
  virtual Bool_t UseCORRFW() const = 0;
  /** @return use AOD tags */
  virtual Bool_t UseAODTAGS() const = 0;
  /** @return use sys info */
  virtual Bool_t UseSysInfo() const = 0;
  /* @} */
  
  /** 
   * @{ 
   * @name Modules 
   *  Analysis modules to be included. Some may not be yet fully implemented.
   */
  /** @return Analysis produces an AOD or dAOD's */
  virtual Bool_t UseAODhandler() const = 0;
  /** @return ESD to AOD filter (barrel + muon tracks) */
  virtual Bool_t UseESDfilter() const = 0;
  /** @return Use Muon train  */
  virtual Bool_t UsePWGMuonTrain() const = 0;
  /** @return Task that copies only muon events */
  virtual Bool_t UseMUONcopyAOD() const = 0;
  /** @return Jet analysis (PWG4) */
  virtual Bool_t UseJETAN() const = 0;
  /** @return Jet delta AODs */
  virtual Bool_t UseJETANdelta() const = 0;
  /** @return Vertexing HF task (PWG3) */
  virtual Bool_t UsePWGHFvertexing() const = 0;
  /** @return JPSI filtering (PWG3) */
  virtual Bool_t UsePWGDQJPSIfilter() const = 0;
  /** @return D0->2 hadrons (PWG3) */
  virtual Bool_t UsePWGHFd2h() const = 0;
  /** @return PID response */
  virtual Bool_t UsePIDResponse() const = 0;
  /** @return Forward mult task (PWGLF) */
  virtual Bool_t UsePWGLFForward() const = 0;
  /* @} */
};

VirtualAODCfg* aodCfg = 0;

//====================================================================
/** 
 * Load a library/module 
 * 
 * @param module Library/module name 
 * 
 * @return true on success
 */
Bool_t LoadLibrary(const char *module)
{
  // Load a module library in a given mode. Reports success.
  Int_t result = 0;
  TString mod(module);
  ::Info("LoadLibrary", "Loading %s", module);
  gROOT->IncreaseDirLevel();

  if (mod.IsNull()) {
    ::Error("AnalysisTrainNew.C::LoadLibrary", "Empty module name");
    gROOT->DecreaseDirLevel();
    return kFALSE;
  }

  // If a library is specified, just load it
  if (mod.EndsWith(".so")) {
    mod.Remove(mod.Index(".so"));
    ::Info("LoadLibrary", "Loading .so: %s", mod.Data()); 
    result = gSystem->Load(mod);
    if (result < 0) {
      ::Error("AnalysisTrainNew.C::LoadLibrary", 
	      "Could not load library %s", module);
    }
    gROOT->DecreaseDirLevel();      
    return (result >= 0);
  }
  // Check if the library is already loaded
  if (strlen(gSystem->GetLibraries(Form("%s.so", module), "", kFALSE)) > 0) {
    ::Info("LoadLibrary", "Module %s.so already loaded", module);
    gROOT->DecreaseDirLevel();      
    return kTRUE;
  }

  ::Info("LoadLibrary", "Trying to load lib%s.so", module);
  result = gSystem->Load(Form("lib%s.so", module));
  if (result < 0)
    ::Error("AnalysisTrainNew.C::LoadLibrary", 
	    "Could not load module %s", module);
  ::Info("LoadLibrary", "Module %s, successfully loaded", module);
  gROOT->DecreaseDirLevel();      
  return (result >= 0);
}

//====================================================================
/** 
 * Load common libraries 
 * 
 * @return true on sucess 
 */
Bool_t LoadCommonLibraries()
{
  // Load common analysis libraries.
  if (!gSystem->Getenv("ALICE_ROOT")) {
    ::Error("AnalysisTrainNew.C::LoadCommonLibraries", 
	    "Analysis train requires that analysis libraries are "
	    "compiled with a local AliRoot");
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
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Load common libraries:    SUCCESS");
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Include path for Aclic compilation:\n%s",
	   gSystem->GetIncludePath());
  } else {
    ::Info("AnalysisTrainNew.C::LoadCommodLibraries", 
	   "Load common libraries:    FAILED");
  }
  return success;
}

//====================================================================
/** 
 * Load libraries needed by the train 
 * 
 * @param useTender 
 * @param doCDBconnect 
 * @param iESDfilter 
 * @param iPWGMuonTrain 
 * @param iJETAN 
 * @param iJETANdelta 
 * @param iPWGHFvertexing 
 * @param iPWGHFd2h 
 * @param iPWGDQJPSIfilter 
 * 
 * @return true on success 
 */
Bool_t LoadAnalysisLibraries()
{
  // Load common analysis libraries.
  if (aodCfg->UseTender() || aodCfg->UseCDBconnect()) {
    if (!LoadLibrary("TENDER") ||!LoadLibrary("TENDERSupplies")) return kFALSE;
  }
  // CDBconnect
  if ((aodCfg->UseCDBconnect() && !aodCfg->UseTender()) 
      && !LoadLibrary("PWGPP")) return false;
  if ((aodCfg->UseESDfilter() || aodCfg->UsePWGMuonTrain()) 
      && !LoadLibrary("PWGmuon")) return kFALSE;
  // JETAN
  if ((aodCfg->UseJETAN() || aodCfg->UseJETANdelta()) 
      && !LoadLibrary("JETAN")) return kFALSE;
  if (aodCfg->UseJETANdelta()) { // CINT doesn't like long '||' chains
    if (!LoadLibrary("CGAL"))           return false;
    if (!LoadLibrary("fastjet"))        return false;
    if (!LoadLibrary("siscone")) 	return false;
    if (!LoadLibrary("SISConePlugin"))  return false;
    if (!LoadLibrary("FASTJETAN"))      return false;
  }

  // PWG3 Vertexing HF
  if ((aodCfg->UsePWGHFvertexing() || aodCfg->UsePWGHFd2h())) { 
    // CINT doesn't like long '||' chains
    if (!LoadLibrary("PWGflowBase"))      return false;
    if (!LoadLibrary("PWGflowTasks"))     return false;
    // if (!LoadLibrary("PWGTRD"))           return false;
    if (!LoadLibrary("PWGHFvertexingHF")) return false;
  }

  // PWG3 dielectron
  if (aodCfg->UsePWGDQJPSIfilter() && 
      !LoadLibrary("PWGDQdielectron")) return kFALSE;
  
  ::Info("AnalysisTrainNew.C::LoadAnalysisLibraries", 
	 "Load other libraries:   SUCCESS");
  return kTRUE;
}

//====================================================================
/** 
 * Add tasks to the train 
 * 
 * @param cdb_location 
 */
void AddAnalysisTasks(const char *cdb_location)
{
  // === Add all analysis task wagons to the train ===================
  // 
  // --- Some constants ----------------------------------------------
  TString ali   = "$(ALICE_ROOT)";
  TString ana   = ali + "/ANALYSIS";
  TString pwghf = ali + "/PWGHF";
  TString pwglf = ali + "/PWGLF";
  TString pwgje = ali + "/PWGJE";
  TString pwgdq = ali + "/PWGDQ";
  TString pwgpp = ali + "/PWGPP";

  // --- Get the analysis manager ------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisManager::SetCommonFileName("AODQA.root");

  // --- PIDResponse(JENS) -------------------------------------------
  if (aodCfg->UsePIDResponse()) {
    gROOT->LoadMacro(ana+"/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse(kTRUE);
  }

  // --- CDB connection ----------------------------------------------
  if (aodCfg->UseCDBconnect() && !aodCfg->UseTender()) {
    gROOT->LoadMacro(pwgpp+"/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect(cdb_location, grp->run);
    if (!taskCDB) return;

    AliCDBManager *cdb = AliCDBManager::Instance();
    // cdb->SetDefaultStorage(cdb_location);
    cdb->SetDefaultStorageFromRun(grp->run);
  }
  if (aodCfg->UseTender()) {
    gROOT->LoadMacro(ana+"/TenderSupplies/AddTaskTender.C");
    AliAnalysisTaskSE *tender = AddTaskTender(aodCfg->UseV0tender());
  }

  // --- Physics selection -------------------------------------------
  if (aodCfg->UsePhysicsSelection()) {
    // Physics selection task
    gROOT->LoadMacro(ana+"/macros/AddTaskPhysicsSelection.C");
    mgr->RegisterExtraFile("event_stat.root");
    AliPhysicsSelectionTask *physSelTask = 
      AddTaskPhysicsSelection(aodCfg->UseMC());
  }

  // --- Centrality (only Pb-Pb) -------------------------------------
  if (aodCfg->UseCentrality()) {
    gROOT->LoadMacro(ana+"/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    taskCentrality->SetMCInput();
  }

  // --- PWGLF - Forward (cholm@nbi.dk) -----------------------------
  if (aodCfg->UsePWGLFForward() && aodCfg->UsePhysicsSelection()) {
    gROOT->LoadMacro(pwglf+"/FORWARD/analysis2/AddTaskForwardMult.C");
    // Arguments are 
    //   mc         Assume MC input
    //   runNo      Run number to do local initialization - not used
    //   sys        Collision system (1:pp,2:PbPb,3:pPb/Pbp) - not used
    //   sNN        Collision energy in GeV - not used
    //   field      L3 magnetic field strength  - not used
    //   config     Configuration script 
    //   corrdir    Possible directory containing custom OADB corrections
    // HACK load custom corrections 
    Info("", "Adding forward AOD task with mc=%d",aodCfg->UseMC() && aodCfg->UseTR());
    AddTaskForwardMult(aodCfg->UseMC() && aodCfg->UseTR(),0,0,0,0,
		       "ForwardAODConfig.C",".");
    gROOT->LoadMacro(pwglf+"/FORWARD/analysis2/AddTaskCentralMult.C");
    AddTaskCentralMult(aodCfg->UseMC() && aodCfg->UseTR(),0,0,0,0,
		       "CentralAODConfig.C");
    const char* hack2="AliForwardCorrectionManager::Instance().Print(\"R\")";
    gROOT->ProcessLine(hack2);
  }


  // --- ESD filter --------------------------------------------------
  if (aodCfg->UseESDfilter()) {
    //  ESD filter task configuration.
    gROOT->LoadMacro(ana+"/ESDfilter/macros/AddTaskESDFilter.C");
    if (aodCfg->UseMUONcopyAOD()) {
      printf("Registering delta AOD file\n");
      mgr->RegisterExtraFile("AliAOD.Muons.root");
      mgr->RegisterExtraFile("AliAOD.Dimuons.root");
    }
    UInt_t runFlag = (grp->Year()%100)*100;
    AliAnalysisTaskESDfilter *taskesdfilter =
      AddTaskESDFilter(aodCfg->UseKFILTER(),
		       aodCfg->UseMUONcopyAOD(),         // write Muon AOD
		       kFALSE,               // write dimuon AOD
		       kFALSE,               // usePhysicsSelection
		       kFALSE,               // centrality OBSOLETE
		       kTRUE,                // enable TPS only tracks
		       kFALSE,               // disable cascades
		       kFALSE,               // disable kinks
		       runFlag);             // run flag (YY00)
  }

  TString configPWGHFd2h = 
    (grp->IsPP() ?
     pwghf+"/vertexingHF/ConfigVertexingHF.C" :
     pwghf+"/vertexingHF/ConfigVertexingHF_Pb_AllCent.C");

  // --- PWG3 vertexing ----------------------------------------------
  if (aodCfg->UsePWGHFvertexing()) {
    gROOT->LoadMacro(pwghf+"/vertexingHF/macros/AddTaskVertexingHF.C");

    if (!aodCfg->UsePWGHFd2h()) 
      TFile::Cp(gSystem->ExpandPathName(configPWG3d2h.Data()), 
		"file:ConfigVertexingHF.C");

    AliAnalysisTaskSEVertexingHF *taskvertexingHF = AddTaskVertexingHF();
    if (!taskvertexingHF) 
      ::Warning("AnalysisTrainNew", 
		"AliAnalysisTaskSEVertexingHF cannot run for this train "
		"conditions - EXCLUDED");
    else 
      mgr->RegisterExtraFile("AliAOD.VertexingHF.root");

    taskvertexingHF->SelectCollisionCandidates(0);
  }

  // ---- PWG3 JPSI filtering (only pp) ------------------------------
  if (aodCfg->UsePWGDQJPSIfilter() && (iCollision==0)) {
    gROOT->LoadMacro(pwgdq+"/dielectron/macros/AddTaskJPSIFilter.C");
    AliAnalysisTaskSE *taskJPSIfilter = AddTaskJPSIFilter();

    if (!taskJPSIfilter) 
      ::Warning("AnalysisTrainNew", 
		"AliAnalysisTaskDielectronFilter cannot run for this train "
		"conditions - EXCLUDED");
    else 
      mgr->RegisterExtraFile("AliAOD.Dielectron.root");

    taskJPSIfilter->SelectCollisionCandidates(0);
  }

  // --- PWG3 D2h ----------------------------------------------------
  if (aodCfg->UsePWGHFd2h()) {
    gROOT->LoadMacro(pwghf+"/vertexingHF/AddD2HTrain.C");
    TFile::Cp(gSystem->ExpandPathName(configPWGHFd2h.Data()), 
	      "file:ConfigVertexingHF.C");
    AddD2HTrain(kFALSE, 1,0,0,0,0,0,0,0,0,0,0);
  }

  // --- Jet analysis ------------------------------------------------

  // Configurations flags, move up?
  if (aodCfg->UseJETAN()) {
#if 1
    Warning("", "JET analysis disabled - major restructuring ofg JETAN");
#else
    TString jetAOD             = "AliAOD.Jets.root";
    UInt_t  highPtMask         = 768;// from esd filter
    TString subtractBranches   = "";
    UInt_t  psFlag             = 0;

    Info("", "Loading macro %s/macros/AddTaskJets.C", pwgje.Data());
    gROOT->LoadMacro(pwgje+"/macros/AddTaskJets.C");
    // Default jet reconstructor running on ESD's
    // no background subtraction
    AliAnalysisTask* task = AddTaskJets("AOD","UA1",0.4F,highPtMask,1.F,0); 
    if (!task) 
      ::Fatal("AnalysisTrainNew", 
	      "AliAnalysisTaskJets cannot run for this train "
	      "conditions - EXCLUDED");
    
    AliAnalysisTaskJets* taskjets = static_cast<AliAnalysisTaskJets*>(task);
    if(!jetAOD.IsNull()) taskjets->SetNonStdOutputFile(jetAOD);

    if (aodCfg->UseJETANdelta()) {
      // need to modify this accordingly in the add task jets
      // AddTaskJetsDelta("AliAOD.Jets.root"); 
      mgr->RegisterExtraFile(jetAOD);
      TString cTmp("");

      if(grp->IsAA()){
	// UA1 intrinsic background subtraction
	// background subtraction
	taskjets = AddTaskJets("AOD","UA1",0.4,highPtMask,1.,2); 
	if(!jetAOD.IsNull()) taskjets->SetNonStdOutputFile(jetAOD);
      }
      // SICONE
      //no background subtration to be done later....
      taskjets = AddTaskJets("AOD","SISCONE",0.4,highPtMask,0.15,0); 
      if(!jetAOD.IsNull()) taskjets->SetNonStdOutputFile(jetAOD.Data());
      cTmp = taskjets->GetNonStdBranch();
      if(!cTmp.IsNull()) subtractBranches += Form("%s ",cTmp.Data());

      // Add the clusters..
      gROOT->LoadMacro(pwgje+"/macros/AddTaskJetCluster.C");
      AliAnalysisTaskJetCluster *taskCl = 0;
      Float_t fCenUp = 0;
      Float_t fCenLo = 0;
      Float_t fTrackEtaWindow = 0.9;
      // this one is for the background and random jets, random cones
      // with no skip
      taskCl = AddTaskJetCluster("AOD","",highPtMask,
				 psFlag,"KT",0.4,0,1, 
				 jetAOD,0.15,
				 fTrackEtaWindow,0); 
      taskCl->SetBackgroundCalc(kTRUE);
      taskCl->SetNRandomCones(10);
      taskCl->SetCentralityCut(fCenLo,fCenUp);
      taskCl->SetGhostEtamax(fTrackEtaWindow);
      TString bkgBranch = Form("%s_%s",
			       AliAODJetEventBackground::StdBranchName(),
			       taskCl->GetJetOutputBranch());

      taskCl = AddTaskJetCluster("AOD", "", highPtMask, psFlag, "ANTIKT",
				 0.4, 2, 1, jetAOD, 0.15);
      taskCl->SetCentralityCut(fCenLo,fCenUp);
      if(grp->IsAA()) taskCl->SetBackgroundBranch(bkgBranch.Data());

      taskCl->SetNRandomCones(10);
      subtractBranches += Form("%s ",taskCl->GetJetOutputBranch());

      taskCl = AddTaskJetCluster("AOD", "", highPtMask, psFlag, "ANTIKT",
				 0.2, 0 , 1, jetAOD, 0.15);
      taskCl->SetCentralityCut(fCenLo,fCenUp);
      if(grp->IsAA())taskCl->SetBackgroundBranch(bkgBranch);
      
      subtractBranches += Form("%s ",taskCl->GetJetOutputBranch());

      // DO THE BACKGROUND SUBTRACTION
      if(grp->IsAA() && !subtractBranches.IsNull()) {
	gROOT->LoadMacro(pwgje+"/macros/AddTaskJetBackgroundSubtract.C");
	AliAnalysisTaskJetBackgroundSubtract *taskSubtract = 0;
	taskSubtract = AddTaskJetBackgroundSubtract(subtractBranches,1,
						    "B0","B%d");
	taskSubtract->SetBackgroundBranch(bkgBranch);
	if(!jetAOD.IsNull()) taskSubtract->SetNonStdOutputFile(jetAOD.Data());
      }
    }
#endif
  }
}



//====================================================================
/** 
 * Create the input chain
 * 
 * 
 * @return Pointer to newly allocated train 
 */
TChain *CreateChain()
{
  // Create the input chain
  chain = new TChain("esdTree");
  if (gSystem->AccessPathName("AliESDs.root"))
    ::Error("AnalysisTrainNew.C::CreateChain", 
	    "File: AliESDs.root not in ./data dir");
  else
    chain->Add("AliESDs.root");
  if (chain->GetNtrees()) return chain;
  return NULL;
}

//====================================================================
/** 
 * Merge AOD output 
 * 
 */
void AODMerge()
{
  // Merging method. No staging and no terminate phase.
  TStopwatch  timer; timer.Start();
  TString     outputDir     = "wn.xml";
  TString     outputFiles   = ("EventStat_temp.root," +
			       "AODQA.root," +
			       "AliAOD.root," +
			       "AliAOD.VertexingHF.root," +
			       "AliAOD.Muons.root," +
			       "AliAOD.Jets.root," +
			       "pyxsec_hists.root");
  TString     mergeExcludes = "";
  TObjArray*  tokens        = outputFiles.Tokenize(",");
  TIter       iter(tokens);
  TObjString* str           = 0;
  TString     outputFile;
  Bool_t      merged = kTRUE;
  while ((str = static_cast<TObjString*>(iter()))) {
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
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  timer.Print();
}

//====================================================================
/** 
 * Set-up and run AOD train 
 * 
 * @param run      Run number 
 * @param xmlfile  Collection 
 * @param stage    Stage 
 */
void AOD(UInt_t run, const char* xmlfile=0, Int_t stage=0)
{
  // -----------------------------------------------------------------
  // 
  // Get GRP parameters.  Defines global "grp" as a pointer to GRPData
  //
  gROOT->Macro(Form("GRP.C(%d)", run));
  gROOT->Macro("AODConfig.C");

  // --- Some settings -----------------------------------------------
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

  // --- Friendly message --------------------------------------------
  printf("===================================================\n");
  printf("===========    RUNNING FILTERING TRAIN   ==========\n");
  printf("===================================================\n");
  printf("=  Configuring analysis train for:\n");
  if (aodCfg->UsePhysicsSelection())   printf("=  Physics selection");
  if (aodCfg->UseTender())             printf("=  TENDER");
  if (aodCfg->UseESDfilter())            printf("=  ESD filter");
  if (aodCfg->UseMUONcopyAOD())          printf("=  MUON copy AOD");
  if (aodCfg->UseJETAN())                printf("=  Jet analysis");
  if (aodCfg->UseJETANdelta())           printf("=     Jet delta AODs");
  if (aodCfg->UsePWGHFvertexing())       printf("=  PWGHF vertexing");
  if (aodCfg->UsePWGDQJPSIfilter())      printf("=  PWGDQ j/psi filter");
  if (aodCfg->UsePWGHFd2h())             printf("=  PWGHF D0->2 hadrons QA");

  // Load common libraries and set include path
  if (!LoadCommonLibraries()) {
    ::Error("AnalysisTrain", "Could not load common libraries");
    return;
  }

  // === Make the analysis manager and connect event handlers ========
  // 
  // --- Analysis manager and load libraries -------------------------
  AliAnalysisManager *mgr = new AliAnalysisManager("Filter","Production train");
  if (aodCfg->UseSysInfo()) mgr->SetNSysInfo(100);
  if (!LoadAnalysisLibraries()) {
    ::Error("AnalysisTrain", "Could not load analysis libraries");
    return;
  }

  // --- Create ESD input handler ------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);

  // --- Monte Carlo handler -----------------------------------------
  if (aodCfg->UseMC()) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetPreReadMode(1);
    mcHandler->SetReadTR(aodCfg->UseTR());
  }

  // --- AOD output handler ------------------------------------------
  if (aodCfg->UseAODhandler()) {
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("AliAOD.root");
    mgr->SetOutputEventHandler(aodHandler);
  }

  // --- Debugging if needed -----------------------------------------
  if (aodCfg->UseDBG()) mgr->SetDebugLevel(3);

  // === Set up tasks ================================================
  //
  // --- Create tasks ------------------------------------------------
  AddAnalysisTasks(cdbPath);

  // --- If merging, do so here and exit -----------------------------
  if (stage > 0) {
    AODMerge();
    mgr->InitAnalysis();
    mgr->SetGridHandler(new AliAnalysisAlien);
    mgr->StartAnalysis("gridterminate",0);
    return;
  }
  // === Run the analysis ============================================
  //
  // --- Make our chain ----------------------------------------------
  TChain *chain = CreateChain();
  if (!chain) return;

  // --- Run the thing -----------------------------------------------
  TStopwatch timer;
  timer.Start();
  mgr->SetSkipTerminate(kTRUE);
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  timer.Print();
}

// 
// EOF
// 

