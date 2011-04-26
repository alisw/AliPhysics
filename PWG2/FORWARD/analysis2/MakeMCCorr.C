/**
 * @file   MakedNdeta.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 09:41:56 2011
 * 
 * @brief  Run second pass analysis - make @f$ dN/d\eta@f$
 * 
 * @ingroup pwg2_forward_scripts_makers
 */
/** 
 * Run second pass analysis - make @f$ dN/d\eta@f$
 * 
 * If the ROOT AliEn interface library (libRAliEn) can be loaded, 
 * and the parameter @a name is not empty, then use the plugin to do
 * the analysis.  Note that in this case, the output is placed 
 * in a sub-directory named by @a name after escaping spaces and special 
 * characters 
 * 
 * @param aoddir     AOD input directory. Any file matching the pattern 
 *                   *AliAODs*.root are added to the chain 
 * @param nEvents    Number of events to process.  If 0 or less, then 
 *                   all events are analysed
 * @param trig       Trigger to use 
 * @param useCent    Whether to use centrality or not 
 * @param scheme     Normalisation scheme 
 * @param vzMin      Least @f$ v_z@f$ (centimeter)
 * @param vzMax      Largest @f$ v_z@f$ (centimeter)
 * @param proof      If larger then 1, run in PROOF-Lite mode with this 
 *                   many number of workers. 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 *
 * @ingroup pwg2_forward_dndeta
 */
void MakeMCCorr(const char* esddir   = ".", 
	        Int_t       nEvents  = -1, 
		Double_t    vzMin    = -10,
		Double_t    vzMax    = +10,
	        Int_t       proof    = 0,
		const char* name     = 0)
{
  if ((name && name[0] != '\0') && gSystem->Load("libRAliEn") >= 0) {
    Error("MakeMCCorr", "Plug-in mode not implemented yet!");
    return;

    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			     "$ALICE_ROOT/ANALYSIS/macros",
			     gROOT->GetMacroPath()));
    gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gROOT->LoadMacro("TrainSetup.C+");
    MakeMCCorrTrain t(name, vzMin, vzMax);
    t.SetDataDir(esddir);
    t.SetDataSet("");
    t.SetAllowOverwrite(true);
    t.SetProofServer(Form("workers=%d",proof));
    t.Run(proof > 0 ? "PROOF" : "LOCAL", "FULL", nEvents, proof > 0);
    return;
  }
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    LoadPars(proof);
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD", esddir,true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager(name, "Forward MC corr");
  AliAnalysisManager::SetCommonFileName("forward_mccorr.root");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);      
       
  // --- Monte Carlo handler -----------------------------------------
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mcHandler->SetReadTR(true);    

  // --- Add Physics Selection ---------------------------------------
  // Physics selection 
  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(kTRUE, kTRUE, kFALSE);
  // --- Fix up physics selection to give proper A,C, and E triggers -
  AliInputEventHandler* ih =
    static_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  // Ignore trigger class when selecting events.  This mean that we
  // get offline+(A,C,E) events too
  // ps->SetSkipTriggerClassSelection(true);

  // --- Add our task ------------------------------------------------
  gROOT->LoadMacro("AddTaskForwardMCCorr.C");
  AddTaskForwardMCCorr();

  gROOT->LoadMacro("AddTaskCentralMCCorr.C");
  AddTaskCentralMCCorr();
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakedNdeta", "Failed to initialize analysis train!");
    return;
  }
  // Skip terminate if we're so requested and not in Proof or full mode
  mgr->SetSkipTerminate(false);
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !proof) 
    mgr->SetUseProgressBar(kTRUE,100);
  
  // Run the train 
  t.Start();
  Printf("=== RUNNING ANALYSIS ==================================");
  mgr->StartAnalysis(proof ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
