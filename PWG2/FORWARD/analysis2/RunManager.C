/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Flags for the analysis
 *
 * @ingroup pwg2_forward_scripts
 */
enum {
  kMC        = 0x01, // MC input 
  kProof     = 0x02, // Run in proof mode 
  kFull      = 0x04, // Do full analysis - including terminate 
  kAnalyse   = 0x08, // Run the analysis 
  kTerminate = 0x10  // Run only terminate 
};

/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * @param esddir    ESD input directory. Any file matching the pattern 
 *                  *AliESDs*.root are added to the chain 
 * @param nEvents   Number of events to process.  If 0 or less, then 
 *                  all events are analysed
 * @param flags     Job flags. A bit mask of 
 *  - 0x01 (MC)        Monte-Carlo truth handler installed 
 *  - 0x02 (PROOF)     Proof mode 
 *  - 0x04 (FULL)      Run full analysis - including terminate 
 *  - 0x08 (ANALYSE)   Run only analysis - not terminate 
 *  - 0x10 (TERMINATE) Run no analysis - just terminate.  
 * 
 * of these, PROOF, FULL, ANALYSE, and TERMINATE are mutually exclusive. 
 *
 * If PROOF mode is selected, then Terminate will be run on the master node 
 * in any case. 
 * 
 * If FULL is selected, then the full analysis is done.  Note that
 * this can be combined with PROOF but with no effect.
 *
 * ANALYSE cannot be combined with FULL, PROOF, or TERMINATE.  In a
 * local job, the output AnalysisResults.root will still be made, but
 * terminate is not called.
 *
 * In TERMINATE, the file AnalysisResults.root is opened and all
 * containers found are connected to the tasks.  The terminate member
 * function is then called
 * 
 *
 * @ingroup pwg2_forward_scripts
 */
void RunManager(const char* esddir, 
		Int_t       nEvents=1000, 
		UShort_t    flags=kFull)
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (flags & kProof) 
    gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeESDChain.C");
  TChain* chain = MakeESDChain(esddir);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", 
						    "FMD analysis train");

  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("AliESDACORDE "
				  "AliRawDataErrorLogs "
				  "CaloClusters "
				  "Cascades "
				  "EMCALCells "
				  "EMCALTrigger "
				  "Kinks "
				  "Cascades "
				  "MuonTracks "
				  "TrdTracks "
				  "CaloClusters "
				  "HLTGlobalTrigger");
  mgr->SetInputEventHandler(esdHandler);      
       
  // Monte Carlo handler
  if (flags & kMC) { 
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(readTR);    
  }
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");

  // --- Add tasks ---------------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliAnalysisTask* task = AddTaskFMD();
  // mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  task = AddTaskPhysicsSelection((flags & kMC), kTRUE, kTRUE);
  // mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Skip terminate if we're so requested and not in Proof or full mode
  mgr->SetSkipTerminate(!(flags & kProof) &&
			!(flags & kFull)  && 
			(flags & kAnalyse));
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !(flags & kProof)) 
    mgr->SetUseProgressBar(kTRUE);

  // Run the train 
  t.Start();
  if (!(flags & kTerminate)) {
    Printf("=== RUNNING ANALYSIS ==================================");
    mgr->StartAnalysis((flags & kProof) ? "proof" : "local", chain, nEvents);
  }
  else {
    Printf("=== RUNNING TERMINATE =================================");
    // mgr->ImportWrappers(0);
    mgr->Terminate();
  }
  t.Stop();
  t.Print();
}
//
// EOF
//
