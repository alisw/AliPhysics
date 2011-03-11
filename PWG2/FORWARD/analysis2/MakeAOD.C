/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */

/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * @param esddir    ESD input directory. Any file matching the pattern 
 *                  *AliESDs*.root are added to the chain 
 * @param nEvents   Number of events to process.  If 0 or less, then 
 *                  all events are analysed
 * @param proof     Run in proof mode 
 * @param mc        Run over MC data
 *
 * If PROOF mode is selected, then Terminate will be run on the master node 
 * in any case. 
 * 
 *
 * @ingroup pwg2_forward_scripts
 */
void MakeAOD(const char* esddir, 
	     Int_t       nEvents=-1, 
	     Int_t       proof=0,
	     Bool_t      mc=false,
	     Bool_t      centrality=true )
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    if (!LoadPars(proof)) { 
      Error("MakeAOD", "Failed to load PARs");
      return;
    }
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD", esddir,true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Train", 
						    "Forward multiplicity");
  AliAnalysisManager::SetCommonFileName("forward.root");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches(// "AliESDRun " 
				  // "AliESDHeader "
				  // "AliESDZDC "
				  // "AliESDFMD "
				  // "AliESDVZERO " 
				  "AliESDTZERO " 
				  "TPCVertex " 
				  // "SPDVertex "
				  // "PrimaryVertex "
				  // "AliMultiplicity "
				  "PHOSTrigger "
				  "EMCALTrigger "
				  "SPDPileupVertices " 
				  "TrkPileupVertices " 
				  // "Tracks "
				  "MuonTracks " 
				  "PmdTracks "
				  "TrdTracks "
				  "V0s " 
				  "Cascades " 
				  "Kinks " 
				  "CaloClusters "
				  "EMCALLCells "
				  "PHOSCells "
				  "AliRawDataErrorLogs "
				  "ALIESDCACORDE " 
				  "HLTGlobalTrigger");
  mgr->SetInputEventHandler(esdHandler);      
       
  // --- Monte Carlo handler -----------------------------------------
  if (mc) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(true);    
  }

  // --- AOD output handler ------------------------------------------
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAOD.root");

  // --- Add tasks ---------------------------------------------------
  // Physics selection 
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  // test HHD
  AddTaskPhysicsSelection(mc, kTRUE, kFALSE);

#if 0
  // Centrality 
  if (!proof) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/Compile.C");
    Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskCopyHeader.C","");
    AddTaskCopyHeader();
  }
#endif
  if(centrality) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AddTaskCentrality();
  }
  // FMD 
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  AddTaskFMD(mc);

  // Central 
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskCentral.C");
  AddTaskCentral();
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeAOD", "Failed to initialize analysis train!");
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
  Printf("=== RUNNING ANALYSIS on %9d events ==========================",
	 nEvents);
  mgr->StartAnalysis(proof ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
