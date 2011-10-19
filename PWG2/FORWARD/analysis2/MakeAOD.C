/** 
 * @defgroup pwg2_forward_scripts_makers Maker scripts 
 * @ingroup pwg2_forward_scripts
 */
//====================================================================
/**
 * @file   MakeAOD.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 09:40:10 2011
 * 
 * @brief  Run first pass of the analysis - AOD generation
 * 
 * @ingroup pwg2_forward_scripts_makers
 */
//====================================================================
/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * If the ROOT AliEn interface library (libRAliEn) can be loaded, 
 * and the parameter @a name is not empty, then use the plugin to do
 * the analysis.  Note that in this case, the output is placed 
 * in a sub-directory named by @a name after escaping spaces and special 
 * characters 
 *
 * If PROOF mode is selected, then Terminate will be run on the master node 
 * in any case. 
 *
 * @param esddir     ESD input directory. Any file matching the pattern 
 *                   *AliESDs*.root are added to the chain 
 * @param nEvents    Number of events to process.  If 0 or less, then 
 *                   all events are analysed
 * @param proof      If larger then 1, run in PROOF-Lite mode with this 
 *                   many number of workers. 
 * @param mc         Data is assumed to be from simulations  
 * @param centrality Whether to use centrality or not 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 *
 * @ingroup pwg2_forward_aod
 */
void MakeAOD(const char* esddir, 
	     Int_t       nEvents    = -1, 
	     Int_t       proof      = 0,
	     Bool_t      mc         = false,
	     Bool_t      centrality = true,
	     const char* name       = 0,
	     bool        debug      = false)
{
  // --- Possibly use plug-in for this -------------------------------
  if ((name && name[0] != '\0') && gSystem->Load("libRAliEn") >= 0) {
    const char* builder = 
      "$(ALICE_ROOT)/PWG2/FORWARD/analysis2/trains/BuildTrain.C" 
    gROOT->LoadMacro(builder);

    BuildTrain("MakeAODTrain");

    MakeAODTrain t(name, 0, 0, 0, centrality, false);
    t.SetDataDir(esddir);
    t.SetDataSet("");
    t.SetProofServer(Form("workers=%d",proof));
    t.SetUseGDB(debug);
    t.Run(proof > 0 ? "PROOF" : "LOCAL", "FULL", nEvents, mc, proof > 0);
    return;
  }

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
  
  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager(name, 
						    "Forward multiplicity");
  AliAnalysisManager::SetCommonFileName("forward.root");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
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
  aodHandler->SetNeedsHeaderReplication();
  aodHandler->SetOutputFileName("AliAOD.root");

  // --- Add tasks ---------------------------------------------------
  // Physics selection 
  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(mc, kTRUE, kFALSE);
  // --- Fix up physics selection to give proper A,C, and E triggers -
  AliInputEventHandler* ih =
    static_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  // Ignore trigger class when selecting events.  This mean that we
  // get offline+(A,C,E) events too
  // ps->SetSkipTriggerClassSelection(true);
  
  // Centrality 
  if(centrality) {
    gROOT->LoadMacro("AddTaskCentrality.C");
    AliCentralitySelectionTask* centTask = AddTaskCentrality();
    centTask->SetPass(1);
    if(mc)centTask->SetMCInput();
  }

  // Copy header information 
  gROOT->LoadMacro("AddTaskCopyHeader.C");
  AddTaskCopyHeader();

  // FMD 
  gROOT->LoadMacro("AddTaskForwardMult.C");
  AddTaskForwardMult(mc);

  // Central 
  gROOT->LoadMacro("AddTaskCentralMult.C");
  AddTaskCentralMult(mc);
  
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
  if (proof) mgr->SetDebugLevel(3);
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
