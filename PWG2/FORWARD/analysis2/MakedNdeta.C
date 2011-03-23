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
 *
 * @ingroup pwg2_forward_dndeta
 */
void MakedNdeta(const char* aoddir   = ".", 
	        Int_t       nEvents  = -1, 
		const char* trig     = "INEL",
		Bool_t      useCent  = false,
		const char* scheme   = 0,
		Double_t    vzMin    = -10,
		Double_t    vzMax    = +10,
	        Int_t       proof    = 0)
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    LoadPars(proof);
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("AOD", aoddir,true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Train", 
						    "Forward dN/deta");
  AliAnalysisManager::SetCommonFileName("forward_dndeta.root");

  // --- ESD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler);      
       
  // --- Add tasks ---------------------------------------------------
  // Forward 
  gROOT->LoadMacro("AddTaskForwarddNdeta.C");
  AddTaskForwarddNdeta(trig, vzMin, vzMax, useCent, scheme, true);
  // Central
  gROOT->LoadMacro("AddTaskCentraldNdeta.C");
  AddTaskCentraldNdeta(trig, vzMin, vzMax, useCent, scheme);

  
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
