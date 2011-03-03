/**
 * Script to analyse AOD input for flow 
 * 
 * @par Inputs: 
 * - 
 * 
 * @par Outputs: 
 * - 
 * 
 */
void MakeFlow(TString path      = "", 
	      Bool_t  recursive = true,
	      Int_t   nevents   = 100, 
	      TString type      = "", 
	      Int_t   etabins)
{
  Bool_t proof = false;

  // --- Load libs -------------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    if (!LoadPars(proof)) { 
      Error("MakeAOD", "Failed to load PARs");
      return;
    }
  }

  // --- Add to chain either ESD or AOD ----------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("AOD", path, recursive);
  // If 0 or less events is select, choose all 
  if (nevents <= 0) nevents = chain->GetEntries();

  // --- Initiate the event handlers -------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Flow", 
						    "Flow in forward region");
  mgr->SetUseProgressBar(kTRUE);

  // --- AOD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler);      

  // --- Create output handler ---------------------------------------
  AliAODHandler* aodOut = new AliAODHandler();
  mgr->SetOutputEventHandler(aodOut);
  aodOut->SetOutputFileName("AliAODs.Flow.root");

  // --- Add the tasks ---------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskForwardFlow.C");
  AddTaskForwardFlow(type.Data(), etabins);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeFlow", "Failed to initialize analysis train!");
    return;
  }
  mgr->PrintStatus();

  t.Start();
  if (nevents == 0) mgr->StartAnalysis("local", chain);
  if (nevents != 0) mgr->StartAnalysis("local", chain, nevents);
  t.Stop();
  t.Print();
}
//
// EOF
//
