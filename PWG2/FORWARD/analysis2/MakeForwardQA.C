/**
 * @file   MakeForwardQA.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:08:14 2011
 * 
 * @brief  Generate energy loss fits 
 * 
 * @ingroup pwg2_forward_scripts_makers
 */
/** 
 * Run a pass on ESD data to produce the energ loss fits 
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
 * @param cent       Whether to use centrality or not 
 * @param name       Name of train - free form.  This will be the name
 *                   of the output directory if the plug-in is used 
 *
 * @ingroup pwg2_forward_eloss
 */
void MakeForwardQA(const char* esddir, 
		   Int_t       nEvents = 1000, 
		   Int_t       proof   = 0,
		   Bool_t      mc      = false,
		   Bool_t      cent    = false,
		   const char* name    = 0)
{
  // --- Possibly use plug-in for this -------------------------------
  if ((name && name[0] != '\0') && gSystem->Load("libRAliEn") >= 0) {
    const char* builder = 
      "$(ALICE_ROOT)/PWG2/FORWARD/analysis2/trains/BuildTrain.C";
    gROOT->LoadMacro(builder);

    BuildTrain("MakeQATrain");

    MakeQATrain t(name, cent, false);
    t.SetDataDir(esddir);
    t.SetDataSet("");
    t.SetProofServer(Form("workers=%d",proof));
    t.Run(proof > 0 ? "PROOF" : "LOCAL", "FULL", nEvents, mc, proof > 0);
    return;
  }

  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof > 0) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    LoadPars(proof);
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD",esddir, true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();
  Info("MakeForwardQA", "Will analyse %d events", nEvents);

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward QA Train", 
						    "Forward QA");
  AliAnalysisManager::SetCommonFileName("forward_qa.root");

  // --- ESD input handler -------------------------------------------
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);      
       
  // --- Monte Carlo handler -----------------------------------------
  if (mc) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    mcHandler->SetReadTR(true);    
  }

  // --- Add tasks ---------------------------------------------------
  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(mc, kTRUE, kFALSE);
  
  // Centrality 
  if(cent) {
    gROOT->LoadMacro("AddTaskCentrality.C");
    AddTaskCentrality();
  }
  
  // FMD ELoss fitter 
  gROOT->LoadMacro("AddTaskForwardQA.C");
  AddTaskForwardQA(mc, cent);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && proof <= 0) 
    mgr->SetUseProgressBar(kTRUE,100);

  // Run the train 
  t.Start();
  Printf("=== RUNNING ANALYSIS on %9d events ========================",nEvents);
  mgr->StartAnalysis(proof > 0 ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
