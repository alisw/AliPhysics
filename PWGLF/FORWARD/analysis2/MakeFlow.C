/**
 * @file   MakeFlow.C
 * @author Alexander Hansen 
 * @date   Wed Mar 23 12:11:33 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_scripts_makers
 * 
 */
/**
 * Script to analyse AOD input for flow
 * 
 * Takes either a single (AOD) .root file as input or a .txt
 * The .txt file is expected to contain the path to the files 
 * from the current directory or the absolute path.
 * 
 * @par Inputs: 
 *  
 * 
 * @par Outputs: 
 * - 
 *
 * @ingroup pwglf_forward_flow
 */
void MakeFlow(const char* data  = "", 
	      Int_t   nEvents   = -1, 
	      TString type      = "", 
              Bool_t  mc        = kFALSE,
	      const char* name  = 0,
	      Int_t   proof     = 0,
	      Bool_t  dispVtx   = kFALSE,
	      TString addFlow   = "",
              Int_t   addFType  = 0,
              Int_t   addFOrder = 0,
              Bool_t  gdb       = kFALSE)
{
  // --- Load libs ---------------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");

  Int_t debug = 0;

  // --- Possibly use plug-in for this -------------------------------
  if ((name && name[0] != '\0') && gSystem->Load("libRAliEn") >= 0) {
 
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/TrainSetup.C+");
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/MakeFlowTrain.C+");

    char* tok = strtok(name, " ,:");
    TString server = "", dataset = "", datadir = "", jobname = "";
    if (tok[0] == NULL) Fatal("Input name invalid!");
    if (strcmp(tok, "hehi00.nbi.dk") == 0) {
      server = tok;
      dataset = data;
      datadir = "";
      tok = strtok(NULL, " ,:");
      jobname = tok;
    } else {
      server = Form("workers=%d", proof);
      dataset = "";
      datadir = data;
      jobname = tok;
    }
    printf("Making flow train on server: %s\t with datadir: %s\t on dataset: %s\n", 
	    server.Data(), datadir.Data(), dataset.Data());

    MakeFlowTrain t(jobname.Data(), type.Data(), mc, addFlow.Data(), addFType, addFOrder, false);
    t.SetUseDispVtx(dispVtx);
    t.SetProofServer(server.Data());
    t.SetDataDir(datadir.Data());
    t.SetDataSet(dataset.Data());
    t.SetDebugLevel(debug);
    t.SetUseGDB(gdb);
    t.Run("proof", "full", nEvents, proof > 0);

    return;
  }

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Add to chain either AOD ------------------------------------
  if (data == '\0') {
    AliError("You didn't add a data file");
    return;
  }
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("AOD", data, true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Initiate the event handlers --------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Flow", 
						    "Flow in the forward region");

  // --- AOD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler); 

  // --- Add the tasks ---------------------------------------------
  gROOT->LoadMacro("AddTaskForwardFlow.C");
  AddTaskForwardFlow(type, mc, dispVtx, addFlow, addFType, addFOrder);
  mgr->SetDebugLevel(debug);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeFlow", "Failed to initialize analysis train!");
    return;
  }
  mgr->PrintStatus();
  Printf("****************************************");
  Printf("Doing flow analysis on %d Events", nEvents);
  Printf("****************************************");
  // 
  if (mgr->GetDebugLevel() < 1) 
    mgr->SetUseProgressBar(kTRUE, nEvents < 10000 ? 100 : 1000);

//  mgr->SetSkipTerminate(true);

  t.Start();
  mgr->StartAnalysis("local", chain, nEvents);
  t.Stop();
  t.Print();
}
//----------------------------------------------------------------
//
// EOF
//
