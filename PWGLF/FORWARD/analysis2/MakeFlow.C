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
 * @ingroup pwg2_forward_flow
 */
void MakeFlow(TString data      = "", 
	      Int_t   nevents   = 0, 
	      TString type      = "", 
	      Int_t   etabins   = 48,
              Bool_t  mc        = kFALSE,
	      TString addFlow   = "",
              Int_t   addFType  = 0,
              Int_t   addFOrder = 0,
	      Bool_t  proof     = kFALSE)
{
  Bool_t proof = kFALSE;

  // --- Load libs ---------------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof> 0) { 
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    if (!LoadPars(proof)) { 
      AliError("MakeFlow", "Failed to load PARs");
      return;
    }
  }

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Add to chain either AOD ------------------------------------
  if (data.IsNull()) {
    AliError("You didn't add a data file");
    return;
  }
  TChain* chain = new TChain("aodTree");

  if (data.Contains(".txt")) MakeChain(data, chain);

  if (data.Contains(".root")) {
    TFile* test = TFile::Open(data.Data());
    if (!test) {
      AliError(Form("AOD file %s not found", data.Data()));
      return;
    }
    test->Close(); // Remember to close!
    chain->Add(data.Data());
  }

  // --- Initiate the event handlers --------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Flow", 
						    "Flow in the forward region");

  // --- AOD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler); 

  // --- Add the tasks ---------------------------------------------
  gROOT->LoadMacro("AddTaskForwardFlow.C");
  AddTaskForwardFlow(type, etabins, mc, addFlow, addFType, addFOrder);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeFlow", "Failed to initialize analysis train!");
    return;
  }
  mgr->PrintStatus();
  Printf("****************************************");
  Printf("Doing flow analysis on %d Events", nevents == 0 ? chain->GetEntries() : nevents);
  Printf("****************************************");
  // 
  if (proof) mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !proof) 
    mgr->SetUseProgressBar(kTRUE,chain->GetEntries() < 10000 ? 100 : 1000);

//  mgr->SetSkipTerminate(true);

  t.Start();
  if (nevents == 0) mgr->StartAnalysis("local", chain);
  if (nevents != 0) mgr->StartAnalysis("local", chain, nevents);
  t.Stop();
  t.Print();
}
//----------------------------------------------------------------
void MakeChain(TString data = "", TChain* chain = 0)
{
  // creates chain of files in a given directory or file containing a list.
  
  // Open the input stream
  ifstream in;
  in.open(data.Data());

  // Read the input list of files and add them to the chain
  TString line;
  TFile* file;
  while(in.good()) 
  {
    in >> line;
      
    if (line.Length() == 0)
      continue;      
    if (!(file = TFile::Open(line))) 
      gROOT->ProcessLine(Form(".!rm %s", line.Data()));
    else {
      chain->Add(line);
      file->Close();
    }
  }

  in.close();

  return;
}
//
// EOF
//
