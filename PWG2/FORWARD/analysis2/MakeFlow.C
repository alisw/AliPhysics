/**
 * @file   MakeFlow.C
 * @author Alexander Hansen 
 * @date   Wed Mar 23 12:11:33 2011
 * 
 * @brief Analyse AODs for flow 
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
 * @param data         Input data (directory, list, or ROOT file)
 * @param nevents      Number of events to scan 
 * @param type         Type of analysis (v<n>)
 * @param etabins      How may eta bins to make 
 * @param addFlow      If true, add flow to MC particles
 * @param addFType     Which type of flow to add to MC particles
 * @param addFOrder    Order of flow to add to MC particles
 * @param proof 
 *
 * @ingroup pwg2_forward_flow
 */
void MakeFlow(TString data      = "", 
	      Int_t   nevents   = 0, 
	      TString type      = "", 
	      Int_t   etabins   = 40,
	      TString addFlow   = "",
              Int_t   addFType  = 0,
              Int_t   addFOrder = 0,
	      Bool_t  proof     = false)
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
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("AOD", data.Data(), true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();

  // --- Initiate the event handlers --------------------------------
  AliAnalysisManager *mgr  = 
    new AliAnalysisManager("Forward Flow", 
			   "Flow in the forward region");

  // --- AOD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler);      

  // --- Add the tasks ---------------------------------------------
  gROOT->LoadMacro("AddTaskForwardFlow.C");
  AddTaskForwardFlow(type, etabins, addFlow, addFType, addFOrder);

  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("MakeFlow", "Failed to initialize analysis train!");
    return;
  }
  mgr->PrintStatus();
  // 
  if (proof) mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !proof) 
    mgr->SetUseProgressBar(kTRUE,chain->GetEntries() < 10000 ? 10 : 100);

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
  while(in.good()) 
  {
    in >> line;
      
    if (line.Length() == 0)
      continue;      
    
    if (TFile::Open(line))
      chain->Add(line);
  }

  in.close();

  return;
}
//
// EOF
//
