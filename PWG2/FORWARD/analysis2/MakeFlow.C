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
	      Int_t   etabins   = 40,
	      Int_t   zVertex   = 2,
	      TString addFlow   = "",
              Int_t   addFType  = 0,
              Int_t   addFOrder = 0)
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

  // --- Add to chain either AOD ------------------------------------
  if (data.Length() <= 1) {
    AliError("You didn't add a data file");
    return;
  }
  TChain* chain = new TChain("aodTree");

  if (data.Contains(".txt"))
    MakeChain(data, chain);

  if (data.Contains(".root")) {
    if (!TFile::Open(data.Data())) {
      AliError(Form("AOD file %s not found", data.Data()));
      return;
    }
    chain->Add(data.Data());
  }

  // --- Initiate the event handlers --------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward Flow", 
						    "Flow in the forward region");
  mgr->SetUseProgressBar(kTRUE, 10);

  // --- AOD input handler -------------------------------------------
  AliAODInputHandler *aodInputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodInputHandler);      

  // --- Create output handler ---------------------------------------
  AliAODHandler* aodOut = new AliAODHandler();
  mgr->SetOutputEventHandler(aodOut);
  aodOut->SetOutputFileName("AliAODs.Flow.root");

  // --- Add the tasks ---------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskForwardFlow.C");
  AddTaskForwardFlow(type, etabins, zVertex, addFlow, addFType, addFOrder);

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
