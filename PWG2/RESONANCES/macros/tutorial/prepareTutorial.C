//
// *** prepareTutorial.C ***
//
// This example macro does no analysis.
// It is useful to clean our PROOF session from all old
// versions of all required support libraries,
// and to check if the connection to CAF works properly.
//
// The libraries required for analysis are: STEERBase, ESD, AOD, ANALYSIS and ANALYSISalice,
// plus the core package library which is PWG2resonances
// This macro cleans previous compilation of them, and recompile each library.
// Then, the compiled libraries are loaded into PROOF.
//
void prepareTutorial(const char *proofToConnect = 0) {

  // check argument
  if (!proofToConnect) {
    cerr << "This macro requires an argument: <userID>@alicecaf.cern.ch" << endl;
    cerr << "You can edit the macro and add it as default value, " << endl;
    cerr << "or you must run it as 'root -q -b 'prepareTutorial.C(\"<userID>@alicecaf.cern.ch\")'" << endl;
    return;
  }

  // Boolean if the process is running without errors
  Bool_t isProcessOK = kTRUE;

  // Sets up timer
  TStopwatch timer;
  timer.Start();

  // Loads "PWG2resonancesUtils.C" macro, where AliRsnUtils class is defined
  gROOT->LoadMacro("PWG2resonancesUtils.C");

  // Setting up the proof connection (in this example we will connect to 'alicecaf.cern.ch' with 'mvala' username)
  TString proof = proofToConnect;

  // Seting up the par-files we will use(without .par extension) and divided by ':'
  TString pars = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice:PWG2resonances";

  // Create utils object which will connect to the proof, loads libraries
  // and in next tutorials it will sets up inpit files and run analysis
  AliRsnUtils *utils = new AliRsnUtils(AliRsnUtils::kProof,proof);

  // Cleans all package (leave it uncommented only if you have problem to load parfiles in proof)
  // so in future tutorial we will comment it, because we don't want to recompile same parfile all the time
  utils->CleanPackages();

  // Loading pars
  isProcessOK = utils->LoadPars(pars);

  // if one of the par-file was not load correctly we will stop the macro
  if (!isProcessOK) return;

  // Prints process time
  timer.Stop();
  timer.Print();
}