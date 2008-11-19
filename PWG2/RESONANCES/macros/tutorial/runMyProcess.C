//
// *** runMyProcess.C ***
//
// This example macro is the one which the users must run from the prompt
// for all tutorials (and it can be used event for real analysis).
// Different tutorials change in the script which configures the analysis,
// but in all cases this is the macro to be runned as front-end.
//
// It loads the PAR libraries and initialize the resonance analysis,
// by creating some tasks and setting up informations about data to be used.
//
// It is structured very similar to "prepareTutorial.C", with some more lines,
// and here it is also provided a line to be un-commented if the user wants
// to clear the PAR libraries and recompile them once more (but normally it is not needed)
//
void runMyProcess(const char *proofToConnect = 0) {

  // check argument
  TString proof("");
  if (proofToConnect) {
    proof.Append(proofToConnect);
  }
  else {
    Warning("runMyProcess.C", "No connection string specified. Retrieving user from system");
    proof = getenv("USER");
    proof.Append("@alicecaf.cern.ch");
  }
  Info("runMyProcess.C", "Connection string = \"%s\"", proof.Data());

  // Boolean if the process is running without errors
  Bool_t isProcessOK = kTRUE;

  // Sets up timer
  TStopwatch timer;
  timer.Start();

  // Loads "PWG2resonancesUtils.C" macro, where AliRsnUtils class is defined
  gROOT->LoadMacro("PWG2resonancesUtils.C");

  // Seting up the par-files we will use(without .par extension) and divided by ':'
  TString pars = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice:PWG2resonances";

  // input data type (we will use DATASET)
  AliRsnUtils::EDataType dataType = AliRsnUtils::kDataSet;
  // DATASET's name
  TString inputFileName = "/COMMON/COMMON/LHC08c11_10TeV_0.5T";
  // input tree (we will read ESD so it is esdTree)
  TString treeName = "esdTree";

  // my Resonance Analysis macro
  TString taskMacro = "myRsnAnalysis.C";
  // number of evenst to be processed
  Long64_t numEvents = 100000;
  // number of events to be skipped
  Long64_t eventsSkip = 0;

  // Create utils object which will connect to the proof, loads libraries
  // and in next tutorials it will sets up inpit files and run analysis
  AliRsnUtils *utils = new AliRsnUtils(AliRsnUtils::kProof,proof);

  // Cleans all package (leave it uncommented only if you have problem to load parfiles in proof)
  // so in future tutorial we will comment it, because we don't want to recompile same parfile all the time
  // utils->CleanPackages();

  // Loading pars
  isProcessOK = utils->LoadPars(pars);
  // if one of the par-file was not load correctly we will stop the macro
  if (!isProcessOK) return;

  // Input Data
  isProcessOK = utils->SetInputData(dataType,inputFileName,treeName);
  // if there was problem with input data we will stop the macro
  if (!isProcessOK) return kFALSE;

  // runs macro
  utils->Run(taskMacro,numEvents,eventsSkip);

  // Prints process time
  timer.Stop();
  timer.Print();
}
