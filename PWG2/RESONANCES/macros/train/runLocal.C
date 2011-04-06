//
// This is an example steering macro for running RSN analysis task
// locally with a collection of files specified in a text file:
//
// Inputs:
//   - nReadFiles  = number of files to process from the list
//   - nSkipFiles  = how many lines to be skipped when reading the list
//   - addTaskName = name of the macro to add the RSN analysis task
//                   (assumed to have inside it a function named like the file)
//   - inputSource = name of the file containing all the inputs
//                   ---> to run on a local collection, the collection file
//                        must contain on each line the full path
//                        of one input file and it must have the ".txt" extension
//                   ---> to run on an AliEn collection, the collection file must be an XML
//                        file collection like those built from the "find -x" method in aliensh.
//   - options     = a label which is used to know what kind of data are being read
//                   (it is propagated to the 'addTask' macro for eventual setting up of something
//   - outName     = name for the file with RSN package outputs (without ROOT extension)
//
// Notes:
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
//
//
// In principle, the user should never modify this macro.
//
void runLocal
(
   Int_t       nReadFiles      = 0,
   Int_t       nSkipFiles      = 0,
   Int_t       nmix            = 0,
   const char *inputSource     = "pbpb_data.txt",
   const char *runOptions      = "esd_data_phys_cent",
   const char *analysisOptions = "tpcpid_tofpid_mult",
   const char *outName         = "test.root",
   const char *taskList        = "AddRsnAnalysisTask.C",
   //const char *taskPath        = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train"
   const char *taskPath        = "$(HOME)/code/resonances/alice-rsn-package/PWG2resonances/RESONANCES/macros/test/pulvir"
)
{
   //
   // === PREPARATION ==============================================================================
   //
   
   // this option is not needed when using plugin
   // gEnv->SetValue("XSec.GSI.DelegProxy","2");
   
   // stopwatch
   TStopwatch timer;
   timer.Start();
   
   // some options
   TString opt(runOptions);
   opt.ToUpper();
   Bool_t useTender = opt.Contains("TENDER");
   Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
   Bool_t isESD     = opt.Contains("ESD");
   ::Info("runPlugin.C", "useTender = %d", useTender);
   ::Info("runPlugin.C", "isMC      = %d", isMC     );
   ::Info("runPlugin.C", "runOpts   = %s", runOptions);
   ::Info("runPlugin.C", "anaOpts   = %s", analysisOptions);
   
   // basic configurations
   gROOT->LoadMacro(Form("%s/AnalysisSetup.C", taskPath));
   AnalysisSetup(isMC, nmix, runOptions, outName, taskPath);
   
   // define tree name
   Char_t treeName[200];
   if (isESD) sprintf(treeName, "esdTree"); else sprintf(treeName, "aodTree");
   
   //
   // === BUILD INPUT LIST =========================================================================
   //

   // check extension of input to distinguish between XML and TXT
   TString sInput(inputSource);
   sInput.ToLower();
   Bool_t isTXT = (!strcmp(sInput(sInput.Length() - 3, 3).Data(), "txt"));
   cout << "Input = " << (isTXT ? "TXT" : "XML") << endl;

   // if input is XML, connect to AliEn
   if (!isTXT) TGrid::Connect("alien://");

   // create TChain of input events
   TChain *analysisChain = 0x0;
   if (isTXT) analysisChain = CreateChainFromText(inputSource, treeName, nReadFiles, nSkipFiles);
   else       analysisChain = CreateChainFromXML (inputSource, treeName, nReadFiles, nSkipFiles);
   if (!analysisChain) {
      Error("runLocal", "Analysis chain not properly initialized");
      return;
   }
   
   //
   // === CONFIGURATION ============================================================================
   //

   //AliLog::SetGlobalDebugLevel(AliLog::kDebug + 3);
   //AliLog::SetClassDebugLevel("AliRsnCutESD2010", AliLog::kDebug+3);
   //AliLog::SetClassDebugLevel("AliRsnCutValue", AliLog::kDebug+3);
   //AliLog::SetClassDebugLevel("AliRsnCutESDCutMultiplicity", AliLog::kDebug+3);
   //AliLog::SetClassDebugLevel("AliRsnValue", AliLog::kDebug+3);
   //AliLog::SetClassDebugLevel("AliRsnCutTrackQuality", AliLog::kDebug+3);
   //AliLog::SetGlobalDebugLevel(AliLog::kDebug+3);
   
   //
   // === ANALYSIS TASK CREATION AND INCLUSION =====================================================
   //
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return;
   
   gROOT->LoadMacro(Form("%s/AddRsnAnalysisTask.C", taskPath));
   AddRsnAnalysisTask(isMC, (nmix > 0), runOptions, analysisOptions);

   // initialize and start analysis
   if (!mgr->InitAnalysis()) {
      ::Error("runPlugin.C", "Failed to init analysis");
      return;
   }
   mgr->PrintStatus();
   if (isTXT) mgr->StartAnalysis("local", analysisChain);
   else       mgr->StartAnalysis("alien", analysisChain);
   
   // final operation
   // gObjectTable->Print();
   timer.Stop();
   timer.Print();
}

//_________________________________________________________________________________________________
TChain* CreateChainFromXML
(const char *xmlFileName, const char *treeName, Int_t nread, Int_t nskip)
{
//
// Create a TChain with all required files listed into an XML collection.
// Necessary to run analysis in AliEn jobs.
// ---
// Arguments:
//  - xmlFileName = input list
//  - treeName    = "esdTree" or "aodTree"
//  - nread       = how many files to read (0 = all)
//  - nskip       = how many files to skip from beginning
//

   // if nread argument is 0, it is disabled
   if (nread == 0) nread = 1000000000;

   // initialize output object
   TChain *chain = new TChain(treeName);

   // initialize the AliEn collection
   TAlienCollection *myCollection = TAlienCollection::Open(xmlFileName);
   if (!myCollection) {
      Error("CreateChainFromXML", "Cannot create an AliEn collection from %s", xmlFileName);
      return 0x0;
   }

   // loop on collection
   myCollection->Reset();
   while (myCollection->Next()) {
      // skip until reached required number of offset
      if (nskip > 0) {--nskip; continue;}

      // stop if required number of read files is reached
      // otherwise update the counter
      if (nread <= 0) break;
      nread--;

      // recovery file and add it
      Info("CreateChainFromXML", Form("Adding: %s", myCollection->GetTURL("")));
      chain->Add(myCollection->GetTURL(""));
   }

   return chain;
}

//_________________________________________________________________________________________________
TChain* CreateChainFromText(const char *fileName, const char *treeName, Int_t nread, Int_t nskip)
{
//
// Create a TChain with all required files listed into a text file.
// Necessary to run analysis in local jobs.
// ---
// Arguments:
//  - xmlFileName = input file list
//  - treeName    = "esdTree" or "aodTree"
//  - nread       = how many files to read (0 = all)
//  - nskip       = how many files to skip from beginning
//

   // if third argument is 0, it is interpreted
   // as "read all lines"
   Bool_t readAll = (nread <= 0);

   // initialize output object
   TChain* target = new TChain(treeName);

   // open text file
   ifstream fileIn(fileName);

   // loop on collection
   TString line;
   while (fileIn.good()) {
      fileIn >> line;
      if (line.IsNull()) continue;

      // skip until reached required number of offset
      if (nskip > 0) {--nskip; continue;}

      // stop if required number of read files is reached
      // otherwise update the counter
      if (!readAll && nread <= 0) break;
      nread--;

      // add file
      Info("CreateChainFromText", "Adding '%s'", line.Data());
      target->Add(line.Data());
   }

   return target;
}
