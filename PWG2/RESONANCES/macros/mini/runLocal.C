void runLocal
(
   Int_t       nReadFiles      =  0,
   Int_t       nSkipFiles      =  0,
   Int_t       nmix            = 10,
   //const char *inputSource     = "file-collections/AOD048_LHC11a10b.txt",
   //const char *options         = "aod",
   //const char *inputSource     = "file-collections/AOD049_LHC10h_pass2.txt",
   //const char *options         = "aod",
   const char *inputSource     = "file-collections/ESD_LHC10d1.txt",
   const char *options         = "esd_mc",
   const char *outName         = "test.root",
   const char *macroPath       = ".",
   const char *setupName       = "AnalysisSetupRsnMini.C"
)
{
   //
   // === PREPARATION ==============================================================================
   //
   
// AliLog::SetClassDebugLevel("AliRsnMiniOutput"      , 2);
// AliLog::SetClassDebugLevel("AliRsnMiniAnalysisTask", 1);
   
   // execute the general setup from the apposite macro
   // it returns also a TString value with the input tree name
   gROOT->LoadMacro(setupName);
   TString out = Setup(nmix, options, outName, macroPath);
   cout << "OUT = " << out.Data() << endl;
   if (out.Length() < 1) return;

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
   if (isTXT) analysisChain = CreateChainFromText(inputSource, out.Data(), nReadFiles, nSkipFiles);
   else       analysisChain = CreateChainFromXML (inputSource, out.Data(), nReadFiles, nSkipFiles);
   if (!analysisChain) {
      Error("runLocal", "Analysis chain not properly initialized");
      return;
   }
   
   //
   // === ANALYSIS TASK CREATION AND INCLUSION =====================================================
   //
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return;
   
   // initialize and start analysis
   if (!mgr->InitAnalysis()) {
      ::Error("runLocal.C", "Failed to init analysis");
      return;
   }
   mgr->PrintStatus();
   if (isTXT) mgr->StartAnalysis("local", analysisChain);
   else       mgr->StartAnalysis("alien", analysisChain);
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

//__________________________________________________________________________________________________
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
