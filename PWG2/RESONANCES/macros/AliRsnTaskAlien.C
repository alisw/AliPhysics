//
// Example macro to run a full resonance analysis task in a local machine
// using data taken from AliEn. This macro can be used as a template to create
// an AliEn job which runs this analysis on a collection of source data files (ESD/AOD/MC).
//
// This macro does the following:
//  - loads all required PAR libraries
//  - creates a TChain of events to be analyzed (from an XML collection)
//  - creates the AnalysisTask macro which does the job
//
// Since it is possible that all required libraries/macroes are stored in different
// directories, for each one a variable must be initialized which gives the path,
// but, in order to allow to explode (if necessary) some environment variables,
// these variables are NOT arguments, but they are hard-coded in the macro at its
// beginning, to allow a user to customize them.
//
// Arguments:
//  - name of XM collection for source files
//  - a flag to know it source data is ESD or AOD
//  - maximum number of files to be analyzed
//  - the number of files to be skipped starting from first one in collection.
//
Bool_t AliRsnTaskAlien
(
  const char *collectionFile = "wn.xml",
  Bool_t      useESDsource   = kTRUE,
  Int_t       maxFiles       = 1,
  Int_t       skipFiles      = 0
)
{
    // path for the "AliRsnLoad.C" macro which loads all libraries
    TString strLoadPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
    strLoadPath.Append("/AliRsnLoad.C");
    
    // path for the "CreateAnalysisManager.C" macro which creates the AnalysisManager
    TString strTaskPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
    strTaskPath.Append("/CreateAnalysisManager.C");
    
    // connect to grid
    TGrid::Connect("alien://");

    // load the macro for uploading packages;
    // the unique assumption which must be done here is that the PAR libraries
    // are in the working directory, which is the case in all kinds of analysis
    gROOT->LoadMacro(strLoadPath.Data());
    AliRsnLoad();
    
    // create the TChain of files to read
    // its name depends on the kind of source to be read
    char treeName[200];
    if (useESDsource) {
        sprintf(treeName, "esdTree");
        Info("AliRsnSimpleTaskAlien", "Using ESD source data");
    }
    else {
        sprintf(treeName, "aodTree");
        Info("AliRsnSimpleTaskAlien", "Using AOD source data");
    }
    TChain* analysisChain = new TChain(treeName);
    TAlienCollection *myCollection = TAlienCollection::Open(collectionFile);
    if (!myCollection) {
        Error("AliRsnSimpleTaskAlien", Form("Cannot create an AliEn collection from %s", collectionFile));
        return kFALSE;
    }
    // add files to the TChain, keeping trace of their number with a counter
    Int_t nfiles = 0;
    char fileName[255];
    myCollection->Reset();
    while (myCollection->Next()) {
        if (skipFiles) {
            skipFiles--;
            continue;
        }
        sprintf(fileName, "%s", myCollection->GetTURL(""));
        Info("AliRsnSimpleTaskAlien", Form("Adding file '%s'", fileName));
        analysisChain->Add(fileName);
        nfiles++;
        if (maxFiles > 0 && nfiles >= maxFiles) break;
    }
    Info("AliRsnSimpleTaskAlien", Form("# processed events = %d", (Int_t)analysisChain->GetEntries()));

    // load and execute the macro which creates the AnalysisManager
    // this macro is expected to define a function with the standard name
    // *** "CreateAnalysisManager()" *** [without arguments],
    // irrespectively of the filename of the macro itself,
    // which returns an AliAnalysisManager object
    gROOT->LoadMacro(strTaskPath.Data());
    AliAnalysisManager *mgr = CreateAnalysisManager(kFALSE);
    
    // initialize analysis and run it in "local" mode
    if (mgr->InitAnalysis()) {
        mgr->PrintStatus();
        return mgr->StartAnalysis("local", analysisChain);
    }
    return kTRUE;
}
