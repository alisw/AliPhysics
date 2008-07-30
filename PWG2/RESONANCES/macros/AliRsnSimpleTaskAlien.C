//
// Example macro to run a full resonance analysis task in AliEn
// This macro does the following:
//  - loading the RSN package library and all other libraries required for it
//  - creating the TChain of events to be analyzed from an XML collection
//  - loading and executing the "universal" AnalysisTask macro which does the job
//
// Since it is possible that all required libraries/macroes are stored in different
// directories, for each one a variable must be initialized which gives the path,
// but, in order to allow to explode (if necessary) some environment variables,
// these variables are NOT arguments, but they are hard-coded in the macro at its
// beginning, to allow a user to customize them.
//
// The unique argument is the collection file name, which is always in the working dir,
// and the options which tell what kind of source data is being used or other working options.
//

Bool_t AliRsnSimpleTaskAlien
(
    const char *collectionFile = "wn.xml",
    const char *configFile     = "RsnConfig.C",
    Bool_t      useESDsource   = kTRUE,
    Int_t       maxFiles       = 200
)
{
    // path for the "AliRsnLoad.C" macro which loads all libraries
    TString strLoadPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
    //TString strLoadPath(".");
    strLoadPath.Append("/AliRsnLoad.C");
    
    // path for the "AliRsnSimpleAnalysisTask.C" macro which executes the job
    TString strTaskPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
    //TString strTaskPath(".");
    strTaskPath.Append("/AliRsnSimpleTask.C");
    
    // path for the configuration macro macro which executes the job
    //TString strConfigPath(Form("%s/PWG2/RESONANCES/macros", getenv("ALICE_ROOT")));
    TString strConfigPath(".");
    strConfigPath += '/';
    strConfigPath.Append(configFile);
    
    // message
    cout << "Path for PAR library loading macro: " << strLoadPath.Data() << endl;
    cout << "Path for configuration macro      : " << strConfigPath.Data() << endl;
    cout << "Path for job execution macro      : " << strTaskPath.Data() << endl;
    
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
    // initialize a counter to check the number of read files
    Int_t nfiles = 0;
    myCollection->Reset();
    while (myCollection->Next()) {
        char fileName[255];
        sprintf(fileName, "%s", myCollection->GetTURL(""));
        Info("AliRsnSimpleTaskAlien", Form("Adding file '%s'", fileName));
        analysisChain->Add(fileName);
        nfiles++;
        if (maxFiles > 0 && nfiles >= maxFiles) break;
    }
    Info("AliRsnSimpleTaskAlien", Form("# processed events = %d", (Int_t)analysisChain->GetEntries()));

    // load read macro
    gROOT->LoadMacro(strTaskPath.Data());
    return AliRsnSimpleTask(analysisChain, strConfigPath.Data(), "output.root", useESDsource);
}
