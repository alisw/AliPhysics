void RunGridHLT(const char* pluginmode = "full", Bool_t theMCon=kFALSE, Bool_t UseTree=kFALSE) {
    // Load Libraries.
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -g");
    // load base root libraries
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");
    gSystem->Load("libMinuit");
    // Load analysis framework libraries
    gSystem->Load("libSTEERBase");
    gSystem->Load("libSTEER.so");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libGui.so");
    gSystem->Load("libProof.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libRAWDatabase.so");
    gSystem->Load("libRAWDatarec.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libXMLIO.so");
    gSystem->Load("libXMLParser.so");
    gSystem->Load("libCDB.so");
    gSystem->Load("libSTAT.so");
    gSystem->Load("libTOFbase.so");
    gSystem->Load("libTOFrec.so");
    //gSystem->Load("libTOFcalib.so");
    //__________________________________________________________________________
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    //__________________________________________________________________________
    // Create and configure the alien handler plugin
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode, theMCon);
    if (!alienHandler) return;
    alienHandler->Print();
    //__________________________________________________________________________
    // Create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");
    //mgr->SetDebugLevel(0);
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
    // Input handlers
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(kTRUE);
    mgr->SetInputEventHandler(esdH);
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    // Compile libraries and Analysis Task.
    gROOT->LoadMacro("./AliAnalysisBGMonitorQAHLT.cxx++g");
    // Add Analysis Task.
    gROOT->LoadMacro("./AddTaskBGMonitorQAHLT.C");
    AliAnalysisBGMonitorQAHLT *task = AddTaskBGMonitorQAHLT(UseTree);
    // Proceed with running.
    //mgr->SetDebugLevel(10);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
    if(pluginmode=="local") mgr->StartAnalysis("local",chain);
}
// -------------------------------------------------------------------------
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",
Bool_t theMCon=kFALSE)
{
    // Check if user has a valid token, otherwise make one. This has limitations.
    // One can always follow the standard procedure of calling alien-token-init then
    //   source /tmp/gclient_env_$UID in the current shell.
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(pluginmode.Data());
    plugin->SetUser("blim");
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    //plugin->SetAliROOTVersion("v5-08-09-1");
    plugin->SetAliPhysicsVersion("vAN-20160501-1");
    //plugin->SetAliROOTVersion("v5-06-34");
    //plugin->SetAliPhysicsVersion("v5-06-34-01");
    plugin->SetNtestFiles(1);
    //************************************************
    // Set data search pattern for DATA
    //************************************************
    plugin->SetGridDataDir("/alice/data/2015/LHC15j/"); // specify LHC period
    //225582   plugin->SetDataPattern("pass1/*ESDs.root"); // specify reco pass and AOD set
    plugin->SetDataPattern("muon_calo_pass2/*ESDs.root"); // specify reco pass and AOD set
    plugin->SetRunPrefix("000");   // real data
    Int_t nruns = 0;
    //plugin->AddRunNumber(235226); // data only
    //plugin->AddRunNumber(238468); // data only
    plugin->AddRunNumber(237356); // data only
    nruns++;
    plugin->SetOutputToRunNo(kTRUE);
    plugin->SetNrunsPerMaster(nruns);
    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir("BGMonitorQA_HLT/20160509_02");
    plugin->SetGridOutputDir("output");
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TOF -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include");
    plugin->SetAdditionalLibs("AliAnalysisBGMonitorQAHLT.h AliAnalysisBGMonitorQAHLT.cxx libGui.so libProof.so libMinuit.so libRAWDatabase.so libRAWDatarec.so libANALYSIS.so  libANALYSISalice.so libXMLIO.so libXMLParser.so libCDB.so libSTEERBase.so libSTEER.so libSTAT.so"); //libOADB.so
    plugin->SetAnalysisSource("AliAnalysisBGMonitorQAHLT.cxx");
    plugin->SetDefaultOutputs(kFALSE);
    // merging via
    plugin->SetMergeViaJDL(kTRUE);
    plugin->SetOneStageMerging(kFALSE);
    //plugin->SetMaxMergeStages(5);
    plugin->SetSplitMaxInputFileNumber(20); //20
    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro("TaskAAbg.C");
    plugin->SetExecutable("TaskAAbg.sh");
    plugin->SetValidationScript("TaskAAbg_validation.sh");
    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName("TaskAAbg.jdl");
    plugin->SetOverwriteMode(kFALSE);
    plugin->SetTTL(72000);
    //Optionally modify job price (default 1)
    //plugin->SetPrice(1);
    //Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");
    // plugin->SetOutputFiles("myTree.root, histos.root");
    plugin->SetOutputFiles("AnalysisResults.root");
    return plugin;
}
TChain *CreateChain(const char *fileName)
{
    TString *treename;
    treename = new TString("esdTree");
    TChain* chainGood = new TChain(treename->Data());
    TChain* chainNull = new TChain(treename->Data());
    int counter=0;
    char line[500] = new char*;
    TString *name;
    ifstream inputstream(fileName);
    if(!inputstream) {cout<<"input file not found"<<endl; return chainNull;}
    while(!inputstream.eof()) {
        inputstream.getline(line,500);
        name = new TString(line);
        name->Append("/AliESDs.root");
        if(!inputstream.eof()) chainGood->Add(name->Data());
        counter++;
    }
    inputstream.close();
    return chainGood;
}