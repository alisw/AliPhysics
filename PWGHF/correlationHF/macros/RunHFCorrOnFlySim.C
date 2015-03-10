
/*_____________________________________________________________
 
 Run analysis macro for HF Correlations OnFlySim task
 Jitendra Kumar (jitendra.kumar@cern.ch)
 Andrea Rossi   (andrea.rossi@cern.ch)
 _____________________________________________________________*/


void SourceEnv_Libs(){
    
    gSystem->Setenv("GEN_TOTAL_EVENTS" , "50");
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/PYTHIA6 -g");
    
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");
    
    // Analysis-specific
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libEMCALUtils.so");
    gSystem->Load("libJETAN.so");
    
    // MC generator libraries
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia6.4.25.so");
    gSystem->Load("libEGPythia6.so");
    gSystem->Load("libAliPythia6.so");
    gSystem->Load("libgeant321");
    
}


//_____________________| Run Analysis
void RunHFCorrOnFlySim(Long64_t nEvents = 2000){
    
    SourceEnv_Libs();
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/PYTHIA6 -g");
    
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    
    // Seeting for analysis run
    TString      analysisMode = "local"; // "local", "grid", or "proof"
    TString        pluginmode = "test"; // full/ terminate/test mode
    Bool_t     useAlienPlugin =  0;
    
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode);
    if(!alienHandler) return;
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
    mgr->SetDebugLevel(11);
    mgr->SetFileInfoLog("fileinfo.log");
    if(useAlienPlugin)mgr->SetGridHandler(alienHandler);
    
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    
    //const Int_t seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    TString PYTHIA_TUNE = "Perugia2011";
    AliGenerator *genPythia = NULL;
    genPythia = AddPythiaGenSettings(PYTHIA_TUNE);
    mcInputHandler->SetGenerator(genPythia);
     //mcInputHandler->SetSeedMode(2);
    mcInputHandler->SetSeedMode(1);
    mcInputHandler->SetSeed(10);
    
    mgr->SetInputEventHandler(new AliDummyHandler());
    mgr->SetMCtruthEventHandler(mcInputHandler);
    
    gROOT->LoadMacro("AliAnalysisHFCorrOnFlySim.cxx++g");
    gROOT->LoadMacro("AddTaskHFCorrOnFlySim.C");
    AliAnalysisHFCorrOnFlySim *task = AddTaskHFCorrOnFlySim(PYTHIA_TUNE);
    
    if(!mgr->InitAnalysis())return;
    mgr->PrintStatus();
    
    if(!useAlienPlugin)mgr->EventLoop(nEvents);
    if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
    if(analysisMode=="grid")mgr->StartAnalysis(analysisMode.Data());
    return;
    
}

//_____________________| Alien Handler function and plugin
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test")
{
    
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    
    //plugin->SetProductionMode();
    plugin->AddDataFile("dummy.xml"); // TODO to be removed
    plugin->SetRunMode(pluginmode.Data());
    plugin->SetUser("jikumar");
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-08-6");
    plugin->SetAliROOTVersion("vAN-20140529");
    
    plugin->SetAdditionalRootLibs("libVMC.so libPhysics.so libTree.so libMinuit.so libProof.so libSTEERBase.so libESD.so libAOD.so");
    
    plugin->SetDefaultOutputs(kTRUE);
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
    
    
    //Method 1: To create automatically xml through plugin
    plugin->SetAnalysisSource("AliAnalysisHFCorrOnFlySim.cxx");
    plugin->SetAdditionalLibs("AliAnalysisHFCorrOnFlySim.cxx AliAnalysisHFCorrOnFlySim.h liblhapdf.so libpythia6.4.25.so libEGPythia6.so libAliPythia6.so");
    
    plugin->SetMaxMergeFiles(15);
    plugin->SetTTL(70000);
    plugin->SetAnalysisMacro("test.C");
    plugin->SetInputFormat("xml-single");
    plugin->SetValidationScript("validation.sh");
    TString excludeFiles = "";
    plugin->SetRegisterExcludes(excludeFiles + " AliAOD.root");
    TString additionalpackages = "";
    plugin->AddExternalPackage(additionalpackages);
    
    plugin->SetSplitMaxInputFileNumber(100);
    plugin->SetJDLName("testsim.jdl");
    plugin->SetExecutable("testSim.sh");
    plugin->SetSplitMode("se");
    plugin->SetGridWorkingDir("Sim");
    
    plugin->SetMergeViaJDL();
    
    
    //plugin->GenerateTest("Jitendra","");
    Long64_t totalEvents = 5;
    Int_t splitMaxInputFileNumber = 2;
    Int_t nFiles = 0;
    Long64_t neededJobs = totalEvents / splitMaxInputFileNumber;
    plugin->SetMCLoop(true);
    plugin->SetSplitMode(Form("production:1-%d", neededJobs));
    plugin->SetNMCjobs(neededJobs);
    plugin->SetNMCevents(nFiles);
    plugin->SetExecutableCommand("aliroot -b -q");
    plugin->SetKeepLogs(kTRUE);
    plugin->SetNtestFiles(nFiles);
    
    return plugin;
    
}


//__________________________________________| Setting for..Generator..
AliGenerator* AddPythiaGenSettings(TString PythiaTune = "Perugia2011"){
    
    Float_t e_cms = 7000.0;
    Bool_t cr = kTRUE;
    
    AliGenPythia* genP = new AliGenPythia(-1);
    genP->SetMomentumRange(0, 999999.);
    genP->SetThetaRange(0., 180.);
    genP->SetYRange(-20.,20.);
    genP->SetPtRange(0,1000.);
    genP->SetVertexSmear(kPerEvent);
    genP->SetProcess(kPyMbDefault);
    genP->SetEnergyCMS(e_cms);
    
    if(PythiaTune == "Perugia0"){
        genP->SetTune(320);
        if(!cr) genP->SetTune(324);
        cout << "---Perugia0 tunes--" << endl;
    
    }else if(PythiaTune == "Perugia2010"){
        genP->SetTune(327);
        if(!cr) genP->SetTune(324);
        cout << "---Perugia2010 tunes--" << endl;
        
    }else if(PythiaTune == "Perugia2011"){
        genP->SetTune(350);
        if(!cr) genP->SetTune(354);
        cout << "---Perugia2011 tunes--" << endl;
        
    }else if(PythiaTune == "Perugia2012"){
        genP->SetTune(370);
        if(!cr)genP->SetTune(375);
        cout << "---Perugia2012 tunes--" << endl;

    }else cout << "Select Proper tunes" << endl;
    
    
    //genP->UseNewMultipleInteractionsScenario();
        
    //genP->SetCrossingAngle(0,0.000280); //! WARNING for crossing angle seeting (Phi distributions weird)

    //Float_t sigmaz = 6.245; // [cm]
    //Float_t sigmax = 0.0025;
    //Float_t sigmay = 0.0029;
    //genP->SetOrigin(0., 0., 0.); // Taken from OCDB
    //genP->SetSigma(sigmax, sigmay, sigmaz);

    genP->Init();
    genP->Print();
    return genP;
}


