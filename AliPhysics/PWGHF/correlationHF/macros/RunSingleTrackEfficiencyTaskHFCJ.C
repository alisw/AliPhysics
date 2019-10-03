

class AliAnalysisGrid;
class AliAnalysisAlien;

//Modified on: Feb16, 2016
//By: Jitendra



//Parameter to be Set
const Bool_t readAOD = 1;
TString fUsername = "jikumar";
TString  fAnalysisMode    =  "local";
TString    fPluginMode    =  "full";
TString     fInputMode    =  "list";
Long64_t      nEntries    =  123567890, firstentry=0;
Bool_t        IsPlugin    =  kTRUE;
TString TestFilesWPlugin  =  "filesAOD.txt";
Bool_t ifTaskPIDQA = kFALSE;


void RunSingleTrackEfficiencyTaskHFCJ()
{
    
    TBenchmark fBenchMark;
    fBenchMark.Start("AliCFSingleTrackEfficiencyTask");
    
    Load();
    
    if(fAnalysisMode=="grid"){
        //gSystem->Exec(Form("alien-token-init %s",fUsername.Data()));
        TGrid::Connect("alien://") ;
    }
    
    if(IsPlugin) {
        AliAnalysisGrid *alienHandler = CreateAlienHandler();
        if(!alienHandler) return;
    }
    
    printf("CREATE ANALYSIS MANAGER\n");
    AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");
    mgr->SetDebugLevel(10);
    if(IsPlugin) mgr->SetGridHandler(alienHandler);
    
    AliMCEventHandler*  mcHandler = new AliMCEventHandler();
    if (!readAOD) mgr->SetMCtruthEventHandler(mcHandler);
    
    AliInputEventHandler* dataHandler;
    if   (readAOD) dataHandler = new AliAODInputHandler();
    else           dataHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(dataHandler);
    
    if (!readAOD) {
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE);
    }
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE, kTRUE, kTRUE, 2, kFALSE, "", kTRUE, kFALSE, -1);
    
    if(ifTaskPIDQA){
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
        AliAnalysisTaskPIDqa *pidQA = AddTaskPIDqa();
    }
    
    printf("Prepare to create the task\n");
    gROOT->LoadMacro("./AddSingleTrackEfficiencyTaskDhCorrelations.C");
    
    //NCharge Filterbin0 + kFAST configuration + w/o external cut file
    AliCFSingleTrackEfficiencyTask *taskNchF0wFile  = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "NchFbit0woFile", AliPID::kPion, 0, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID, "", "");
    
    //NCharge Filterbin0 + kFAST configuration + w external cut file
    AliCFSingleTrackEfficiencyTask *taskNchF0  = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "NchFbit0", AliPID::kPion, 0, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID, "AssocPartCuts_Std_NewPoolsAndCode_10000Tr.root", "AssociatedTrkCuts");
    
    //Pion Filterbin0+ kFAST configuration + w external cut file
    AliCFSingleTrackEfficiencyTask *taskPionF0 = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "PionFbit0", AliPID::kPion, 211, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID,"AssocPartCuts_Std_NewPoolsAndCode_10000Tr.root", "AssociatedTrkCuts");
    
    //Kaons Filterbin0+ kFAST configuration + w external cut file
    AliCFSingleTrackEfficiencyTask *taskKaonF0 = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "KaonFbit0", AliPID::kKaon, 321, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID,"AssocPartCuts_Std_NewPoolsAndCode_10000Tr.root", "AssociatedTrkCuts");
    
    //Protons Filterbin0+ kFAST configuration + w external cut file
    AliCFSingleTrackEfficiencyTask *taskProtF0 = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "ProtonFbit0", AliPID::kProton, 2212, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID,"AssocPartCuts_Std_NewPoolsAndCode_10000Tr.root", "AssociatedTrkCuts");
    
    //Electron Filterbin0+ kFAST configuration + w external cut file
    AliCFSingleTrackEfficiencyTask *taskElecF0 = AddSingleTrackEfficiencyTaskDhCorrelations(kTRUE, "ElectronFbit0", AliPID::kElectron,11, AliVEvent::kAnyINT, kFALSE, AliCFSingleTrackEfficiencyTask::kFast, AliSingleTrackEffCuts::kNoBayesianPID, "AssocPartCuts_Std_NewPoolsAndCode_10000Tr.root", "AssociatedTrkCuts");
    
    
    
    // Run the analysis
    TChain * analysisChain=0;
    if(analysisChain) printf("CHAIN HAS %d ENTRIES\n",(Int_t)analysisChain->GetEntries());
    if(!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    if(fAnalysisMode=="grid" && !IsPlugin) fAnalysisMode="local";
    if(fAnalysisMode!="proof") {
        mgr->StartAnalysis(fAnalysisMode.Data(),analysisChain,nEntries,firstentry);
    }
    
    fBenchMark.Stop("AliCFSingleTrackEfficiencyTask");
    fBenchMark.Show("AliCFSingleTrackEfficiencyTask");
    
    return;
    
}


//_______________________________| CreateAlienHandler |________________________________
AliAnalysisGrid* CreateAlienHandler()
{
    
    AliAnalysisAlien *Plugin = new AliAnalysisAlien();
    Plugin->SetRunMode(fPluginMode.Data());
    Plugin->SetUser(fUsername.Data());
    Plugin->SetAPIVersion("V1.1x");
    Plugin->SetROOTVersion("v5-34-30-alice-12");
    Plugin->SetAliROOTVersion("v5-07-20-4");
    Plugin->SetAliPhysicsVersion("vAN-20160129-1");
    Plugin->SetNtestFiles(1);
    Plugin->SetFileForTestMode(TestFilesWPlugin.Data());
    
    // Set data search pattern for DATA grid Mode
    Plugin->SetGridDataDir("/alice/sim/2013/LHC13b2_efix_p1"); // specify MC sample
    if(readAOD) Plugin->SetDataPattern("AOD/*AliAOD.root");// specify AOD set
    else Plugin->SetDataPattern("*/AliESDs.root");
    
    Int_t totruns=0;
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/AddGoodRuns.C");
    gROOT->LoadMacro("AddGoodRuns.C");
    totruns += AddGoodRuns(Plugin,"LHC13b","LHC13b2_efix_p1"); //Set accordingly
    Plugin->SetNrunsPerMaster(totruns);
    
    Plugin->SetGridWorkingDir("SingleTrkEff/LHC13b2_efix_p1/16Feb2016/");
    Plugin->SetGridOutputDir("out");
    
    Plugin->SetExecutable("STE16Feb2016.sh");
    Plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER/STEER -I$ALICE_PHYSICS/STEER/STEERBase -I$ALICE_PHYSICS/STEER/ESD -I$ALICE_PHYSICS/STEER/AOD -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS -I$ALICE_PHYSICS/OADB -g");
    
    Plugin->SetSplitMaxInputFileNumber(5);
    Plugin->SetDefaultOutputs(kTRUE);
    Plugin->SetMergeViaJDL(kTRUE); // Merge Via JDL
    Plugin->SetOneStageMerging(kFALSE);
    Plugin->SetMaxMergeStages(2);
    Plugin->SetAnalysisMacro("STE16Feb2016.C");
    Plugin->SetJDLName("STE16Feb2016.jdl");
    
    return Plugin;
}


//_____| Loading Libraries
void Load() {
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER/STEER -I$ALICE_PHYSICS/STEER/STEERBase -I$ALICE_PHYSICS/STEER/ESD -I$ALICE_PHYSICS/STEER/AOD -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGPP -g");
    
    //load the required aliroot libraries
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libMinuit");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGPP");
    gSystem->Load("libPWGTools");
}

