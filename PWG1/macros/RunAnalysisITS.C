class AliAnalysisGrid;

void RunAnalysisITS(TString pluginmode,Int_t firstrun,Int_t lastrun,
		    Bool_t readMC=kFALSE,
		    Bool_t runAlign=kTRUE,
		    Bool_t runITS=kTRUE,
		    Bool_t runImpPar=kTRUE,
		    Bool_t runVtx=kFALSE,
		    Bool_t runSPD=kFALSE) 
{
  //
  // Macro to analyze ESDs from raw data reconstruction
  // A.Dainese, andrea.dainese@pd.infn.it
  //

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -g"); 

  //
  TString analysisMode = "grid"; // "local", "grid", or "proof" (not yet)

  Long64_t nentries=1000000000000000,firstentry=0;
  Bool_t useAlienPlugin=kTRUE;
  Bool_t uselibPWG1=kFALSE;
  TString loadMacroPath="../../";
  Bool_t readHLT=kFALSE;
  //

  if(analysisMode=="grid") {
    // Connect to AliEn
    TGrid::Connect("alien://");
  } else if(analysisMode=="proof") {
    // Connect to the PROOF cluster
    printf("PROOF mode not yet functional..\n");
    return;
    TProof::Open("alicecaf");
    //TProof::Reset("alicecaf");
  }

  // Load analysis libraries
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  if(uselibPWG1) {gSystem->Load("libTENDER.so");gSystem->Load("libPWG1.so");}

  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,uselibPWG1,firstrun,lastrun);
    if(!alienHandler) return;
  }

  TChain *chainESD = 0;
  if(!useAlienPlugin) {
    // Prepare input chain
    //    chainESD = CreateESDChain("/home/dainesea/alignData/RAWdata_CosmicsSum09/RecoSPD/chunk.",13,13);
    chainESD=new TChain("esdTree");
    //chainESD->Add("alien:///alice/cern.ch/user/s/sitta/output/000088361/ESDs/pass1/09000088361017.10/AliESDs.root");
    chainESD->Add("AliESDs.root");
  }

  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  // Enable debug printouts
  mgr->SetDebugLevel(10);
  // Connect plug-in to the analysis manager
  if(useAlienPlugin) mgr->SetGridHandler(alienHandler);

  // Add ESD handler
  Bool_t readRP=kFALSE;
  AliESDInputHandler *esdH = 0;
  if(readRP) {
    esdH = new AliESDInputHandlerRP();
  } else {
    esdH = new AliESDInputHandler();
  }
  esdH->SetActiveBranches("ESDfriend");
  if(readHLT) esdH->SetReadHLT();
  mgr->SetInputEventHandler(esdH);
  if(readMC) {
    AliMCEventHandler  *mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH); 
  }
  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  TString taskName;
  
  if(runAlign) {
    if(!uselibPWG1) gROOT->LoadMacro("AliAlignmentDataFilterITS.cxx++g");
    taskName="AddTaskAlignmentDataFilterITS.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAlignmentDataFilterITS *alignTask = AddTaskAlignmentDataFilterITS();
  }
  if(runITS) {
    if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskITSTrackingCheck.cxx++g");
    taskName="AddTaskPerformanceITS.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskITSTrackingCheck *itsTask = AddTaskPerformanceITS(readMC,kFALSE,kFALSE);  
  }
  if(runImpPar) {
    if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskSEImpParRes.cxx++g");
    taskName="AddTaskImpParRes.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEImpParRes *d0Task = AddTaskImpParRes(readMC,-1,kFALSE);  
  }
  if(runVtx) {
    if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskVertexESD.cxx++g");
    taskName="AddTaskVertexESD.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskVertexESD *vtxTask = AddTaskVertexESD(readMC);
  }
  if(runSPD) {
    if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskSPD.cxx++g");
    taskName="AddTaskSPDQA.C"; 
    taskName.Prepend("$ALICE_ROOT/PWG1/PilotTrain/");
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSPD *spdTask = AddTaskSPDQA();
  }
  /*  
  if(!uselibPWG1) gROOT->LoadMacro("AliTrackMatchingTPCITSCosmics.cxx++g");
  taskName="AddTaskTrackMatchingTPCITS.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliTrackMatchingTPCITSCosmics *tpcitsTask = AddTaskTrackMatchingTPCITS();
  if(readHLT) tpcitsTask->SetReadHLTESD(kTRUE);  
  */

  
  // Apply the event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  Bool_t bkgRej=kTRUE;
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(readMC,bkgRej);
  

  //
  // Run the analysis
  //    
  if(chainESD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainESD->GetEntries());

  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
  mgr->StartAnalysis(analysisMode.Data(),chainESD,nentries,firstentry);

  return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",
				    Bool_t uselibPWG1=kFALSE,
				    Int_t firstrun,Int_t lastrun)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(pluginmode.Data());
   plugin->SetUser("dainesea");
   plugin->SetNtestFiles(1);
   // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-26-00b-6");
   plugin->SetAliROOTVersion("v4-19-15-AN");
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   TString wdname="analysisITS_Runs_";
   wdname+=firstrun;
   wdname.Append("_");
   wdname+=lastrun;
   plugin->SetGridWorkingDir(wdname.Data());
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
   // Declare input data to be processed.
   // Method 1: Create automatically XML collections using alien 'find' command.
   // Define production directory LFN
   plugin->SetGridDataDir("/alice/data/2010/LHC10b");
   // Set data search pattern
   plugin->SetDataPattern("pass2/*AliESDs.root");
   //plugin->SetDataPattern("ESD.tag.root");
   Int_t n=0;
   FILE *in = fopen("/home/dainesea/alignData/RAWdata_pp10/goodruns_pp10.txt","r");
   if(!in) printf("run file not found\n");
   Int_t lines=0; 
   Float_t runnumber; 
   while(1) {
     Int_t ncol = fscanf(in,"%f",&runnumber);
     if(ncol<1) break;
     if(runnumber<firstrun || runnumber>lastrun) continue;
     TString runnumberstring="000";
     Int_t runnumberint=(Int_t)runnumber;
     runnumberstring+=runnumberint;
     n++; plugin->AddRunNumber(runnumberstring.Data());
   }
   fclose(in);
   plugin->SetNrunsPerMaster(n);
   // Method 2: Declare existing data files (raw collections, xml collections, root file)
   // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
   // XML collections added via this method can be combined with the first method if
   // the content is compatible (using or not tags)
   // e.g.: find -z -x 80015 /alice/data/2009/LHC09c/000080015/ESDs/ ESD.tag.root > 80015.xml
   //plugin->AddDataFile("79876.xml");
   //plugin->AddDataFile("80015.xml");
   //plugin->AddDataFile("80261.xml");
   //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   if(!uselibPWG1) {
     plugin->SetAnalysisSource("AliAnalysisTaskITSTrackingCheck.cxx AliAlignmentDataFilterITS.cxx AliAnalysisTaskSEImpParRes.cxx AliAnalysisTaskVertexESD.cxx");
     //plugin->SetAnalysisSource("AliAnalysisTaskVertexESD.cxx");
     //plugin->SetAnalysisSource("AliAlignmentDataFilterITS.cxx");
     //plugin->SetAnalysisSource("AliTrackMatchingTPCITSCosmics.cxx");
   }
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   //plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -g");
   if(!uselibPWG1) {
     //plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
     plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx AliAnalysisTaskITSTrackingCheck.h AliAnalysisTaskITSTrackingCheck.cxx AliAnalysisTaskSEImpParRes.h AliAnalysisTaskSEImpParRes.cxx AliAnalysisTaskVertexESD.h AliAnalysisTaskVertexESD.cxx libGui.so libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
     //plugin->SetAdditionalLibs("AliAnalysisTaskVertexESD.h AliAnalysisTaskVertexESD.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
   } else {
     plugin->SetAdditionalLibs("libGui.so libProof.so libMinuit.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so libTPCbase.so libTPCrec.so libTRDbase.so libTRDrec.so libTENDER.so libPWG1.so");
   }
   // Declare the output file names separated by blancs.
   // (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs(kTRUE);
   // Optionally define the files to be archived.
   //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
   // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("AnalysisITS.C");
   plugin->SetExecutable("analysisITS.sh");
   plugin->SetExecutableCommand("root.exe -b -q");
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(5);
   // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   //plugin->SetMaxInitFailed(5);
   // Optionally resubmit threshold.
   //plugin->SetMasterResubmitThreshold(90);
   // Optionally set time to live (default 30000 sec)
   plugin->SetTTL(80000);
   // Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
   // Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("TaskAnalysisITS.jdl");
   // Optionally modify job price (default 1)
   //plugin->SetPrice(1);      
   // Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   // Optionally set the preferred SE    
   plugin->SetPreferedSE("ALICE::CNAF::SE");
   
   return plugin;
}
//-----------------------------------------------------------------------------
TChain *CreateESDChain(TString esdpath=".",Int_t ifirst=-1,Int_t ilast=-1) {


  TChain *chainESD = new TChain("esdTree");

  if(ifirst<0) {
    chainESD->Add("AliESDs.root");
  } else {
    for(Int_t i=ifirst; i<=ilast; i++) {
      TString esdfile=esdpath; esdfile+=i; esdfile.Append("/AliESDs.root");
      chainESD->Add(esdfile.Data());
    }
  }
  
  return chainESD;
}
