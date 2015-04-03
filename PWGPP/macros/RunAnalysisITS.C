class AliAnalysisGrid;

void RunAnalysisITS(TString pluginmode="",Int_t firstrun=177173,Int_t lastrun=177173,
		    Bool_t readMC=kFALSE,
		    Bool_t runAlign=kFALSE,
		    Bool_t runITS=kFALSE,
		    Bool_t runImpPar=kTRUE,
		    Bool_t runVtx=kFALSE,
		    Bool_t runSPD=kFALSE) 
{
  //
  // Macro to analyze ESDs from raw data reconstruction
  // A.Dainese, andrea.dainese@pd.infn.it
  //

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS -g"); 

  //
  TString analysisMode = "grid"; // "local", "grid", or "proof" (not yet)

  Long64_t nentries=1000000000000000,firstentry=0;
  Bool_t useAlienPlugin=kFALSE;
  Bool_t uselibPWGPP=kTRUE;
  TString loadMacroPath="./";
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
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  if(uselibPWGPP) {gSystem->Load("libTender");gSystem->Load("libPWGPP");}

  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,uselibPWGPP,firstrun,lastrun);
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
  esdH->SetReadFriends(1);
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
    if(!uselibPWGPP) gROOT->LoadMacro("AliAlignmentDataFilterITS.cxx++g");
    taskName="AddTaskAlignmentDataFilterITS.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAlignmentDataFilterITS *alignTask = AddTaskAlignmentDataFilterITS();
  }
  if(runITS) {
    if(!uselibPWGPP) gROOT->LoadMacro("AliAnalysisTaskITSTrackingCheck.cxx++g");
    taskName="AddTaskPerformanceITS.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskITSTrackingCheck *itsTask = AddTaskPerformanceITS(readMC,kFALSE,kFALSE,0,1000000,1);  
  }
  if(runImpPar) {
    if(!uselibPWGPP) gROOT->LoadMacro("AliAnalysisTaskSEImpParRes.cxx++g");
    taskName="AddTaskImpParRes.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEImpParRes *d0Task = AddTaskImpParRes(readMC,-1,kFALSE,kFALSE,0,1000000,0);  
  }
  if(runVtx) {
    if(!uselibPWGPP) gROOT->LoadMacro("AliAnalysisTaskVertexESD.cxx++g");
    taskName="AddTaskVertexESD.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskVertexESD *vtxTask = AddTaskVertexESD(readMC);
  }
  if(runSPD) {
    if(!uselibPWGPP) gROOT->LoadMacro("AliAnalysisTaskSPD.cxx++g");
    taskName="AddTaskSPDQA.C"; 
    taskName.Prepend("$ALICE_PHYSICS/PWGPP/PilotTrain/");
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSPD *spdTask = AddTaskSPDQA();
  }
  /*  
  if(!uselibPWGPP) gROOT->LoadMacro("AliTrackMatchingTPCITSCosmics.cxx++g");
  taskName="AddTaskTrackMatchingTPCITS.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliTrackMatchingTPCITSCosmics *tpcitsTask = AddTaskTrackMatchingTPCITS();
  if(readHLT) tpcitsTask->SetReadHLTESD(kTRUE);  
  */

  
  // Apply the event selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  Bool_t bkgRej=kTRUE;
  //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(readMC,bkgRej);
  

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
				    Bool_t uselibPWGPP=kFALSE,
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
   if(!uselibPWGPP) {
     plugin->SetAnalysisSource("AliAnalysisTaskITSTrackingCheck.cxx AliAlignmentDataFilterITS.cxx AliAnalysisTaskSEImpParRes.cxx AliAnalysisTaskVertexESD.cxx");
     //plugin->SetAnalysisSource("AliAnalysisTaskVertexESD.cxx");
     //plugin->SetAnalysisSource("AliAlignmentDataFilterITS.cxx");
     //plugin->SetAnalysisSource("AliTrackMatchingTPCITSCosmics.cxx");
   }
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   //plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TPC -I$ALICE_PHYSICS/CONTAINERS -I$ALICE_PHYSICS/STEER -I$ALICE_PHYSICS/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS -g");
   if(!uselibPWGPP) {
     //plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
     plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx AliAnalysisTaskITSTrackingCheck.h AliAnalysisTaskITSTrackingCheck.cxx AliAnalysisTaskSEImpParRes.h AliAnalysisTaskSEImpParRes.cxx AliAnalysisTaskVertexESD.h AliAnalysisTaskVertexESD.cxx libGui.so libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
     //plugin->SetAdditionalLibs("AliAnalysisTaskVertexESD.h AliAnalysisTaskVertexESD.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
   } else {
     plugin->SetAdditionalLibs("libGui.so libProof.so libMinuit.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so libTPCbase.so libTPCrec.so libTRDbase.so libTRDrec.so libTender.so libPWGPP.so");
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
