class AliAnalysisGrid;

void RunAnalysisITS() {
  //
  // Macro to analyze ESDs from raw data reconstruction
  // A.Dainese, andrea.dainese@pd.infn.it
  //
  //gSystem->Setenv("alien_CLOSE_SE","ALICE::CNAF::SE");

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -g"); 

  //
  TString analysisMode = "local"; // "local", "grid", or "proof" (not yet)

  Long64_t nentries=1000000000000000,firstentry=0;
  Bool_t useAlienPlugin=kFALSE;
  Bool_t uselibPWG1=kTRUE;
  TString pluginmode="full";
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
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  if(uselibPWG1) {gSystem->Load("libTENDER.so"); gSystem->Load("libPWG1.so");}

  // Create Alien plugin, if requested
  if(useAlienPlugin) {  
    AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,uselibPWG1);
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
  AliESDInputHandler *esdH = new AliESDInputHandler();
  if(readHLT) esdH->SetReadHLT();
  mgr->SetInputEventHandler(esdH);

  //-------------------------------------------------------------------

  
  //-------------------------------------------------------------------
  // Analysis tasks (wagons of the train)   
  //
  TString taskName;
  
  if(!uselibPWG1) gROOT->LoadMacro("AliAlignmentDataFilterITS.cxx++g");
  taskName="AddTaskAlignmentDataFilterITS.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAlignmentDataFilterITS *itsTask = AddTaskAlignmentDataFilterITS();
    
  if(!uselibPWG1) gROOT->LoadMacro("AliTrackMatchingTPCITSCosmics.cxx++g");
  taskName="AddTaskTrackMatchingTPCITS.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliTrackMatchingTPCITSCosmics *tpcitsTask = AddTaskTrackMatchingTPCITS();
  if(readHLT) tpcitsTask->SetReadHLTESD(kTRUE);  
  /*
  Bool_t readMC=kTRUE;

  if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskVertexESD.cxx++g");
  taskName="AddTaskVertexESD.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskVertexESD *vtxTask = AddTaskVertexESD(readMC);
    
  if(!uselibPWG1) gROOT->LoadMacro("AliAnalysisTaskITSTrackingCheck.cxx++g");
  taskName="AddTaskPerformanceITS.C"; 
  taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskITSTrackingCheck *itsTask = AddTaskPerformanceITS(readMC);

  if(readMC) {
    AliMCEventHandler  *mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH); 
  }
  */
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
				    Bool_t uselibPWG1=kFALSE)
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
   plugin->SetAPIVersion("V2.4");
   plugin->SetROOTVersion("v5-24-00");
   plugin->SetAliROOTVersion("v4-18-07-AN");
   // Declare input data to be processed.
   // Method 1: Create automatically XML collections using alien 'find' command.
   // Define production directory LFN
   plugin->SetGridDataDir("/alice/data/2009/LHC09c");
   //plugin->SetGridDataDir("/alice/cern.ch/user/s/sitta/output/000088361/");
   // Set data search pattern
   //plugin->SetDataPattern("AliESDs.root");
   plugin->SetDataPattern("ESD.tag.root");
   Int_t n=0;
   n++; plugin->AddRunNumber("000080015");
   n++; plugin->AddRunNumber("000080261");
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
   // Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("analysisITS");
   // Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output151009"); // In this case will be $HOME/work/output
   // Declare the analysis source files names separated by blancs. To be compiled runtime
   // using ACLiC on the worker nodes.
   if(!uselibPWG1) {
     plugin->SetAnalysisSource("AliAlignmentDataFilterITS.cxx");
     plugin->SetAnalysisSource("AliTrackMatchingTPCITSCosmics.cxx");
   }
   // Declare all libraries (other than the default ones for the framework. These will be
   // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   //plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -g");
   if(!uselibPWG1) {
     plugin->SetAdditionalLibs("AliAlignmentDataFilterITS.h AliAlignmentDataFilterITS.cxx AliTrackMatchingTPCITSCosmics.h AliTrackMatchingTPCITSCosmics.cxx libProof.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libITSbase.so libITSrec.so");
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
   // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(1);
   // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   //plugin->SetMaxInitFailed(5);
   // Optionally resubmit threshold.
   //plugin->SetMasterResubmitThreshold(90);
   // Optionally set time to live (default 30000 sec)
   //plugin->SetTTL(20000);
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
