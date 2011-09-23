void runAnalysisTaskITSTPCalignment()
{
  TStopwatch timer;
  timer.Start();

  //runProof("/ITS/dainesea/run104070#esdTree");
  //runProof("/ITS/dainesea/run104070_newTPCalign#esdTree");
  //runProof("/ALIREC/aliprod/run104671#esdTree");
  //runProof("/ITS/dainesea/run104070_newTPCcalib#esdTree");
  //runProof("/COMMON/COMMON/LHC09d9a_0.9TeV_0.5T#esdTree");
  //runLocal("find /data/alice3/mikolaj/runs_pp/pass1/ -name AliESDs.root","tree");
  //runLocal("find /data/alice3/mikolaj/ITSmisal -path */0.001/* -name AliESDs.root","tree");
  //runLocal("find /data/alice3/mikolaj/TPCfullmisalignmentB0 -name AliESDs.root","tree");
  //runLocal("find /data/alice3/mikolaj/TPCfullmisalignmentB5/ -name AliESDs.root -path \"*187824*\" ","tree");
  //runLocal("find /data/alice3/mikolaj/TPCfullmisalignmentB5/ -name AliESDs.root","tree");
  runLocal("fileList");
  //runAlienPlugin("terminate");

  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
void runLocal(TString inputFile = "fileList", TString options="")
{

  TString outputFilename = "outputITSTPCalignment.root";

  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libITSbase");
  gSystem->Load("libITSrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDrec");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libHMPIDbase");
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWG1");

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  //gROOT->LoadMacro("AliRelAlignerKalman.cxx++");
  //gROOT->LoadMacro("AliRelAlignerKalmanArray.cxx++g");
  //gROOT->LoadMacro("AliAnalysisTaskITSTPCalignment.cxx++g");
  gROOT->LoadMacro("AddTaskITSTPCalignment.C");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(inputFile,100000);

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ITSTPCalignmentAnalysisManager");

  // input handlers
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kTRUE);
  esdH->SetInactiveBranches("Calo FMD");
  mgr->SetInputEventHandler(esdH);
  
  //AliMCEventHandler* MCH = new AliMCEventHandler();
  //mgr->SetMCtruthEventHandler(MCH);

  //add the task
  AliAnalysisTaskITSTPCalignment *task = AddTaskITSTPCalignment();

  //start analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("local",chain);
}

//______________________________________________________________________________
void runProof(const char* dataset, TString options="" )
{

  TString outputFilename = "outputITSTPCalignment.root";

  printf("****** Connect to PROOF *******\n");
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("mkrzewic@alicecaf.cern.ch");
  //gProof->SetParallel();
  gProof->ClearPackages();

  // Enable the Analysis Package
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  gProof->UploadPackage("CORRFW.par");
  gProof->EnablePackage("CORRFW");

  //gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-18-12-AN/AF-v4-18-12-AN");
  //gProof->EnablePackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-18-12-AN/AF-v4-18-12-AN");

  gProof->Load("AliRelAlignerKalman.cxx++g");
  gProof->Load("AliRelAlignerKalmanArray.cxx++g");
  gProof->Load("AliAnalysisTaskITSTPCalignment.cxx++g");
  gProof->Load("AddTaskITSTPCalignment.C");

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ITSTPCalignmentAnalysisManager");

  // input handlers
  AliVEventHandler* esdH = new AliESDInputHandler();
  ((AliESDInputHandler*)esdH)->SetReadFriends(kTRUE);
  mgr->SetInputEventHandler(esdH);

  //add the task
  AliAnalysisTaskITSTPCalignment *task = AddTaskITSTPCalignment();

  //start analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("proof",dataset,10000,0);
}

//______________________________________________________________________________
void runAlienPlugin(const char* pluginmode="full")
{
  //must be configured in CreateAlienHandler.C file

  TString outputFilename = "outputITSTPCalignment.root";

  // Load common libraries
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libITSbase");
  gSystem->Load("libITSrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDrec");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libHMPIDbase");
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWG1");

  // Use AliRoot includes to compile our task
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Create and configure the alien handler plugin
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(pluginmode);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-26-00b-4");
  plugin->SetAliROOTVersion("v4-19-11-AN");
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  plugin->SetGridDataDir("/alice/data/2010/LHC10b");
  // Set data search pattern
  plugin->SetDataPattern("*/pass1/*AliESDs.root");
  // ...then add run numbers to be considered
  plugin->SetRunPrefix("000");
  plugin->AddRunNumber(117222);
  plugin->AddRunNumber(117220);
  plugin->AddRunNumber(117112);
  plugin->AddRunNumber(117099);
  plugin->AddRunNumber(117048);
  plugin->AddRunNumber(117048);
  plugin->AddRunNumber(116288);
  plugin->AddRunNumber(115322);
  plugin->AddRunNumber(114931);
  //plugin->AddRunNumber(117051);
  //plugin->SetRunRange(117029,117121);
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //   plugin->AddDataFile("/alice/cern.ch/user/m/mkrzewic/pp2009pass2.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("analysisGridLHC10bPass1_3runs");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliRelAlignerKalman.cxx AliRelAlignerKalmanArray.cxx AliAnalysisTaskITSTPCalignment.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("AliRelAlignerKalman.h AliRelAlignerKalman.cxx AliRelAlignerKalmanArray.h AliRelAlignerKalmanArray.cxx AliAnalysisTaskITSTPCalignment.h AliAnalysisTaskITSTPCalignment.cxx");// AddTaskITSTPCalignment.C");
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  plugin->SetOutputFiles("outputITSTPCalignment.root");
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("AnalysisITSTPCalignmentGenerated.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(20);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(100000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("taskITSTPCalignment.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ITSTPCalignmentAnalysisManager");

  // input handlers
  mgr->SetGridHandler(plugin);

  //input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetInactiveBranches("Calo FMD");
  esdH->SetReadFriends(kTRUE);
  mgr->SetInputEventHandler(esdH);

  gROOT->LoadMacro("AddTaskITSTPCalignment.C");

  //add the task
  AliAnalysisTaskITSTPCalignment *task = AddTaskITSTPCalignment();

  //mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
}

//______________________________________________________________________________
Int_t setupPar(const char* pararchivename)
{
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename)
  {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
    {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh"))
      {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
    {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }

    gSystem->ChangeDirectory("../");
  }
  return 1;
}


