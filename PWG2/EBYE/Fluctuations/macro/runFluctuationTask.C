/****************************************************************

Lunch the Analysis Task for Fluctuation Study 
Author: Satyajit Jena
Email:  sjena@cern.ch

Task : MF, MFT, CF, OF
Running in Local = analmode = local, proof, grid 

*******************************************************************/


TString tasktype = "MF"; //OF, MF, MFT, CF

//CAF-AAF Related
TString dataset      = "/alice/data/LHC10h_000137162_p1_3plus";
TString proofCluster = "alice-caf.cern.ch";
TString alirootVer   = "VO_ALICE@AliRoot::v4-21-12-AN";


Bool_t UsePar        = kFALSE; // Do you want to use par
Bool_t UseTrigger    = kFALSE; // activate trigger seletction
Bool_t UseCentrality = kFALSE; // activate centrality slection


void runFluctuationTask(TString analmode = "local")
{
 
  if(analmode.CompareTo("local") == 0 ) {
    if(!LoadLibraries(analmode)) {
      printf("Library Not loaded\n");
      return;
    }
    runLocal("file.txt",1000);
    
 }
  else  if(analmode.CompareTo("proof") == 0) {
    if(!LoadLibraries(analmode)) return;
    runproof(dataset);
 }
  else printf("load error\n");
}

//___________________________________________________________________________________
Bool_t LoadLibraries(TString mode = "local")
{
  
  if(mode.CompareTo("local") == 0) {
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    //  gSystem->Load("libPWG2ebye.so");
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    printf("Library is Loaded \n");
    return kTRUE;
  }
  else if(mode.CompareTo("proof") == 0) {
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    TProof *p = TProof::Open(proofCluster.Data());
    p->ShowPackages();
    p->EnablePackage(alirootVer.Data());
    return kTRUE;
  }
  else if(mode.CompareTo("proof") == 0) {
    printf("Satya FIXME \n");
  }
  else {
    printf(" ERROR: Give proper running mode \n");
    return kFALSE;
  }
  
}

//______________________________________________________________________
Bool_t LoadIfParAnalysis(TString mode = "local") {
  // FIXME: For par file upload methods..
  


}

//______________________________________________________________________
// Running local
void runLocal(TString input,int nfile)
{

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(input, nfile);
  
  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);


  //  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  //  gROOT->LoadMacro("AddTaskCentrality.C");
  
  //  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();
  //  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();


  if(!UsePar) gROOT->LoadMacro("../AliEbyEEventBase.cxx++g");  
  
  if(tasktype.CompareTo("CF") == 0) {
    if(!par) {
      gROOT->LoadMacro("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
      gROOT->LoadMacro("../AliEbyEChargeFluctuationAnalysis.cxx++g");  
    }
   
    gROOT->LoadMacro("AddFluctuationTask.C");
    AddTaskCF();
  }
  else if(tasktype.CompareTo("MF") == 0) {
    if(!UsePar) {
      gROOT->LoadMacro("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
      gROOT->LoadMacro("../AliEbyEMFAnalysisTask.cxx++g"); 
    }
    gROOT->LoadMacro("AddFluctuationTask.C");
    AddTaskMF();
  }
  else if(tasktype.CompareTo("OF") == 0) {
    if(!UsePar) {
      gROOT->LoadMacro("../AliEbyEFluctuationAnalysis.cxx++g"); 
      gROOT->LoadMacro("../AliEbyEFluctuationTask.cxx++g"); 
    }
    gROOT->LoadMacro("AddFluctuationTask.C");
    AddTaskF();
  }
  else if (tasktype.CompareTo("MFT") == 0) {
    if(!UsePar) {
      gROOT->LoadMacro("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
      gROOT->LoadMacro("../AliEbyEMFAnalysisTaskT.cxx++g"); 
    }
    gROOT->LoadMacro("AddFluctuationTask.C");
    AddTaskMFT();
  }
  else return;

   
  

  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain);


}


//________________________________________________________
//Running in proof

void runproof(TString dataset)
{

     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) mgr = new AliAnalysisManager("FirstCheck");
     
     AliVEventHandler* esdH = new AliESDInputHandler();
     mgr->SetInputEventHandler(esdH);
     
     if(UseTrigger) {
       gProof->Load("AddTaskPhysicsSelection.C");
       AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection();
     }

     if(UseCentrality) {
       gProof->Load("AddTaskCentrality.C");    
       AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
     }

     if(!UsePar) gProof->Load("../AliEbyEEventBase.cxx++g");  
     
     if(tasktype.CompareTo("CF") == 0) {
       if(!UsePar) {
	 gProof->Load("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
	 gProof->Load("../AliEbyEChargeFluctuationAnalysis.cxx++g"); 
       } 
       gProof->Load("AddFluctuationTask.C");
       AddTaskCF();
     }
     else if(tasktype.CompareTo("MF") == 0) {
       if(!UsePar) {
	 gProof->Load("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
	 gProof->Load("../AliEbyEMFAnalysisTask.cxx++g"); 
       }
       gProof->Load("AddFluctuationTask.C");
       AddTaskMF();
     }

     else if(tasktype.CompareTo("OF") == 0) {
       if(!UsePar) {
	 gProof->Load("../AliEbyEFluctuationAnalysis.cxx++g"); 
	 gProof->Load("../AliEbyEAnalysisTask.cxx++g"); 
       }
       gProof->Load("AddFluctuationTask.C");
       AddTaskF();
     }
     
     else if (tasktype.CompareTo("MFT") == 0) {
       if(!UsePar) {
	 gProof->Load("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
	 gProof->Load("../AliEbyEMFAnalysisTaskT.cxx++g"); 
       }
       gProof->Load("AddFluctuationTask.C");
       AddTaskMFT();
     }
     else return;
          
     mgr->SetDebugLevel(0);
     if (!mgr->InitAnalysis()) 
       return;
     
     mgr->PrintStatus();
        
     mgr->StartAnalysis("proof",dataset.Data());
  
}


//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, 
				    const char *proofcluster, const char *proofdataset)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-27-06b");
    plugin->SetAliROOTVersion("v4-21-08-AN");

    // Declare input data to be processed.

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/data/2010/LHC10b");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
    plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetRunPrefix("000");   // real data, comment out for MC
    // ...then add run numbers to be considered
    plugin->AddRunNumber(115514);
    //plugin->SetRunRange(114917,115322);
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    plugin->SetMergeViaJDL();

    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(taskname);

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliAnalysisTaskEx01.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskEx01.h AliAnalysisTaskEx01.cxx");

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    plugin->SetDefaultOutputs();
    //plugin->SetOutputFiles("list.root");

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(Form("%s.C",taskname));

    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(100);

    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s.sh",taskname));

    // set number of test files to use in "test" mode
    plugin->SetNtestFiles(10);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);      

    // Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");
    
    //----------------------------------------------------------
    //---      PROOF MODE SPECIFIC SETTINGS         ------------
    //---------------------------------------------------------- 
    // Proof cluster
    plugin->SetProofCluster(proofcluster);
    // Dataset to be used   
    plugin->SetProofDataSet(proofdataset);
    // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
    plugin->SetProofReset(0);
    // May limit number of workers
    plugin->SetNproofWorkers(0);
    // May limit the number of workers per slave
    plugin->SetNproofWorkersPerSlave(1);   
    // May use a specific version of root installed in proof
    plugin->SetRootVersionForProof("current");
    // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
    plugin->SetAliRootMode("default"); // Loads AF libs by default
    // May request ClearPackages (individual ClearPackage not supported)
    plugin->SetClearPackages(kFALSE);
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode("file1.text"); // file should contain path name to a local directory containg *ESDs.root etc
    // Request connection to alien upon connection to grid
    plugin->SetProofConnectGrid(kFALSE);

    return plugin;
}

