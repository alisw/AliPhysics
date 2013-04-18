//-*- Mode: C++ -*-
// $Id$

/// @file   run-single-task.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-04-06
/// @brief  Run a single task
///
/// Helper macro to run a single task either locally or on Grid
/// Usage:
/// aliroot -b -q -l run-single-task.C'("mode", "input", "tasks", "name", "options", events, "path", "pattern", "friendPattern", "outputDir", "user")'
///  arguments
///   mode:    local, full, test, compile, merge, collect
///   input:   In grid mode a list of run numbers. 
///            In local mode:
///             - ESD file AliESDs.root
///             - a list of ESD files in a text file <choose-file-name>.txt
///             - AOD file AliAOD.root
///             - AODs with predefined list/AODs in same folder, specifiy as "AOD"
///
///   tasks:   list of class names/source or header files of tasks, or AddTask-macros
///
///  optional arguments
///   name:    analysis name (default 'myanalysis')
///   options: optional arguments passed onto the task, e.g. 'mc', 'event-mixing'
///            options of run-single-task:
///              'mcData' -> the run numbers indicate MC Data, no '000' prepended
///   events:  number of events to be processed (default -1 -> all)
///   path:    data search path for grid analysis (default from configuration file)
///   pattern: data search pattern (default from configuration file)
///   friend pattern: friend file search pattern (default from configuration file)
///   output dir: output directory in Grid home (default gridwork/yyyy-mm-dd_hh-mm)
///   user:    default NULL, using user of active token
///
/// Examples:
/// aliroot -b -q -l run-single-task.C'("full", "146860", "AliAnalysisTaskSample", "myanalysis_LHC11a")'
///
/// aliroot -b -q -l run-single-task.C'("local", "$ALICE_ROOT/test/ppbench/AliESDs.root", "AliAnalysisTaskSample")'
///
/// aliroot -b -q -l run-single-task.C'("local", "AOD", "AddTaskSample.C")'
///
/// aliroot -b -q -l run-single-task.C'("full", "146860", "AliAnalysisTaskSample", "correlation3p_LHC11a", 0, -1, "/alice/data/2011/LHC11a", "*/pass2_without_SDD/AOD*/*/AliAOD.root")'
///
/// aliroot -b -q -l run-single-task.C'("merge", "gridwork/mydir", "AliAnalysisTaskSample", "myanalysis_LHC11a")'
///
/// aliroot -b -q -l run-single-task.C'("collect", "gridwork/mydir", "AliAnalysisTaskSample", "myanalysis_LHC11a")'
///
/// Data input:
/// depending on the format of the search pattern either the ESD or AOD input handler is used.
///
/// Source files:
/// If the task and classes used by the task are not in an AliRoot library available, e.g.
/// for the purpose of development, all header and source files need to be in the local
/// directory. The macro searches automatically for dependencies, compiles those and
/// copies files to the Grid working directory. In order to make the files accessible in
/// the local directory, the files can be just linked.
/// <pre>
/// for f in <search path>; do ln -s $f; done
/// </pre>
/// If there are dependencies (include headers) which are not available in the working
/// directory, the macro searches for the aliroot libraries implementing them, and adds
/// the libraries if found to the setup.
///
/// Local analysis:
/// requires only the path to the input file and the task class name. If the specified file is
/// a text file (.txt) each line can contain an input ESD file path, all files are chained.
/// Analysis on local AOD files needs to be setup prior to this macro. gDirectory must contain
/// a TChain object of name 'aodTree'. This is for example created by macros like
/// $ALICE_ROOT/PWGHF/vertexingHF/MakeAODInputChain.C
/// Set $ALICE_ROOT/PWGHF/correlationHF/macros/setupDxHFE.C for an example.
///
/// Grid analysis:
/// All modes provided by the AliAnalysisAlien plugin can be used, e.g. full, test, offline
/// A couple of settings need to be defined in a configuration file 'grid-config.C' which can be
/// either in the local directory or home directory. The file can look like
/// <pre>
///   const char* alienAPIVersion="V1.1x";
///   const char* alienROOTVersion="v5-34-01";
///   const char* alienAliROOTVersion="v5-03-61-AN";
///   const char* defaultGridDataDir="/alice/data/2011/LHC11a";
///   const char* defaultDataPattern="*/pass2_without_SDD/*/AliESDs.root";
///   const char* defaultFriendDataPattern="";
///   {} // note this empty body
/// </pre>
/// Data path and pattern can also be specified as command line arguments.
/// The working directory in the grid home directory of the user is set to
/// gridwork/<date>_<time> (gridwork/yyyy-mm-dd_hh-mm), can be overridden by command line
/// parameter.
///
/// Options:
/// Options to the macro can be propagated via the parameter 'arguments', the known options
/// are parsed and filtered out from the arguments, which are than further propagated to
/// AdTask macros and task.
/// --mcData        switch indicates that the input data is mc data, the run numbers have
///                 a different format in real data
/// --nTestFiles=   number of test files to be used in test mode (default 10)
/// --merge=        merging mode 'local', 'grid', 'collect' (default Grid)
///
/// Merging of output:
/// The result of Grid analysis is merged remotely by using the ability of AliAnalysisAlien
/// to create stages of merging jobs. Being in the directory where analysis was originally
/// launched, the macro can be run in mode 'merge' with the remote working directory,
/// tasknames and the analysis name as arguments. Final result can be collected in mode
/// 'collect' with the same arguments, see examples.
/// Optionally, within the 'options' argument a list of output directories in the
/// remote grid directory can be specified, e.g. "000 002".
///
/// Merging of output in mode 'terminate':
/// The output files can be merged locally by using the argument '--merge=local'. In that
/// case all files are downloaded to the local machine and merged. The merging on grid
/// requires to launch the analysis with argument '--merge=grid'. After the jobs are done
/// further steps are required, 
/// 1) run in mode "terminate" with argument '--merge=grid' and working directory on grid,
/// 2) run in mode "terminate" with argument '--merge=collect' and working directory on grid.
/// Step 1) can be repeated  multiple times, the AnalysisManager will in each stage merge
/// several files of the previous stage, it will notify you when the final result has
/// already been merged.
/// 
/// Suggestions:
/// Feedback appreciated: Matthias.Richter@cern.ch
/// If you find this macro useful but not applicable without patching it, let me know
/// your changes and use cases.


#if defined(__CINT__) && !defined(__MAKECINT__)
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment
//
int macroVerbosity=0;
const char* defaultAnalysisName="myanalysis";
const char* includePath="-I. -I$ROOTSYS/include -I$ALICE_ROOT/include";
const char* libraryDependencies=
  "libSTEERBase.so "
  "libESD.so "
  "libAOD.so "
  "libANALYSIS.so "
  "libANALYSISalice.so "
  ;

TString BuildJDLName(const char* analysisName);
AliAnalysisManager* LoadAnalysisManager(const char* filename);
TString GetIncludeHeaders(const char* filename, TString& headers, TString& libs, bool loadClass=true);
TObject* BuildCodeTree(const char* filename, TObject* pTree);
int  ProcessCodeTree(TObject* tree, TString& sources, TString& headers, TString& libs);
void ErrorConfigurationFile(const char* fileName);

void run_single_task(const char* mode,
		     const char* input,
		     const char* tasknames=NULL,
		     const char* analysisName=defaultAnalysisName,
		     const char* arguments="",
		     int nevents=-1,
		     const char* gridDataDir=NULL,
		     const char* dataPattern=NULL,
		     const char* friendDataPattern=NULL,
		     TString odir="",
		     const char* user=NULL
		     )
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // defaults
  //
  if (analysisName==defaultAnalysisName && gDirectory!=NULL) {
    // NOTE: the direct pointer comparison is on purpose
    // string comparison not necessary in this special case

    // fetch the analysis name from the setup file
    const char* confObjectName="analysis_name";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      analysisName=confObject->GetTitle();
    }
  }

  bool bRunLocal=strcmp(mode, "local")==0;
  bool bCompileOnly=strcmp(mode, "compile")==0;
  if (bCompileOnly) {
    bRunLocal=true;
    if (tasknames==NULL) {
      // short form with tasknames as the second argument
      tasknames=input;
    }
  }
  int mergeMode=0;
  if ((strcmp(mode, "merge")==0 && (mergeMode=1)>0) ||
      (strcmp(mode, "collect")==0 && (mergeMode=2)>0)) {
    mode="terminate";
    odir=input;
    input="0";
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // argument settings
  //
  bool bRunAnalysis=true;
  bool bDefaultOutput=true;
  bool bMergeOnGrid=mergeMode==2?false:true; // default true in all cases but 'collect'
  bool mcData=false;
  int nTestFiles=10;
  int nMaxInputFiles=100;
  TString strArguments(arguments);
  TString mergeDirs;
  TObjArray* tokens=strArguments.Tokenize(" ");
  if (tokens) {
    for (int iToken=0; iToken<tokens->GetEntriesFast(); iToken++) {
      TObject* token=tokens->At(iToken);
      if (!token) continue;
      TString arg=token->GetName();
      const char* key=NULL;

      if (arg.CompareTo("--help")==0 ||
          arg.CompareTo("-h")==0 ||
          arg.CompareTo("help")==0 ||
          arg.CompareTo("options")==0) {
	// printing help when called without arguments
	run_single_task();
	return;
      }
      key="--mcData";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	// switch indicates that the input data is mc data, the run numbers have
	// a different format in real data
	// NOTE: not to be confused with option 'mc' which is propagated to tasks
	// and switches processing and output modes inside tasks
	mcData=true;
	continue;
      }
      key="--merge=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	if (arg.CompareTo("local")==0) {
	  // download all files and merge locally
	  bMergeOnGrid=false;
	} else if (arg.CompareTo("collect")==0) {
	  // download the output of merging on Grid
	  // macro must have been called in mode "terminate" with option
	  // --merge=grid and the correct grid working directory
	  bMergeOnGrid=false;
	} else if (arg.CompareTo("grid")==0) {
	  // merge output on grid,  the correct grid working directory
	  // must be provided
	  bMergeOnGrid=true;
	}
	continue;
      }
      key="--nTestFiles=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	nTestFiles=arg.Atoi();
	continue;
      }
      key="--noDefaultOutput";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	bDefaultOutput=false;
	continue;
      }
      key="--stopBeforeRunning";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	bRunAnalysis=false;
	continue;
      }
      key="--maxInputFiles=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	nMaxInputFiles=arg.Atoi();
	continue;
      }
      if (!arg.BeginsWith("-") && mergeMode>0) {
	// treat as subdirectories in the remote grid dir
	mergeDirs+=" "; mergeDirs+=arg;
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
      }
    }
    delete tokens;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // load task classes and find and load all dependencies
  //
  gSystem->AddIncludePath(includePath);
  TString libraries=libraryDependencies;
  if (gDirectory!=NULL) {
    // fetch the analysis libraries from the setup file
    const char* confObjectName="analysis_libraries";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      TString analysisLibraries(confObject->GetTitle());
      TObjArray* pTokens=analysisLibraries.Tokenize(" ");
      if (pTokens) {
	for (int i=0; i<pTokens->GetEntriesFast(); i++) {
	  if (libraries.Contains(pTokens->At(i)->GetName())==0) {
	    libraries+=" ";
	    libraries+=pTokens->At(i)->GetName();
	    TString library=pTokens->At(i)->GetName();
	    if (!library.EndsWith(".so")) libraries+=".so";
	  }
	}
	delete pTokens;
      }
    }
  }
  TObjArray* pTokens=libraries.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      TString library=pTokens->At(i)->GetName();
      if (!library.EndsWith(".so")) {
	cerr << "libraries need to have ending '.so' in order to be correctly processed by alien plugin, please correct library name '" << library << "'" << endl;
      }
      if (gSystem->Load(pTokens->At(i)->GetName())==0) {
	cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }

  TString taskNames=tasknames;
  TString taskClasses="";
  TString taskSources="";
  TString taskHeaders="";
  TString addTaskMacros="";
  TString dependencyHeader;
  TString dependencySource;
  TString parPackages="";
  TString delimiter(" ");
  TStringToken taskNameTokens(taskNames, delimiter);
  TObject* pCodeTree=NULL;
  {
    while (taskNameTokens.NextToken()) {
      TString taskSource(taskNameTokens);
      TString taskHeader(taskNameTokens);
      bool bIsAddTask=false;
      if (taskSource.EndsWith(".C")) {
	// suppose that's an 'AddTask' macro
	taskHeader="";
	bIsAddTask=true;
      } else if (taskSource.EndsWith(".par")) {
	// par file
	if (gSystem->AccessPathName(taskSource)!=0) {
	  ::Error("run_single_task", Form("par file '%s' not found in current directory, you might want to set a symbolic link", taskSource.Data()));
	  return;
	}
	parPackages+=" ";
	parPackages+=taskSource;
	continue;
      } else if (taskSource.EndsWith(".h")) {
	taskSource.ReplaceAll(".h", "");
	taskClasses+=" ";
	taskClasses+=taskSource;
	taskSource+=".cxx";
      } else if (taskSource.EndsWith(".cxx")) {
	taskHeader.ReplaceAll(".cxx", "");
	taskClasses+=" ";
	taskClasses+=taskHeader;
	taskHeader+=".h";
      } else {
	taskClasses+=" ";
	taskClasses+=taskSource;
	taskSource+=".cxx";
	taskHeader+=".h";
      }
      if (gSystem->AccessPathName(taskSource)==0) {
	pCodeTree=BuildCodeTree(taskSource, pCodeTree);
	if (!bIsAddTask) {taskSources+=" "; taskSources+=taskSource;}
	else {addTaskMacros+=" "; addTaskMacros+=taskSource;}
      }
      if (gSystem->AccessPathName(taskHeader)==0) {
	pCodeTree=BuildCodeTree(taskHeader, pCodeTree);
	taskHeaders+=" "; taskHeaders+=taskHeader;
      }
    }
  }
  ProcessCodeTree(pCodeTree, dependencySource, dependencyHeader, libraries);
  if (bCompileOnly) return;

  cout << "Tasks: " << taskClasses << endl;
  cout << "Task files: " << taskSources << addTaskMacros << taskHeaders << endl;
  cout << "Dependency classes: " << dependencySource << endl;
  cout << "Dependency headers: " << dependencyHeader << endl;
  cout << "Dependency libraries: " << libraries << endl;
  cout << "Packages: " << parPackages << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // grid defaults
  //
  const char* gridConfigFile="grid-config.C";
  TString strGridConfigFile=gridConfigFile;
  if (gSystem->AccessPathName(strGridConfigFile)!=0) {
    strGridConfigFile.Prepend("/");
    strGridConfigFile.Prepend(gSystem->Getenv("HOME"));
    if (gSystem->AccessPathName(strGridConfigFile)!=0) {
      if (!bRunLocal) {
	ErrorConfigurationFile(gridConfigFile);
	return;
      }
      strGridConfigFile="";
    }
  }

  if (strGridConfigFile.IsNull()==0 && !bRunLocal) {
    cout << "loading grid configuration from file '" << strGridConfigFile << "':" << endl;
    gROOT->LoadMacro(strGridConfigFile);
    cout << " alienAPIVersion          =" << alienAPIVersion     << endl;
    cout << " alienROOTVersion         =" << alienROOTVersion    << endl;
    cout << " alienAliROOTVersion      =" << alienAliROOTVersion << endl;
    cout << " defaultGridDataDir       =" << defaultGridDataDir  << endl;
    cout << " defaultDataPattern       =" << defaultDataPattern  << endl;
    cout << " defaultFriendDataPattern =" << defaultFriendDataPattern  << endl;

    if (gridDataDir==NULL) gridDataDir=defaultGridDataDir;
    if (dataPattern==NULL) dataPattern=defaultDataPattern;
    if (friendDataPattern==NULL) friendDataPattern=defaultFriendDataPattern;
  } else if (bRunLocal) {
    if (dataPattern==NULL) {
      // thats a very crude logic, I guess it can fail in some special cases
      TString strin=input;
      if (strin.Contains("AOD"))
	dataPattern="AOD";
      else if (strin.Contains("ESD"))
	dataPattern="ESD";
    }
  }

  if(!bRunLocal) {
    // Connect to AliEn
    TGrid::Connect("alien://");
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // make the analysis manager
  //
  AliAnalysisManager *pManager=NULL;
  pManager=new AliAnalysisManager("AnalysisManager");
  if (!pManager) {
    cerr << "failed to create AnalysisManager" << endl;
    return;
  }
  AliInputEventHandler *pInputHandler = NULL;
  TString strDataPattern(dataPattern);
  if (strDataPattern.Contains("AOD")) pInputHandler=new AliAODInputHandler;
  else if (strDataPattern.Contains("ESD")) pInputHandler=new AliESDInputHandler;
  else {
    cerr << "can not determine input type from data pattern '" << dataPattern << "'" << endl;
    return;
  }
  if (!pInputHandler) {
    cerr << "failed to created input handler" << endl;
    return;
  }
  //pInputHandler->SetReadFriends(kFALSE);
  pManager->SetInputEventHandler(pInputHandler);  
  pManager->SetNSysInfo(1000);

  TString ofile=Form("%s.root", analysisName);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init for local or GRID analysis
  //
  AliAnalysisAlien *alienHandler = NULL; // for grid analysis
  TChain *chain=NULL; // for local analysis
  TString strInput=input;
  if (bRunLocal) {
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // local analysis
    //
    if (strInput.BeginsWith("alien://")) {
      // file on Grid -> connect to AliEn
      TGrid::Connect("alien://");
    }
    if(strInput.EndsWith("AliESDs.root")){
      // suppose it's a single ESD file
      chain = new TChain("esdTree"); 
      chain->Add(strInput);
    } else if(strInput.EndsWith(".txt")) {
      // Constructs chain from filenames in *.txt
      // in the form $DIR/AliESDs.root  
      gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
      // chain can contain up to 200 files, value can be modified to 
      // include a subset of what the *txt file contains
      chain = CreateESDChain(strInput.Data(),200); 

      // check if the files are on grid
      TIter next(chain->GetListOfFiles());
      TChainElement *chEl = 0;
      while(( chEl = (TChainElement*)next() )){
	TString tmp = chEl->GetTitle();	    
	if(tmp.BeginsWith("alien://")) {
	  TGrid::Connect("alien://");
	  break;
	}
      }
    } else if(strInput.EndsWith("ESD")){
      // fetch esd tree from the setup macro
      const char* esdTreeName="esdTree";
      if (gDirectory!=NULL) {
	TObject* chainObject=gDirectory->FindObject(esdTreeName);
	if (chainObject) {
	  chain=dynamic_cast<TChain*>(chainObject);
	}
      }
      if (!chain) {
	::Error("run_single_task", Form("failed to fetch esd tree object from setup; the chain with name '%s' has to be created before calling this macro", esdTreeName));
	return;
      }
    } else if(strInput.EndsWith("AliAOD.root")){
      // single local AOD file
      chain = new TChain("aodTree"); 
      chain->Add(strInput);
    } else if(strInput.EndsWith("AOD")){
      // fetch aod tree from the setup macro
      const char* aodTreeName="aodTree";
      if (gDirectory!=NULL) {
	TObject* chainObject=gDirectory->FindObject(aodTreeName);
	if (chainObject) {
	  chain=dynamic_cast<TChain*>(chainObject);
	}
      }
      if (!chain) {
	::Error("run_single_task", Form("failed to fetch aod tree object from setup; the chain with name '%s' has to be created before calling this macro", aodTreeName));
	return;
      }
    } else {
      ::Error("run_single_task", Form("invalid input"));
      return;
    }
  } else {
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // grid analysis
    //
    bool bSetRun=true;
    if (!strInput.IsDigit()) {
      // support for external macros specifying the the runs to be
      // analyzed
      // the input is expected to be an external plugin with name 'input'
      // and all run numbers being set
      TObject* pObj=gDirectory->FindObject(input);
      if (pObj) alienHandler=dynamic_cast<AliAnalysisAlien*>(pObj);
      if (!alienHandler) {
	::Error("run_single_task", Form("can not find plugin of name '%s', please setup alien handler with name and run numbers before calling this macro", input));
	return;
      }
      bSetRun=false;
    } else {
      alienHandler=new AliAnalysisAlien();
    }
    if (!alienHandler) {
      ::Error("run_single_task", Form("failed to create alien handler"));
      return;
    }

    // do not check for copying to grid (CLOSE_SE)
    alienHandler->SetCheckCopy(kFALSE);

    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    alienHandler->SetRunMode(mode);

    // number of files in test mode configurable via argument '--nTestFiles='
    if(mode=="test") alienHandler->SetNtestFiles(nTestFiles);
  
    // check the versions available on alien with the command 'packages'
    alienHandler->SetAPIVersion(alienAPIVersion);
    alienHandler->SetROOTVersion(alienROOTVersion);
    alienHandler->SetAliROOTVersion(alienAliROOTVersion);

    // using only default output
    // the alien plugin automatically recognizes all output files associated to output
    // containers, all files are treated in the standard output and added to the
    // root-archieve.root, which also seems to be needed for merging on Grid
    // see further down for using non-default output
    alienHandler->SetDefaultOutputs(bDefaultOutput);

    if (user && user[0]!=0) alienHandler->SetUser(user);

    // data alien directory
    alienHandler->SetGridDataDir(gridDataDir);
  
    // Set data search pattern
    alienHandler->SetDataPattern(dataPattern);
    alienHandler->SetFriendChainName(friendDataPattern);

    TObjArray* packageTokens=parPackages.Tokenize(" " );
    if (packageTokens) {
      for (int iPackageToken=0; iPackageToken<packageTokens->GetEntriesFast(); iPackageToken++) {
    	alienHandler->EnablePackage(packageTokens->At(iPackageToken)->GetName());
      }
      delete packageTokens;
    }

    if (bSetRun) {
      // only set if input is a run number
      if (!mcData && !strInput.BeginsWith("000"))
	alienHandler->SetRunPrefix("000");   // real data

      alienHandler->AddRunNumber(input);
    }

    if (mergeMode>0) {
      // the merge and collect modes have been added to simplify merging on grid
      // the treatment of arguments are a bit different in order to reduce list
      // of required arguments
      TString delimiter(" ");
      if (mergeDirs.IsNull()) mergeDirs="000";
      TStringToken dir(mergeDirs, delimiter);
      while (dir.NextToken()) {
	alienHandler->AddDataFile(dir->Data());
      }
      // use the specified directory names rather than a counter
      alienHandler->SetOutputToRunNo();
    }

    // define working and output directories
    TDatime dt;
    if(odir.IsNull())
      odir=(Form("gridwork/%04d-%02d-%02d_%02d-%02d", dt.GetYear(), dt.GetMonth(), dt.GetDay(), dt.GetHour(), dt.GetMinute()));
    cout << odir << endl;
    alienHandler->SetGridWorkingDir(odir); // relative to $HOME
    alienHandler->SetGridOutputDir("output");   // relative to working dir
    //alienHandler->SetOverwriteMode();                // overwrites the contents of the working and output directory

    // workaround for a Root feature: GetIncludePath() appends always
    // the current Root include path including escaped quotes. Those
    // quotes make it difficult to pass the output directly. Search for the
    // last appended include path and truncate
    TString strIncludePath(gSystem->GetIncludePath());
    Int_t pos=strIncludePath.Index(includePath);
    if (pos>=0) {
      Int_t cut=0;
      do {
	cut=pos+strlen(includePath);
      } while ((pos=strIncludePath.Index(includePath, cut))>cut);
      strIncludePath.Resize(cut);
    }
    alienHandler->AddIncludePath(strIncludePath);

    // Note: there is no extra source or header file to be transferred if 'AddTask' macros are used
    alienHandler->SetAnalysisSource(Form("%s %s %s %s", dependencySource.Data(), dependencyHeader.Data(), taskSources.Data(), taskHeaders.Data()));
    alienHandler->SetAdditionalLibs(Form("%s %s %s %s %s", libraries.Data(), dependencySource.Data(), dependencyHeader.Data(), taskSources.Data(), taskHeaders.Data()));

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    TString macroName; macroName.Form("run_%s.C",analysisName); macroName.ReplaceAll("-","_");
    alienHandler->SetAnalysisMacro(macroName);
  
    //alienHandler->SetExecutable("comparison.sh");
    alienHandler->SetExecutable(Form("run_%s.sh",analysisName));

    alienHandler->SetSplitMaxInputFileNumber(nMaxInputFiles);
  
    // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
    //    alienHandler->SetMaxInitFailed(10);
  
    // Optionally resubmit threshold.
    alienHandler->SetMasterResubmitThreshold(90); // in %

    alienHandler->SetTTL(86400);// in sec
  
    // Optionally set input format (default xml-single)
    alienHandler->SetInputFormat("xml-single");
 
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    alienHandler->SetJDLName(BuildJDLName(analysisName));
 
    // Optionally modify job price (default 1)
    alienHandler->SetPrice(1);
  
    // Optionally modify split mode (default 'se')
    alienHandler->SetSplitMode("se");
  
    // configure merging on grid,
    // argument '--merge=collect' sets 'false' for fetching the merged output
    alienHandler->SetMergeViaJDL(bMergeOnGrid); 

    alienHandler->SetOneStageMerging(kFALSE);
    alienHandler->SetMaxMergeStages(2);
  }

  // Connect plugin to the analysis manager
  if (alienHandler) {
    pManager->SetGridHandler(alienHandler);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // create task from the name, create output container, connect slots
  //
  TObjArray* taskClassTokens=taskClasses.Tokenize(" ");
  if (taskClassTokens) {
    for (int iTaskClassToken=0; iTaskClassToken<taskClassTokens->GetEntriesFast(); iTaskClassToken++) {
      AliAnalysisTaskSE *pTask=NULL;
      TString taskName=taskClassTokens->At(iTaskClassToken)->GetName();
      taskName.ReplaceAll(".cxx", "");
      TClass* pCl=TClass::GetClass(taskName);
      if (!pCl) {
	cerr << "can not load class " << taskName << endl;
	return;
      }
      void* p=pCl->New();
      if (!p) {
	cerr << "failed to instantiate class " << taskName << endl;
	return;
      }
      pTask=reinterpret_cast<AliAnalysisTaskSE*>(p);
      pManager->AddTask(pTask);
      AliAnalysisDataContainer *pContainer=pManager->CreateContainer(analysisName ,TObject::Class(), AliAnalysisManager::kOutputContainer, ofile);       
      pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
      pManager->ConnectOutput(pTask,1,pContainer);
    }
    delete taskClassTokens;
  }
  TObjArray* taskMacroTokens=addTaskMacros.Tokenize(" ");
  if (taskMacroTokens) {
    for (int iTaskMacroToken=0; iTaskMacroToken<taskMacroTokens->GetEntriesFast(); iTaskMacroToken++) {
      TString taskSource= taskMacroTokens->At(iTaskMacroToken)->GetName();

      taskSource+="+g";
      TString configuration;
      if(!strArguments.Contains("file=")) configuration+=Form(" file=%s",ofile.Data()); 
      if(!strArguments.Contains("name=")) configuration+=Form(" name=%s",analysisName); 
      configuration+=" "; configuration+=strArguments.Data();
      if (gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", configuration.Data()));
      gROOT->Macro(taskMacroTokens->At(iTaskMacroToken)->GetName());
    }
    delete taskMacroTokens;
  }

  if (!bDefaultOutput) {
    // fetch all output files from the output containers
    TString ofiles;
    TIter nextcontainer(pManager->GetContainers());
    TObject* objContainer=NULL;
    while ((objContainer=nextcontainer())!=NULL) {
      AliAnalysisDataContainer* container=dynamic_cast<AliAnalysisDataContainer*>(objContainer);
      if (!container) continue;
      ofiles+=container->GetFileName();
      ofiles+=" ";
    }

    alienHandler->SetOutputFiles(ofiles);
    // Optionally define the files to be archived.
    alienHandler->SetOutputArchive("log_archive.zip:stdout,stderr");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // run
  //

  if (!pManager->InitAnalysis()) {
    cerr << "failed to initialize analysis" << endl;
    return;
  }
  if (nevents<0) nevents=1000000000;
  pManager->PrintStatus();
  if (!bRunAnalysis) return;
  if (bRunLocal) {
    pManager->StartAnalysis("local", chain, nevents);
  } else {
    pManager->StartAnalysis("grid", nevents);
  }
}

void run_single_task()
{
  // Print help
  cout << "\nrun-single-task.C: Helper macro to run a single task either locally or on Grid"
       << "\nUsage:"
       << "\naliroot -b -q -l run-single-task.C'(\"mode\", \"input\", \"tasks\", \"name\", \"options\", events, \"path\", \"pattern\", \"friendPattern\", \"outputDir\", \"user\")' "
       << "\n arguments"
       << "\n  mode:    local, full, test, compile, merge, collect"
       << "\n  input:   In grid mode a list of run numbers. "
       << "\n           In local mode:"
       << "\n            - ESD file AliESDs.root"
       << "\n            - a list of ESD files in a text file <choose-file-name>.txt"
       << "\n            - AOD file AliAOD.root"
       << "\n            - AODs with predefined list/AODs in same folder, specifiy as \"AOD\""
       << "\n"
       << "\n  tasks:   list of class names/source or header files of tasks, or AddTask-macros"
       << "\n"
       << "\n optional arguments"
       << "\n  name:    analysis name (default 'myanalysis')"
       << "\n  options: optional arguments passed onto the task, e.g. 'mc', 'event-mixing'"
       << "\n           options of run-single-task:"
       << "\n             'mcData' -> the run numbers indicate MC Data, no '000' prepended"
       << "\n  events:  number of events to be processed (default -1 -> all)"
       << "\n  path:    data search path for grid analysis (default from configuration file)"
       << "\n  pattern: data search pattern (default from configuration file)"
       << "\n  friend pattern: friend file search pattern (default from configuration file)"
       << "\n  output dir: output directory in Grid home (default gridwork/yyyy-mm-dd_hh-mm)"
       << "\n  user:    default NULL, using user of active token"
       << "\n" << endl;
  cout << "Examples:"
       << "\naliroot -b -q -l run-single-task.C'(\"full\", \"146860\", \"AliAnalysisTaskSample\", \"myanalysis_LHC11a\")'"
       << "\n"
       << "\naliroot -b -q -l run-single-task.C'(\"local\", \"$ALICE_ROOT/test/ppbench/AliESDs.root\", \"AliAnalysisTaskSample\")'"
       << "\n"
       << "\naliroot -b -q -l run-single-task.C'(\"local\", \"AOD\", \"AddTaskSample.C\")'"
       << "\n"
       << "\naliroot -b -q -l run-single-task.C'(\"full\", \"146860\", \"AliAnalysisTaskSample\", \"correlation3p_LHC11a\", 0, -1, \"/alice/data/2011/LHC11a\", \"*/pass2_without_SDD/AOD*/*/AliAOD.root\")'"
       << "\n"
       << "\naliroot -b -q -l run-single-task.C'(\"merge\", \"gridwork/mydir\", \"AliAnalysisTaskSample\", \"myanalysis_LHC11a\")'"
       << "\n"
       << "\naliroot -b -q -l run-single-task.C'(\"collect\", \"gridwork/mydir\", \"AliAnalysisTaskSample\", \"myanalysis_LHC11a\")'"
       << "\n" << endl;

  cout << "Further options: \n" 
       << "--merge=local/grid/collect (merge option when running in mode 'terminate', simplified by runniing modes 'merge' and 'collect')\n"
       << "                           (if you want to merge files on grid, run with --merge=grid and --merge=collect to fetch files)\n"
       << "--mcData                   (needed if running on MC)\n"
       << "--nTestFiles=              (number of test files to use from grid, default=10)\n"
       << "--maxInputFiles=           (number of files in each subjob on grid, default=100)\n"
       << "--noDefaultOutput \n"
       << "--stopBeforeRunning        (Testmode for run-single-task, will stop right before starting the analysis)\n\n"
       << "To get the keywords to send to the AddTask-macros, run them individually with argument help\n\n";
}

// method for backward compatibility
void run_single_task(const char* mode,
		     const char* input,
		     const char* tasknames,
		     const char* analysisName,
		     Bool_t useMC,
		     int nevents=-1,
		     const char* gridDataDir=NULL,
		     const char* dataPattern=NULL,
		     const char* friendDataPattern=NULL,
		     TString odir="",
		     const char* user=NULL
		     )
{
  run_single_task(mode,
		  input,
		  tasknames,
		  analysisName,
		  (useMC?"mc":""),
		  nevents,
		  gridDataDir,
		  dataPattern,
		  friendDataPattern,
		  odir,
		  user
		  );
}

// method for calling with a fixed working directory, e.g. in mode terminate 
void run_single_task(const char* mode,
		     const char* input,
		     const char* tasknames,
		     const char* analysisName,
		     const char* arguments,
		     const char* workdir
		     )
{
  TString odir(workdir);
  run_single_task(mode,
		  input,
		  tasknames,
		  analysisName,
		  arguments,
		  -1,
		  NULL,
		  NULL,
		  NULL,
		  odir,
		  NULL
		  );
}

TString GetIncludeHeaders(const char* filename, TString& headers, TString& libs, bool loadClass)
{
  // scan the file and add all include headers found by path
  // to the parameter headers
  ifstream input(filename);
  if (input.bad()) {
    cerr << "failed to open file " << filename << endl;
    return headers;
  }
  TString line; 
  while (!line.ReadLine(input).eof()) {
    if (!line.Contains("#include") || !line.Contains(".h")) continue;
    line=line(0, line.Index(".h"));line+=".h";
    line.Replace(0, line.Index("#include"), "");
    line.ReplaceAll("#include", "");
    line.ReplaceAll(" ", "");
    line.ReplaceAll("\"", "");
    if (!line.BeginsWith("Ali") && !line.BeginsWith("T")) continue;
    if (gSystem->AccessPathName(line)!=0) {
      // not an include file in the current directory, check if class
      // is available or find library
      line.ReplaceAll(".h","");
      //cout << "checking class " << line << endl;
      if (TClass::GetClass(line)==NULL) {
	TString command;
	TString resfilename(gSystem->TempDirectory()); resfilename+="/findlib.txt";
	command.Form("for lib in $ALICE_ROOT/lib/*/lib*.so; do (nm $lib | grep %s | grep ' T ' | grep Class_Name > /dev/null) && echo $lib > %s; done", line.Data(), resfilename.Data());
	gSystem->Exec(command);
	ifstream resfile(resfilename.Data());
	if (resfile.good()) {
	  TString result;
	  if (!result.ReadLine(resfile).eof()) {
	    Ssiz_t haveSlash=-1;
	    while ((haveSlash=result.First('/'))>=0) result.Replace(0, haveSlash+1, "");
	    if (!libs.Contains(result)) {
	      cout << "loading dependency library '" << result << "' for class '" << line << "'" << endl;
	      gSystem->Load(result);
	      if (!libs.IsNull()) libs+=" ";
	      libs+=result;
	    }
	  }
	  command="rm "; command+=resfilename;
	  gSystem->Exec(command);
	}
      }
    } else {
      if (headers.Contains(line)) {
        if (!headers.BeginsWith(line)) {
          headers.ReplaceAll(line, "");
          if (!headers.IsNull()) headers.Insert(0, " ");
	  if (macroVerbosity>0) cout << "moving " << line << endl;
          headers.Insert(0, line);
        }
        continue;
      }
      if (!headers.IsNull()) headers.Insert(0, " ");
      if (macroVerbosity>0) cout << "inserting " << line << endl;
      headers.Insert(0, line);
      TString source=line; source.ReplaceAll(".h", ".cxx");
      if (gSystem->AccessPathName(source)==0) {
	GetIncludeHeaders(source, headers, libs);
      }
      GetIncludeHeaders(line, headers, libs);
      if (loadClass && gSystem->AccessPathName(source)==0) {
	line.ReplaceAll(".h", "");
	if (TClass::GetClass(line)==NULL) {
	  source+="+g";
	  gROOT->LoadMacro(source);
	}
      }
    }
  }

  return headers;
}

TString BuildJDLName(const char* analysisName)
{
  TString jdlname(Form("run_%s.jdl",(analysisName!=NULL?analysisName:"analysis")));
  return jdlname;
}

AliAnalysisManager* LoadAnalysisManager(const char* filename)
{
  // open file and loop through objects to find AnalysisManager
  TFile* file=TFile::Open(infilename);
  if (!file || file->IsZombie()) {
    return;
  }

  TList* keys=file->GetListOfKeys();
  if (!keys || keys->GetEntries()==0) {
    cerr << "can not find any keys in file " << infilename << endl;
    return;
  }

  TObject* pObj=NULL;
  TObject* pKey=NULL;
  TList output;
  TIter nextkey(keys);
  while (pKey=nextkey()) {
    file->GetObject(pKey->GetName(), pObj);
    if (pObj && pObj->IsA()!=AliAnalysisManager::Class())
      return dynamic_cast<AliAnalysisManager*>(pObj);
  }
}

void ErrorConfigurationFile(const char* fileName) {
  cout << endl;
  cout << "/// -------------------------------------------------------------------" << endl;
  cout << "/// Warning: can not find configuration file '" << fileName << "'" << endl;
  cout << "/// please create a configuration file in either local or HOME directory, or in" << endl;
  cout << "/// specified location. Below is an example, fill in your preferred defaults." << endl;
  cout << "/// -------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "const char* alienAPIVersion=\"V1.1x\";" << endl;
  cout << "const char* alienROOTVersion=\"v5-33-02a\";" << endl;
  cout << "const char* alienAliROOTVersion=\"v5-01-Rev-29\";" << endl;
  cout << "const char* defaultGridDataDir=\"/alice/data/2011/LHC11f\";" << endl;
  cout << "const char* defaultDataPattern=\"*ESDs.root\";" << endl;
  cout << "const char* defaultFriendDataPattern=\"\";" << endl;
  cout << "{} // note this empty body";
  cout << endl;
}

class AliCodeNode : public TNamed
{
public:
  AliCodeNode();
  AliCodeNode(const char* filename);
  ~AliCodeNode();

  enum {
    kTypeInvalid = 0,
    kTypeSource = 1,
    kTypeHeader,
    kTypeMacro,
    kNofTypes
  };
    
  const TList& GetParents() const {return fParents;}
  const TList& GetChilds() const {return fChilds;}
  int AddParent(AliCodeNode* pParent);
  int InsertChild(AliCodeNode* pChild);
  bool HasChilds() const {return (GetChilds().GetEntries()-fProcessedChilds.GetEntries())>0;}
  int MarkProcessed();
  int MarkChildProcessed(AliCodeNode* pChild);
  bool HasSourceParent();
  //int DisconnectParents();

  bool IsHeader() const {return fType==kTypeHeader;}
  bool IsSource() const {return fType==kTypeSource;}
  void HaveFile(bool haveFile) {fHaveFile=haveFile;}
  bool HaveFile() const {return fHaveFile;}
  void Print(Option_t *option="") const;

protected:
  //int DisconnectChild(const AliCodeNode& child);

private:
  TList fParents; // list of parents
  TList fChilds;  // list of childs
  TList fProcessedChilds;  // list of processed childs
  int   fType;    // source of header
  bool  fHaveFile;// file is existing in pwd

  ClassDef(AliCodeNode, 1)
};

class AliCodeTree : public TObject
{
public:
  AliCodeTree(short verbosity=0) : fNodes(), fIndentCount(0), fVerbosity(verbosity) {fNodes.SetOwner(kTRUE);}
  ~AliCodeTree() {}

  AliCodeNode* Build(const char* topfile, AliCodeNode* parent=NULL);
  AliCodeNode* FindNode(const char* name);
  int Sort();
  int LoadClasses(TString& libs);
  int LoadClasses() {TString dummy; return LoadClasses(dummy);}
  int GetHeaders(TString& headers);
  int GetSources(TString& sources);

  void Print(Option_t *option="") const;

private:
  TObjArray fNodes; // list of nodes
  short fIndentCount;
  short fVerbosity;

  ClassDef(AliCodeTree, 1)
};

ClassImp(AliCodeNode)

AliCodeNode::AliCodeNode()
 : TNamed()
 , fParents()
 , fChilds()
 , fProcessedChilds()
 , fType(AliCodeNode::kTypeInvalid)
 , fHaveFile(false)
{
}

AliCodeNode::AliCodeNode(const char* filename)
  : TNamed(filename, filename)
  , fParents()
  , fChilds()
  , fProcessedChilds()
  , fType(AliCodeNode::kTypeInvalid)
 , fHaveFile(false)
{
  TString s(filename);
  if (s.EndsWith(".cxx")) fType=kTypeSource;
  else if (s.EndsWith(".h")) fType=kTypeHeader;
  else if (s.EndsWith(".C")) fType=kTypeMacro;
}

AliCodeNode::~AliCodeNode()
{
}

int AliCodeNode::AddParent(AliCodeNode* pParent)
{
  if (!pParent) return -1;
  if (fParents.FindObject(pParent)) return 0;
  if (GetChilds().FindObject(pParent)) {
    cerr << "error: circular dependency: can not add " << pParent->GetName() << " as parent to " << this->GetName() << endl;
    return -2;
  }
  fParents.Add(pParent);
  return 0;
}

int AliCodeNode::InsertChild(AliCodeNode* pChild)
{
  if (!pChild) return -1;
  if (fChilds.FindObject(pChild)) return 0;
  if (pChild->GetChilds().FindObject(this)) {
    cerr << "error: circular dependency: can not add " << pChild->GetName() << " as child to " << this->GetName() << endl;
    return -2;
  }
  fChilds.Add(pChild);
  return 0;
}

int AliCodeNode::MarkProcessed()
{
  TIter parents(&fParents);
  TObject* obj=NULL;
  while ((obj=parents())) {
    AliCodeNode* parent=dynamic_cast<AliCodeNode*>(obj);
    parent->MarkChildProcessed(this);
  }

  return 0;
}

bool AliCodeNode::HasSourceParent()
{
  if (fType!=kTypeHeader) return false;
  TString name(GetName());
  name.ReplaceAll(".h", ".cxx");
  TIter parents(&fParents);
  TObject* obj=NULL;
  while ((obj=parents())) {
    if (name.CompareTo(obj->GetName())==0)
      return true;
  }
  return false;
}

int AliCodeNode::MarkChildProcessed(AliCodeNode* pChild)
{
  if (!pChild) return -1;
  if (fChilds.FindObject(pChild)==NULL) {
    cerr << "node " << GetName() << ": failed to find child node " << pChild->GetName() << endl;
    return -1;
  }
  if (fProcessedChilds.FindObject(pChild)!=NULL) {
    cerr << "node " << GetName() << ": child node " << pChild->GetName() << " already processed" << endl;
    return 0;
  }
  fProcessedChilds.Add(pChild);
  return 0;
}

void AliCodeNode::Print(Option_t */*option*/) const
{
  cout   << "-- " << GetName() << endl;
  TObject* obj=NULL;
  cout   << "    - parents" << endl;
  TIter parents(&fParents);
  while ((obj=parents())) {
    cout << "    |- " << obj->GetName() << endl;
  }
  cout   << "    - childs" << endl;
  TIter childs(&fChilds);
  while ((obj=childs())) {
    cout << "    |- " << obj->GetName() << endl;
  }
}

ClassImp(AliCodeTree)

AliCodeNode* AliCodeTree::Build(const char* topfile, AliCodeNode* parent)
{
  // scan the file and recursively add all include headers found by path
  int iResult=0;
  AliCodeNode* node=FindNode(topfile);
  if (!node) {
    if (fVerbosity>0) cout << setw(2*fIndentCount) << " " << "new node " << topfile << endl;
    fIndentCount++;
    node=new AliCodeNode(topfile);
    fNodes.Add(node);
    ifstream input(topfile);
    if (input.good()) {
      node->HaveFile(true);
      TString line; 
      while (!line.ReadLine(input).eof()) {
	if (!line.Contains("#include") || !line.Contains(".h")) continue;
	line=line(0, line.Index(".h"));line+=".h";
	line.Replace(0, line.Index("#include"), "");
	line.ReplaceAll("#include", "");
	line.ReplaceAll(" ", "");
	line.ReplaceAll("\"", "");
	if (!line.BeginsWith("Ali") && !line.BeginsWith("T")) continue;
	AliCodeNode* child=NULL;
	TString source=line; source.ReplaceAll(".h", ".cxx");
	if (source.CompareTo(topfile)!=0 && gSystem->AccessPathName(source)==0) {
	  child=Build(source, node);
	  node->InsertChild(child);
	}
	child=Build(line, node);
	node->InsertChild(child);
      }
    }
    fIndentCount--;
  }
  if (parent) {
    if ((iResult=node->AddParent(parent))<0)
      return NULL;
  }
  if (fVerbosity>0) cout << setw(2*fIndentCount) << " " << "finished " << topfile << endl;
  return node;
}

int AliCodeTree::Sort()
{
  TObjArray sortedNodes;
  int nNodes=fNodes.GetEntriesFast();
  int nCount=0;
  while (sortedNodes.GetEntriesFast()<nNodes && nCount<nNodes) {
    for (int i=0; i<nNodes; i++) {
      if (fNodes[i]==NULL) continue;
      AliCodeNode* node=dynamic_cast<AliCodeNode*>(fNodes[i]);
      if (node->HasChilds()) {
	continue;
      }
      fNodes[i]=NULL;
      sortedNodes.Add(node);
      node->MarkProcessed();
    }
    nCount++;
  }

  for (int i=0; i<nNodes; i++) {
    fNodes[i]=sortedNodes[i];
  }

  return 0;
}

int AliCodeTree::LoadClasses(TString& libs)
{
  TIter next(&fNodes);
  TObject* obj=NULL;
  while ((obj=next())!=NULL) {
    AliCodeNode* node=dynamic_cast<AliCodeNode*>(obj);
    TString name(node->GetName());
    if (node->IsHeader()) {
      if (node->HasSourceParent()) {
	// nothing to do, class going to be compiled from source
	continue;
      }
      name.ReplaceAll(".h", "");
      if (TClass::GetClass(name)!=NULL) {
	// class available in the system
	continue;
      }

      TString command;
      TString resfilename(gSystem->TempDirectory()); resfilename+="/findlib.txt";
      command.Form("for lib in $ALICE_ROOT/lib/*/lib*.so; do (nm $lib | grep %s | grep ' T ' | grep Class_Name > /dev/null) && echo $lib > %s; done", name.Data(), resfilename.Data());
      gSystem->Exec(command);
      ifstream resfile(resfilename.Data());
      if (resfile.good()) {
	TString result;
	if (!result.ReadLine(resfile).eof()) {
	  Ssiz_t haveSlash=-1;
	  while ((haveSlash=result.First('/'))>=0) result.Replace(0, haveSlash+1, "");
	  if (!libs.Contains(result)) {
	    cout << "loading dependency library '" << result << "' for class '" << name << "'" << endl;
	    gSystem->Load(result);
	    if (!libs.IsNull()) libs+=" ";
	    libs+=result;
	  }
	}
	command="rm "; command+=resfilename;
	gSystem->Exec(command);
      }
      continue;
    }
    if (node->IsSource()) {
      TString classname(name);
      classname.ReplaceAll(".cxx", "");
      if (TClass::GetClass(classname)!=NULL) {
	// class available in the system
	continue;
      }
      name+="+g";
      gROOT->LoadMacro(name);
    }
  }
  return 0;
}

int AliCodeTree::GetHeaders(TString& headers)
{
  TIter next(&fNodes);
  TObject* obj=NULL;
  while ((obj=next())!=NULL) {
    AliCodeNode* node=dynamic_cast<AliCodeNode*>(obj);
    if (!node->IsHeader() || !node->HaveFile()) continue;
    if (!headers.IsNull()) headers+=" ";
    headers+=node->GetName();
  }
  return 0;
}

int AliCodeTree::GetSources(TString& sources)
{
  TIter next(&fNodes);
  TObject* obj=NULL;
  while ((obj=next())!=NULL) {
    AliCodeNode* node=dynamic_cast<AliCodeNode*>(obj);
    if (!node->IsSource() || !node->HaveFile()) continue;
    if (!sources.IsNull()) sources+=" ";
    sources+=node->GetName();
  }
  return 0;
}

AliCodeNode* AliCodeTree::FindNode(const char* name)
{
  TObject* node=fNodes.FindObject(name);
  if (!node) return NULL;
  return dynamic_cast<AliCodeNode*>(node);
}

void AliCodeTree::Print(Option_t *option) const
{
  const char* key=NULL;
  int indent=0;
  TString childOptions;
  const TString delimiter(" ");
  TStringToken options(option, delimiter);
  while (options.NextToken()) {
    key="indent=";
    if (options.BeginsWith(key)) {
      TString arg(options);
      arg.ReplaceAll(key, "");
      indent=arg.Atoi();
      continue;
    }
    childOptions+=" ";
    childOptions+=options;
  }
  childOptions+=Form("indent=%d", indent+1);
  TIter next(&fNodes);
  TObject* obj=NULL;
  while ((obj=next())!=NULL) {
    obj->Print(childOptions);
  }
}

TObject* BuildCodeTree(const char* filename, TObject* useObject)
{
  AliCodeTree* pTree=NULL;
  if (useObject) pTree=dynamic_cast<AliCodeTree*>(useObject);
  if (!pTree) pTree=new AliCodeTree;
  if (!pTree) return NULL;
  
  pTree->Build(filename);
  return pTree;
}

int ProcessCodeTree(TObject* tree, TString& sources, TString& headers, TString& libs)
{
  if (!tree) return -1;
  AliCodeTree* pTree=dynamic_cast<AliCodeTree*>(tree);
  pTree->Sort();
  pTree->LoadClasses(libs);
  pTree->GetHeaders(headers);
  pTree->GetSources(sources);
  return 0;
}

#else
#include "TObject.h"
#include "TNamed.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TPRegexp.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TChain.h"
#include "TChainElement.h"
// #include "AliAnalysisManager.h"
// #include "AliAnalysisAlien.h"
// #include "AliAnalysisTaskSE.h"
// #include "AliInputEventHandler.h"
// #include "AliAODInputHandler.h"
// #include "AliESDInputHandler.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#endif
