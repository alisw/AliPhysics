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
///   mode:    local, full, test
///   input:   In grid mode a list of run numbers. 
///            In local mode:
///             - ESD file AliESDs.root
///             - a list of ESD files in a text file <choose-file-name>.txt
///             - AOD file AliAOD.root
///             - AODs with predefined list/AODs in same folder, specifiy as "AOD"
///
///   tasks:   list of class names/source or header files of task, or AddTask-Macro
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
/// 
/// Suggestions:
/// Feedback appreciated: Matthias.Richter@cern.ch
/// If you find this macro useful but not applicable without patching it, let me know
/// your changes and use cases.


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

TString GetIncludeHeaders(const char* filename, TString& headers, TString& libs, bool loadClass=true);
void ErrorConfigurationFile(const char* fileName);

void run_single_task(const char* mode,
		     const char* input,
		     const char* tasknames,
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
  // argument settings
  //
  bool mcData=false;
  TString strArguments(arguments);
  TObjArray* tokens=strArguments.Tokenize(" ");
  if (tokens) {
    for (int iToken=0; iToken<tokens->GetEntriesFast(); iToken++) {
      TObject* token=tokens->At(iToken);
      if (!token) continue;
      TString arg=token->GetName();
      const char* key=NULL;

      key="mcData";
      if (arg.CompareTo(key)==0) {
	strArguments.ReplaceAll(key, "");
	mcData=true;
      }
    }
    delete tokens;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // make the analysis manager
  //
  AliAnalysisManager *pManager  = new AliAnalysisManager("AnalysisManager");
  if (!pManager) {
    cerr << "failed to created AnalysisManager" << endl;
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
  TString parPackages="";
  TObjArray* pTaskNames=taskNames.Tokenize(" ");
  if (pTaskNames) {
    for (int iTaskName=0; iTaskName<pTaskNames->GetEntriesFast(); iTaskName++) {
      TString taskSource=pTaskNames->At(iTaskName)->GetName();
      TString taskHeader=pTaskNames->At(iTaskName)->GetName();
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
      TString dependencyHeader;
      TString dependencySource;
      if (gSystem->AccessPathName(taskHeader)==0) {
	GetIncludeHeaders(taskHeader, dependencyHeader, libraries);
	taskHeaders+=" "; taskHeaders+=taskHeader;
      }
      if (gSystem->AccessPathName(taskSource)==0) {
	GetIncludeHeaders(taskSource, dependencyHeader, libraries);
	if (!bIsAddTask) {taskSources+=" "; taskSources+=taskSource;}
	else {addTaskMacros+=" "; addTaskMacros+=taskSource;}
      }
      TObjArray* pTokens=dependencyHeader.Tokenize(" ");
      if (pTokens) {
	for (int i=0; i<pTokens->GetEntriesFast(); i++) {
	  TString sourceFile=pTokens->At(i)->GetName();
	  sourceFile.ReplaceAll(".h", ".cxx");
	  if (gSystem->AccessPathName(sourceFile)!=0) continue;
	  if (!dependencySource.IsNull()) dependencySource+=" ";
	  dependencySource+=sourceFile;
	  if (!libraries.IsNull()) libraries+=" ";
	  libraries+=sourceFile;
	}
	delete pTokens;
      }
      dependencySource.ReplaceAll(taskSource, "");
      dependencyHeader.ReplaceAll(taskHeader, "");
    }
    delete pTaskNames;
  }
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
      if (gDirectory!=NULL) {
	const char* esdTreeName="esdTree";
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
      if (gDirectory!=NULL) {
	const char* aodTreeName="aodTree";
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
    TString strInput(input);
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
  
    // check the versions available on alien with the command 'packages'
    alienHandler->SetAPIVersion(alienAPIVersion);
    alienHandler->SetROOTVersion(alienROOTVersion);
    alienHandler->SetAliROOTVersion(alienAliROOTVersion);

    //Allow non-default outputs
    //This is required to set non-default output files with the SetOutputFiles
    //function below
    alienHandler->SetDefaultOutputs(kFALSE);
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
    alienHandler->SetAdditionalLibs(Form("%s %s %s", libraries.Data(), taskHeaders.Data(), dependencyHeader.Data()));

    alienHandler->SetOutputFiles(ofile);

    // Optionally define the files to be archived.
    alienHandler->SetOutputArchive("log_archive.zip:stdout,stderr");
  
    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    TString macroName; macroName.Form("run_%s.C",analysisName); macroName.ReplaceAll("-","_");
    alienHandler->SetAnalysisMacro(macroName);
  
    //alienHandler->SetExecutable("comparison.sh");
    alienHandler->SetExecutable(Form("run_%s.sh",analysisName));

    alienHandler->SetSplitMaxInputFileNumber(100);
  
    // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
    alienHandler->SetMaxInitFailed(10);
  
    // Optionally resubmit threshold.
    alienHandler->SetMasterResubmitThreshold(90); // in %

    alienHandler->SetTTL(86400);// in sec
  
    // Optionally set input format (default xml-single)
    alienHandler->SetInputFormat("xml-single");
 
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    alienHandler->SetJDLName(Form("run_%s.jdl",analysisName));
 
    // Optionally modify job price (default 1)
    alienHandler->SetPrice(1);
  
    // Optionally modify split mode (default 'se')
    alienHandler->SetSplitMode("se");
  
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    if (strcmp(mode, "terminate")==0) {
      alienHandler->SetMergeViaJDL(kFALSE);
    }

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
      TObject* p=pCl->New();
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
      taskSource+="+g";
      TString configuration;
      configuration.Form("name=%s file=%s %s", analysisName, ofile.Data(), strArguments.Data());
      if (gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", configuration.Data()));
      gROOT->Macro(taskMacroTokens->At(iTaskMacroToken)->GetName());
    }
    delete taskMacroTokens;
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
  if (bRunLocal) {
    pManager->StartAnalysis("local", chain, nevents);
  } else {
    pManager->StartAnalysis("grid", nevents);
  }
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
		  taskname,
		  analysisName,
		  useMC?"mc":"",
		  nevents,
		  gridDataDir,
		  dataPattern,
		  friendDataPattern,
		  odir,
		  user
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
