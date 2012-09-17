//-*- Mode: C++ -*-
// $Id$

/// @file   run-single-task.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-04-06
/// @brief  Run a single task
///
/// Helper macro to run a single task either locally or on Grid
/// Usage:
/// aliroot -b -q -l run-single-task.C'("mode", "run", "task", "name", useMC, events, "path", "pattern", "friendPattern", "outputDir", "user")'
///  arguments
///   mode:    local, full, test
///   run:     list of run numbers. Or if using AODs with predefined list/AODs in same folder, specifiy as "AOD"
///   task:    class name of task
///
///  optional arguments
///   name:    analysis name (default 'myanalysis')
///   events:  number of events to be processed (default -1 -> all)
///   path:    data search path for grid analysis (default from configuration file)
///   pattern: data search pattern (default from configuration file)
///
/// Examples:
/// aliroot -b -q -l run-single-task.C'("full", "146860", "AliAnalysisTaskSample", "myanalysis_LHC11a")'
///
/// aliroot -b -q -l run-single-task.C'("local", "$ALICE_ROOT/test/ppbench/AliESDs.root", "AliAnalysisTaskSample")'
///
/// aliroot -b -q -l run-single-task.C'("local", "AOD", "AddTaskSample.C"")'
///
/// aliroot -b -q -l run-single-task.C'("full", "146860", "AliAnalysisTaskSample", "correlation3p_LHC11a", -1, "/alice/data/2011/LHC11a", "*/pass2_without_SDD/AOD*/*/AliAOD.root")'
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
///
/// Local analysis:
/// requires only the path to the input file and the task class name. If the specified file is
/// a text file (.txt) each line can contain an input file path, all files are chained.
/// Note: AOD mode needs to be implemented
///
/// Grid analysis:
/// All modes provided by the AliAnalysisAlien plugin can be used, e.g. full, test, offline
/// A couple of settings need to be defined in a configuration file 'grid-config.C' which can be
/// either in the local directory or home directory.
/// <pre>
///   const char* alienAPIVersion="V1.1x";
///   const char* alienROOTVersion="v5-33-02a";
///   const char* alienAliROOTVersion="v5-01-Rev-29";
///   const char* defaultGridDataDir="/alice/data/2011/LHC11a";
///   const char* defaultDataPattern="*/pass2_without_SDD/*/AliESDs.root";
///   {} // note this empty body
/// </pre>
/// Data path and pattern can also be specified as command line arguments.
/// The working directory in the grid home directory of the user is set to
/// gridwork/<date>_<time>.
///
/// 

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
		     const char* taskname,
		     const char* analysisName=defaultAnalysisName,
		     Bool_t useMC=kFALSE,
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
      TString strin=input;
      if (strin.EndsWith("AOD"))
	dataPattern="AOD";
      else if (strin.EndsWith("ESD"))
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
	  }
	}
	delete pTokens;
      }
    }
  }
  TObjArray* pTokens=libraries.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      if (gSystem->Load(pTokens->At(i)->GetName())==0) {
	cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }

  bool bCreateAndAddTask=true;
  TString taskSource=taskname;
  TString taskHeader=taskname;
  if (taskSource.EndsWith(".C")) {
    // suppose that's an 'AddTask' macro
    taskHeader="";
    bCreateAndAddTask=false;
  } else if (taskSource.EndsWith(".h")) {
    taskSource.ReplaceAll(".h", ".cxx");
  } else if (taskSource.EndsWith(".cxx")) {
    taskHeader.ReplaceAll(".cxx", ".h");
  } else {
    taskSource+=".cxx";
    taskHeader+=".h";
  }
  TString dependencyHeader;
  TString dependencySource;
  GetIncludeHeaders(taskHeader, dependencyHeader, libraries);
  GetIncludeHeaders(taskSource, dependencyHeader, libraries);
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
  cout << "Task files: " << taskSource << " " << taskHeader << endl;
  cout << "Dependency classes: " << dependencySource << endl;
  cout << "Dependency headers: " << dependencyHeader << endl;
  cout << "Dependency libraries: " << libraries << endl;

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
      chain = CreateESDChainf(strInput.Data(),200); 

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
    } else if(strInput.EndsWith("AOD")){
      // fetch aod tree from the setup macro
      if (gDirectory!=NULL) {
	TObject* chainObject=gDirectory->FindObject("aodTree");
	if (chainObject) {
	  chain=dynamic_cast<TChain*>(chainObject);
	}
      }
      if (!chain) {
	cout << "failed to fetch aod tree object from setup" << endl;
	return -1;
      }
    } else {
      cerr << "invalid input" << endl;
      return -1;
    }
  } else {
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // grid analysis
    //
    alienHandler=new AliAnalysisAlien();
    if (!alienHandler) {
      cerr << "failed to create alien handler" << endl;
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

    if (!useMC)
      alienHandler->SetRunPrefix("000");   // real data

    alienHandler->AddRunNumber(input);

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
    alienHandler->SetAnalysisSource(Form("%s %s %s %s", dependencySource.Data(), dependencyHeader.Data(), (bCreateAndAddTask)?taskSource.Data():"", taskHeader.Data()));
    alienHandler->SetAdditionalLibs(Form("%s %s %s", libraries.Data(), taskHeader.Data(), dependencyHeader.Data()));

    alienHandler->SetOutputFiles(ofile);

    // Optionally define the files to be archived.
    //-
    //alienHandler->SetOutputArchive("log_archive.zip:stdout,stderr");
  
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

    alienHandler->SetTTL(30000);// in sec
  
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
  AliAnalysisTaskSE *pTask=NULL;
  if (bCreateAndAddTask) {
    TClass* pCl=TClass::GetClass(taskname);
    if (!pCl) {
      cerr << "can not load class " << taskname << endl;
      return -1;
    }
    TObject* p=pCl->New();
    if (!p) {
      cerr << "failed to instantiate class " << taskname << endl;
      return -1;
    }
    pTask=reinterpret_cast<AliAnalysisTaskSE*>(p);
    pManager->AddTask(pTask);
    AliAnalysisDataContainer *pContainer=pManager->CreateContainer(analysisName ,TObject::Class(), AliAnalysisManager::kOutputContainer, ofile);       
    pManager->ConnectInput(pTask,0,pManager->GetCommonInputContainer());
    pManager->ConnectOutput(pTask,1,pContainer);
  } else {
    taskSource+="+g";
    TString configuration;
    configuration.Form("name=%s file=%s %s", analysisName, ofile.Data(), useMC?"mc":"");
    if (gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", configuration.Data()));
    gROOT->Macro(taskSource);
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
