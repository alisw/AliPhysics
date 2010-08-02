// ======================================
// ===== ALIEN PLUGIN CONFIGURATION =====
// ======================================
//
// This example macro configures an AliAnalysisAlien object
// to be used to run analysis on AliEn without having to put
// all the required code and scripts in the alien shell.
//
// All the possible configuration parameters are arguments
// of the macro function, even if most of them have default
// values which the user will rarely change, and some args
// do exclude some other ones, depending of user choice.
//
// The macro tries to synchronize some output names, using 
// a unique name ('analysisName') to define all files that
// describe the output, the analysis macros/executables/JDL.
//
// Since the run mode can be more variable than the config
// it is not set here, but it is required in the run macro
// which uses the plugin.
//
// Considered that the arguments are many, they are explained
// inside the list of arguments in the macro definition.
// In ALL cases where a list of strings must be provided, its
// elements must be separated by a blank character.
//
AliAnalysisAlien* PluginDataByRun
(
  // common base name for defining those of all created files (run numbers will be appended)
  const char *analysisName = "AnalysisPhi7TeV",

  // input files to be processed:
  // -- runList   --> list of run numbers, separated by a minus (e.g.: "23211-22231-2232") --> will name the output directory
  // -- runPath   --> AliEn path where runs will be searched
  // -- runPrefix --> string prefix prepended to run numbers
  // -- pattern   --> name of the file to be used (ESD/AOD)
  const char *runList     = "117112",
  const char *runPath     = "/alice/data/2010/LHC10b",
  const char *runPrefix   = "000",
  const char *pattern     = "pass2/*/AliESDs.root",

  // working parameters:
  // -- workDir        --> AliEn working directory (w.r. to AliEn $HOME)
  // -- outDir         --> AliEn output directory (w.t. to workDir)
  // -- analysisName   --> common name used for analysis related files (.C, .JDL, exe, log)
  // -- taskSource     --> name of a task compiled on-fly (without any extension)
  // -- exePre         --> executable command part before task macro name
  // -- exePost        --> executable command part after task macro name
  // -- exeArgs        --> executable arguments to be specified in the JDL
  // -- outList        --> lit of expected output files (separated by commas)
  const char *workDir        = "rsnNewTest_v1",
  const char *outDir         = "run",
  const char *taskSource     = "",
  const char *exePre         = "aliroot -q -b",
  const char *exePost        = "",
  const char *exeArgs        = "",
  
  // additional stuff (libs, includes, code):
  // -- addInclude  --> additional include paths
  // -- addLibs     --> additional libraries compiled in ROOT or ALIROOT
  // -- addPars     --> additional PAR libraries from ALIROOT
  // -- addExternal --> external tasks to be compiled on-fly
  const char *addInclude  = "TOF",
  const char *addLibs     = "",
  const char *addPars     = "PWG2resonances.par",
  const char *addExternal = "",
  
  // job parameters:
  // these are all the job parameters which can be inserted in the JDL
  // their names are the same as the corresponding JDL keywords
  Int_t       split          = 50,
  Int_t       maxInitFailed  = 0,
  Int_t       resubmitThr    = 0,
  Int_t       TTL            = 22500,
  const char *inputFormat    = "xml-single",
  const char *splitMode      = "se",
  Int_t       price          = 1,
  const char *jobTag         = "",
  Int_t       maxMergeFiles  = 50,
  Int_t       nTestFiles     = 6,
  
  // standard package versions
  const char *apiVersion     = "V1.1x",
  const char *rootVersion    = "v5-26-00b-6",
  const char *aliVersion     = "v4-19-22-AN"
)
{
  // append the run number to the reference names used 
  // for the executeble and for naming the used files
  Char_t exeName[255];
  sprintf(exeName, "%s_exe", analysisName);
  
  // create plugin object
  AliAnalysisAlien *plugin = new AliAnalysisAlien;
  plugin->SetOutputToRunNo(kTRUE);
  plugin->SetMergeViaJDL();
  plugin->SetNtestFiles(nTestFiles);
  plugin->SetNrunsPerMaster(0);
  
  // define names of analysis files after 'analysisName'
  const char *analysisMacro = Form("%s.C"   , analysisName);
  const char *logFile       = Form("%s.log" , analysisName);
  const char *jdlFile       = Form("%s.jdl" , analysisName);
  const char *outFile       = Form("%s.root", analysisName);
  const char *exeFile       = Form("%s.sh"  , exeName);
  const char *archiveFile   = Form("%s.zip:stdout,stderr,%s,%s@disk=2", analysisName, logFile, outFile);
  
  // package versions
  plugin->SetAPIVersion(apiVersion);
  plugin->SetROOTVersion(rootVersion);
  plugin->SetAliROOTVersion(aliVersion);
  
  // work and output directory in AliEn shell
  plugin->SetGridWorkingDir(workDir);
  plugin->SetGridOutputDir(outDir);
  
  // executable command and arguments
  // in 'post' args it is added the request to log all output
  // into a file which is named after the common analysis name
  plugin->SetExecutableCommand(exePre);
  plugin->SetExecutableArgs(Form("%s >& %s", exePost, logFile));
  plugin->SetExecutable(exeFile);
  if (strlen(exeArgs) > 0) plugin->SetArguments(exeArgs);
  
  // declare (if any) the analysis source file
  // to be compiled runtime separated by blanks
  // and, in this case, add its header and implementation to additional libs
  if (strlen(taskSource) > 0)
  {
    plugin->SetAnalysisSource(Form("%s.h", taskSource));
    plugin->SetAdditionalLibs(Form("%s %s.h %s.cxx", addLibs, taskSource, taskSource));
  }
  else
    plugin->SetAdditionalLibs(addLibs);
  
  // additional PARs
  if (strlen(addPars) > 0) plugin->EnablePackage(addPars);
  
  // external packages
  if (strlen(addExternal) > 0) plugin->AddExternalPackage(addExternal);
    
  // set names of analysis macro, JDL script, output file and archive
  // from the commonly used 'outName' argument
  plugin->SetAnalysisMacro(analysisMacro);
  plugin->SetJDLName(jdlFile);
  plugin->SetDefaultOutputs(kFALSE);
  plugin->SetOutputFiles(outFile);
  plugin->SetOutputArchive(archiveFile);
  
  // add all runs
  TString sList = runList;
  list  = sList.Tokenize("-");
  Int_t i, n = list->GetEntries();
  plugin->SetRunPrefix(runPrefix);
  plugin->SetGridDataDir(runPath);
  plugin->SetDataPattern(pattern);
  for (i = 0; i < n; i++)
  {
    TObjString *os = (TObjString*)list->At(i);
    plugin->AddRunNumber(os->GetString().Atoi());
    cout << "Added run number " << os->GetString().Data() << endl;
  }
  
  // JDL parameters
  plugin->SetSplitMaxInputFileNumber(split);
  plugin->SetMaxInitFailed(maxInitFailed);
  plugin->SetMasterResubmitThreshold(resubmitThr);
  plugin->SetTTL(TTL);
  plugin->SetPrice(price);
  plugin->SetInputFormat(inputFormat);
  plugin->SetSplitMode(splitMode);
  if (strlen(jobTag) > 0) plugin->SetJobTag(jobTag);
  plugin->SetMaxMergeFiles(maxMergeFiles);
    
  // the end!
  return plugin;
}
