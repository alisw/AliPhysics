AliAnalysisGrid* CreateAlienHandler()
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetUser("jonderwa");
  plugin->SetOverwriteMode();
  // --- Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode("terminate"); // use to merge output when jobs are done
  plugin->SetRunMode("test"); // use to test code
  //plugin->SetRunMode("full"); // use to submit code
  plugin->SetNtestFiles(1);

  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");


  plugin->SetROOTVersion("v5-34-30-alice-22");
  plugin->SetAliROOTVersion("v5-08-03a-6");
  plugin->SetAliPhysicsVersion("vAN-20160414-1");
  //plugin->SetCheckCopy(kFALSE);
  //plugin->EnablePackage("PWGPPevcharQn.par");
  //plugin->EnablePackage("PWGPPevcharQnInterface.par");

  //plugin->SetInputFile("/hera/alice/jonderw/gridanalysis/CalibrationFiles/000137161/CalibrationHistograms.root");
  //plugin->SetMergeDirName("/hera/alice/jonderw/gridanalysis/CalibrationFiles/000137162/");

  plugin->SetGridDataDir("/alice/data/2010/LHC10h/");
  //plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
  //plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/");
  //plugin->SetDataPattern("ESDs/pass2/AOD086/*/AliAOD.root"); // all segments for running on Grid
  //plugin->SetDataPattern("ESDs/pass2/AOD086/0001/AliAOD.root"); // single segment for testing
  //plugin->SetDataPattern("ESDs/pass2/AliESDs.root"); // single segment for testing
  plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // 2010|2011
  //plugin->SetDataPattern("/pass1/*/AliESDs.root"); // 2015
  plugin->SetRunPrefix("000");
  //plugin->AddRunNumber(245145); // LHC15o
  //plugin->AddRunNumber(138534);
  plugin->AddRunNumber(137162);
  //plugin->AddRunNumber(138534);
  //gROOT->ProcessLine(".L AddRunNumbers.C");
  //int added = AddRunNumbers(plugin,0,1,"lhc10h"); // adds one run number for testing
  //int added = AddRunNumbers(plugin,0,5,"lhc10h"); // adds multiple run numbers for running on grid
  //int added = AddRunNumbers(plugin,0,90,"lhc10h"); // adds multiple run numbers for running on grid
  //if(added<0) return NULL;

  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("AnalysisTesting");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  plugin->SetOutputToRunNo(); // Write output to runnumber directory
  

  //plugin->SetAdditionalLibs("libSTEERBase libESD libAOD libANALYSIS libANALYSISalice libANALYSISaliceBase libCORRFW libOADB");

  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAnalysisSource("AliAnalyisTaskEventPlaneCalibration.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  //plugin->SetAdditionalLibs("AliAnalyisTaskEventPlaneCalibration.h AliAnalyisTaskEventPlaneCalibration.cxx");
  // Declare the output file names separated by blancs.
  //plugin->SetOutputFiles("dstTree.root");
  //plugin->SetDefaultOutputs();
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MyTask.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable("MyTask.sh");
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("MyTask.jdl");
  plugin->SetMergeViaJDL(kTRUE);
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  return plugin;
}
