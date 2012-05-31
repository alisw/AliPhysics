AliAnalysisGrid* CreateAnalysisPlugin(TString analysisMode="full") 
{ 
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output results from a previous session
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate") 

  plugin->SetRunMode(TString (analysisMode.Data()));

  // Set versions of used packages 
  plugin->SetAPIVersion("V1.1x"); 
  plugin->SetROOTVersion("v5-28-00f");
  plugin->SetAliROOTVersion("v4-21-33-AN"); 
  // Declare input data to be processed. 

  // Method 1: Create automatically XML collections using alien 'find' command. 
  // Define production directory LFN 
  //    plugin->SetGridDataDir("/alice/sim/LHC10a6"); 
  //*   plugin->SetGridDataDir("/alice/data/2010/LHC10b"); 
  // Set data search pattern 
  //   plugin->SetDataPattern("*ESDs.root");  // simulated, tags not used
  //*  plugin->SetDataPattern("*ESDs/pass2/*ESDs.root"); // real data check reco pass and data base directory
  //*  plugin->SetRunPrefix("000");   // real data
  //    plugin->SetDataPattern("*tag.root");  // Use ESD tags (same applies for AOD's)
  // ...then add run numbers to be considered 
  //    plugin->AddRunNumber(125020);    // simulated
							
  //for ESDs
  plugin->SetGridDataDir("/alice/data/2011/LHC11a"); 
  plugin->SetDataPattern("*ESDs/pass1/*ESDs.root");
  plugin->SetRunPrefix("000");

  //ESDs sim
  // plugin->SetGridDataDir("/alice/sim/LHC10f9b");
  //plugin->SetDataPattern("*ESDs.root"); 

  //for AODs
  // plugin->SetGridDataDir("/alice/data/2010/LHC10c"); 
  // plugin->SetRunPrefix("000");   // real data
  //plugin->SetDataPattern("*ESDs/pass2_recovery_900GeV/AOD017/*AOD.root"); 
  
  //sim AODs
  // plugin->SetGridDataDir("/alice/sim/LHC10d4a"); 
  //plugin->SetDataPattern("*AOD012/*AOD.root");
  
  //  TString runs ="120824:120823:120822:120821:120820:120758:120750:120741:120671:120617:120616:120505:120504:120503:120244:120079:120076:120073:120072:120069:120067:119862:119859:119856:119853:119849:119846:119845:119844:119842:119841:119163:119161:119159:119086:119085:119084:119079:119077:119067:119061:119047:119041:119037"; // dont forget last two runs
  
  //TString runs ="120829:120825";
  // TString runs="118506:118507:118512:118518:118556:118558:118560:118561";
  TString runs ="146801";

  TObjArray* array = runs.Tokenize ( ":" );
  TObjString *str;
  TString strr,strr2_1,strr2_2;
  for ( Int_t i = 0;i < array->GetEntriesFast();i++ ) {
    str = ( TObjString * ) array->At ( i );
    strr = str->GetString();
    if ( !strr.IsNull() ) {
      plugin->AddRunNumber(strr.Atoi());
    }
  }  
	 
  // Method 2: Declare existing data files (raw collections, xml collections, root file) 
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir()) 
  // XML collections added via this method can be combined with the first method if 
  // the content is compatible (using or not tags) 
  //   plugin->AddDataFile("tag.xml"); 
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root"); 

  // Define alien work directory where all files will be copied. Relative to alien $HOME. 
  plugin->SetGridWorkingDir("146801"); 
  // Declare alien output directory. Relative to working directory. 
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output 
  // Declare the analysis source files names separated by blancs. To be compiled runtime 
  // using ACLiC on the worker nodes. 
  plugin->SetAnalysisSource("AliAnalysisTaskEfficiency.cxx"); 
  //    plugin->SetAdditionalRootLibs("CORRFW PWG2resonances");
  //    plugin->SetAdditionalRootLibs("PWG2resonances");
  //    plugin->SetAdditionalRootLibs("PWG2resonances");
  // 
  plugin->SetAdditionalLibs("AliAnalysisTaskEfficiency.h AliAnalysisTaskEfficiency.cxx");
  //    plugin->EnablePackage("PWG2resonances");
  //    plugin->EnablePackage("");
  //    plugin->EnablePackage("");
  // Declare all libraries (other than the default ones for the framework. These will be 
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here. 

  // No need for output file names. Procedure is automatic. 
  //   plugin->SetOutputFiles("Pt.ESD.1.root"); 
  //   plugin->SetDefaultOutputs(); 
  // No need define the files to be archived. Note that this is handled automatically by the plugin.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr"); 
  // Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
  plugin->SetAnalysisMacro("AnalysisTest.C"); 
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore). The optimum for an analysis
  // is correlated with the run time - count few hours TTL per job, not minutes !
  plugin->SetSplitMaxInputFileNumber(100); 
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs. 
  plugin->SetMaxInitFailed(5); 
  // Optionally resubmit threshold. 
  plugin->SetMasterResubmitThreshold(90); 
  // Optionally set time to live (default 30000 sec) 
  plugin->SetTTL(20000); 
  // Optionally set input format (default xml-single) 
  plugin->SetInputFormat("xml-single"); 
  // Optionally modify the name of the generated JDL (default analysis.jdl) 
  plugin->SetJDLName("TaskRsn.jdl"); 
  // Optionally modify job price (default 1) 
  plugin->SetPrice(1);  
  // Optionally modify split mode (default 'se')     
  plugin->SetSplitMode("se"); 
   
   
  //++++++++++++++++ PROOF ++++++++++++++++
  //    Proof cluster
   
  plugin->SetProofCluster("alice-caf");
  //plugin->SetProofCluster("skaf.saske.sk");
  //   plugin->SetProofCluster("skaf-test.saske.sk");
  // Dataset to be used
  //     plugin->SetProofDataSet("/alice/sim/LHC10a12_104316#esdTree");
  //       plugin->SetProofDataSet("/alice/sim/LHC10a12_104157#esdTree");
  //     plugin->SetProofDataSet("ds.txt");
  plugin->SetProofDataSet("g4g.txt");
  // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
  plugin->SetProofReset(0);
  // May limit the number of workers per slave. If used with SetNproofWorkers, SetParallel(nproofworkers) will be called after connection
  //plugin->SetNproofWorkersPerSlave(1);
  // May request connection to alien upon connection to grid
  //     plugin->SetProofConnectGrid(kTRUE);
    
  // plugin->SetNproofWorkers(51);
  // May use a specific version of root installed in proof
  //     plugin->SetRootVersionForProof("current");
  // May set the aliroot mode. Check http://aaf.cern.ch/node/83
  plugin->SetAliRootMode("default"); // Loads AF libs by default
  // May request ClearPackages (individual ClearPackage not supported)
  plugin->SetClearPackages(kFALSE);
  // Plugin test mode works only providing a file containing test file locations
  //    plugin->SetFileForTestMode("AOD.txt");
  plugin->SetFileForTestMode("test.txt");
  //++++++++++++++ end PROOF ++++++++++++++++
  return plugin; 
} 
