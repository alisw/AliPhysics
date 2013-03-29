AliAnalysisGrid* CreateAlienHandlerCaloEtSim(TString outputDir, TString outputName, const char * pluginRunMode, int production, Bool_t isPHOS, Bool_t ispp,Bool_t isData)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output results from a previous session
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  // plugin->SetRunMode("full");  // VERY IMPORTANT - DECRIBED BELOW
  // plugin->SetRunMode("test");  // VERY IMPORTANT - DECRIBED BELOW
  plugin->SetRunMode(pluginRunMode);  // VERY IMPORTANT - DECRIBED BELOW
  cout<<"Running in "<<pluginRunMode<<" mode"<<endl;

  // Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-34-02-1");
   plugin->SetAliROOTVersion("v5-04-34-AN");
  // Declare input data to be processed.

  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  //   plugin->SetGridDataDir("/alice/sim/LHC10a18");
  // Set data search pattern
  //   plugin->SetDataPattern("*ESDs.root");  // simulated, tags not used
  //   plugin->SetDataPattern("*ESDs/pass4/*ESDs.root"); // real data check reco pass and data base directory
  //   plugin->SetRunPrefix("000");   // real data
  //   plugin->SetDataPattern("*tag.root");  // Use ESD tags (same applies for AOD's)
  // ...then add run numbers to be considered
  //   plugin->AddRunNumber(125020);    // simulated
  //   plugin->AddRunNumber(104065);  // real data

   //plugin->SetGridDataDir("/alice/sim/LHC10d4");
   //plugin->AddRunNumber("120741");//smallest of the above
   if(ispp && production==0){
     //pp
     plugin->SetGridDataDir("/alice/sim/2012/LHC12a15e/");
     plugin->AddRunNumber(169838);
   }
   else{
    if(isData){//185 jobs
      cout<<"Running over data"<<endl;
      plugin->SetGridDataDir("/alice/data/2010/LHC10h");//PbPb data
      plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
      plugin->SetRunPrefix("000");   // real data
      plugin->AddRunNumber(139465);
      outputDir = outputDir + "LHC10hPass2";
    }
    else{
      if(production==0){
	//Standard
	if(isPHOS){
	  plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	}
	else{
	  outputDir = outputDir + "LHC11a4_bis";
	  plugin->SetGridDataDir("/alice/sim/LHC11a4_bis");
	  plugin->SetGridDataDir("/alice/sim/LHC11a4_bis");
	}
	plugin->AddRunNumber(139465);
	plugin->AddRunNumber(139470);
	plugin->AddRunNumber(137366);
	plugin->AddRunNumber(137161);
      }
      if(production==1){
	//if(!isPHOS){
	  outputDir = outputDir + "LHC11a10a_bis";
	  plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	  //}
// 	        plugin->AddRunNumber(137366);
// 	        plugin->AddRunNumber(137161);
// 	        plugin->AddRunNumber(139470);
		plugin->AddRunNumber(139465);//probably our focus now
      }
      if(production==2){
	if(!isPHOS){
	  outputDir = outputDir + "LHC11b7";
	  plugin->SetGridDataDir("/alice/sim/LHC11b7");
	}
	plugin->AddRunNumber(137549);
	plugin->AddRunNumber(138200);
      }
      if(production==3){
	if(!isPHOS){
	  outputDir = outputDir + "LHC11a10a";
	  plugin->SetGridDataDir("/alice/sim/LHC11a10a");
	}
	plugin->AddRunNumber(139470);
      }
    }
   }

   plugin->SetDataPattern("*ESDs.root");
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("tag.xml");
  // plugin->AddDataFile("wn.xml"); // test
  // file generated with:  find -x tag /alice/sim/LHC10d1/117222/* AliESDs.root > tag.xml

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(outputDir.Data());
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime IN THE SAME ORDER THEY ARE LISTED
  //plugin->SetAdditionalRootLibs("libPHOSUtils.so libEMCALUtils.so libPWG4CaloCalib.so libPWG4PartCorrBase.so libPWG4PartCorrDep.so");
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisEtSelector.cxx AliAnalysisEtSelectorPhos.cxx AliAnalysisEtSelectorEmcal.cxx AliAnalysisEtTrackMatchCorrections.cxx AliAnalysisEtRecEffCorrection.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisEmEtMonteCarlo.cxx AliAnalysisEmEtReconstructed.cxx AliAnalysisTaskTotEt.cxx");
  //plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisEtSelector.cxx AliAnalysisEtSelectorPhos.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskTotEt.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libPHOSUtils.so libTENDER.so libTENDERSupplies.so AliAnalysisEtCuts.cxx AliAnalysisEtCuts.h AliAnalysisHadEtCorrections.cxx AliAnalysisHadEtCorrections.h AliAnalysisEtCommon.cxx AliAnalysisEtCommon.h AliAnalysisEtSelector.cxx AliAnalysisEtSelector.h AliAnalysisEtSelectorPhos.cxx AliAnalysisEtSelectorPhos.h AliAnalysisEtSelectorEmcal.cxx AliAnalysisEtSelectorEmcal.h AliAnalysisEtTrackMatchCorrections.cxx AliAnalysisEtTrackMatchCorrections.h AliAnalysisEtRecEffCorrection.cxx AliAnalysisEtRecEffCorrection.h AliAnalysisEt.cxx AliAnalysisEt.h AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarlo.h AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloPhos.h AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtMonteCarloEmcal.h AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructed.h AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedPhos.h AliAnalysisEtReconstructedEmcal.cxx AliAnalysisEtReconstructedEmcal.h AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskTransverseEnergy.h AliAnalysisEmEtMonteCarlo.cxx AliAnalysisEmEtMonteCarlo.h AliAnalysisEmEtReconstructed.cxx AliAnalysisEmEtReconstructed.h AliAnalysisTaskTotEt.cxx AliAnalysisTaskTotEt.h badchannels.root corrections.root calocorrections.root ConfigEtMonteCarlo.C ConfigEtReconstructed.C");
  plugin->SetExecutableCommand("aliroot -b -q");
  // add extra include files/path
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");

  // No need for output file names. Procedure is automatic. <-- not true
  //plugin->SetDefaultOutputs(kFALSE);
  //plugin->SetOutputFiles(outputName.Data());
  //plugin->SetOutputFiles("Et.ESD.sim.EMCAL.root event_stat.root");
  // No need define the files to be archived. Note that this is handled automatically by the plugin.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
  plugin->SetAnalysisMacro("DavidEtAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore). The optimum for an analysis
  // is correlated with the run time - count few hours TTL per job, not minutes !
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  //plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(20000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskEt.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1); 
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");

  return plugin;
} 
