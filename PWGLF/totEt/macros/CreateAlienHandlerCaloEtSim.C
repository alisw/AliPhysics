AliAnalysisGrid* CreateAlienHandlerCaloEtSim(TString outputDir, TString outputName, const char * pluginRunMode, int production, Bool_t isPHOS, Bool_t ispp,Bool_t isData, Int_t runnum, Bool_t runCompiledVersion)
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
   plugin->SetROOTVersion("v5-34-08-6");
   plugin->SetAliROOTVersion("vAN-20140623");
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


   plugin->SetDataPattern("*ESDs.root");
   //plugin->SetGridDataDir("/alice/sim/LHC10d4");
   //plugin->AddRunNumber("120741");//smallest of the above
   if(ispp){
     if( production==0){
       //pp
       plugin->SetGridDataDir("/alice/sim/2012/LHC12a15e/");
       plugin->AddRunNumber(169838);
     }
     if( production==1 || production==2 || production==3){
       cout<<"I am here! "<<endl;
       //pp
       if(production==1){
	 plugin->SetGridDataDir(" /alice/sim/LHC11b1a/");//nominal
	 outputDir = outputDir + "NominalLHC11b1a";
       }
       if(production==2){
	 plugin->SetGridDataDir(" /alice/sim/LHC11b1b/");//high material budget
	 outputDir = outputDir + "HighLHC11b1b";
       }
       if(production==3){
	 plugin->SetGridDataDir(" /alice/sim/LHC11b1c/");//low material budget
	 outputDir = outputDir + "LowLHC11b1c";
       }
       plugin->AddRunNumber(121040);//all runs in these productions with good EMC and PHOS and global status
       plugin->AddRunNumber(121039);//
       plugin->AddRunNumber(118558);
       plugin->AddRunNumber(118518);
       plugin->AddRunNumber(118506);
       //Additional runs which may be used 118561, 118560, 118556 emc bad parts
       //118512, 118507 status unknown
     }
   }
   else{
    if(isData){//185 jobs
	 cout<<"Running over data"<<endl;
       if(production==1){
	 plugin->SetGridDataDir("/alice/data/2010/LHC10h");//PbPb data
	 plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
	 plugin->SetRunPrefix("000");   // real data
	 outputDir = outputDir + "LHC10hPass2";
	 if(runnum==0){
	   plugin->AddRunNumber(139465);
// 	   plugin->AddRunNumber(138442);
// 	   plugin->AddRunNumber(138364);
// 	   plugin->AddRunNumber(138396);
// 	   plugin->AddRunNumber(137722);
	   //outputDir = outputDir + "Run139465";
	 }
	 if(runnum==1){
	   plugin->AddRunNumber(138442);
	   outputDir = outputDir + "Run138442";
	 }
	 if(runnum==2){
	   plugin->AddRunNumber(138364);
	   outputDir = outputDir + "Run138364";
	 }
	 if(runnum==3){
	   plugin->AddRunNumber(138534);
	   outputDir = outputDir + "Run138534";
	 }
	 if(runnum==4){
	   plugin->AddRunNumber(138275);
	   outputDir = outputDir + "Run138275";
	 }
       }
       if(production==2){
	 plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");//PbPb data
	 plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
	 plugin->SetRunPrefix("000");   // real data
	 //plugin->AddRunNumber(169099);
	 if(runnum==0){
	   plugin->AddRunNumber(168464);
	 }
	 outputDir = outputDir + "LHC11hPass2";
	 if(runnum==1){
	   plugin->AddRunNumber(169588);
	   outputDir +="169588";
	 }
	 if(runnum==2){
	   plugin->AddRunNumber(170268);
	   outputDir +="170268";
	 }
	 if(runnum==3){
	   plugin->AddRunNumber(170207);
	   outputDir +="170207";
	 }
	 if(runnum==4){
	   plugin->AddRunNumber(168512);
	   outputDir +="168512";
	 }
	 if(runnum==5){
	   plugin->AddRunNumber(170311);
	   outputDir +="170311";
	 }
       }
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
	//cout<<"I am here! Line 93 "<<endl;
	//if(!isPHOS){
	  outputDir = outputDir + "LHC11a10a_bis";
	  plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	  //}
// 	        plugin->AddRunNumber(137366);
// 	        plugin->AddRunNumber(137161);
// 	        plugin->AddRunNumber(139470);
		plugin->AddRunNumber(139465);//probably our focus now
// 	 plugin->AddRunNumber(138442);
// 	 plugin->AddRunNumber(138364);
// 	 plugin->AddRunNumber(138534);
// 	 plugin->AddRunNumber(138275);
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
      if(production==4){
	outputDir = outputDir + "LHC13d2";
	plugin->SetGridDataDir(" /alice/sim/2013/LHC13d2");
	plugin->AddRunNumber(139465);//probably our focus now
      }
      if(production==5){//2011 production - warning only 0-10%
	cout<<"I am here setting grid data dir"<<endl;
	outputDir = outputDir + "LHC12d3";
	plugin->SetGridDataDir("/alice/sim/2012/LHC12d3");
	plugin->AddRunNumber(168464);//probably our focus now
      }
      if(production==6){//2011 production - 0-10%
	cout<<"I am here setting grid data dir"<<endl;
	outputDir = outputDir + "LHC13e1a";
	plugin->SetGridDataDir("/alice/sim/2013/LHC13e1a");
	plugin->AddRunNumber(168464);//probably our focus now
      }
      if(production==7){//2011 production - 10-50%
	cout<<"I am here setting grid data dir"<<endl;
	outputDir = outputDir + "LHC13e1b";
	plugin->SetGridDataDir("/alice/sim/2013/LHC13e1b");
	plugin->AddRunNumber(168464);//probably our focus now
      }
      if(production==8){//2011 production - 50-90%
	cout<<"I am here setting grid data dir"<<endl;
	outputDir = outputDir + "LHC13e1c";
	plugin->SetGridDataDir("/alice/sim/2013/LHC13e1c");
	plugin->AddRunNumber(168464);//probably our focus now
      }


    }
   }

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
  if(runCompiledVersion){
    plugin->SetAdditionalLibs("libPHOSUtils.so libTender.so libTenderSupplies.so libPWGTools.so libPWGEMCAL.so badchannels.root libPWGLFtotEt.so corrections.root calocorrections.root ConfigEtMonteCarlo.C ConfigEtReconstructed.C");
  }
  else{
    plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisEtSelector.cxx AliAnalysisEtSelectorPhos.cxx AliAnalysisEtSelectorEmcal.cxx AliAnalysisEtTrackMatchCorrections.cxx AliAnalysisEtRecEffCorrection.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisEmEtMonteCarlo.cxx AliAnalysisEmEtReconstructed.cxx AliAnalysisTaskTotEt.cxx");
  //plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisEtSelector.cxx AliAnalysisEtSelectorPhos.cxx AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedEmcal.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskTotEt.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("libPHOSUtils.so libTender.so libTenderSupplies.so libPWGTools.so libPWGEMCAL.so AliAnalysisEtCuts.cxx AliAnalysisEtCuts.h AliAnalysisHadEtCorrections.cxx AliAnalysisHadEtCorrections.h AliAnalysisEtCommon.cxx AliAnalysisEtCommon.h AliAnalysisEtSelector.cxx AliAnalysisEtSelector.h AliAnalysisEtSelectorPhos.cxx AliAnalysisEtSelectorPhos.h AliAnalysisEtSelectorEmcal.cxx AliAnalysisEtSelectorEmcal.h AliAnalysisEtTrackMatchCorrections.cxx AliAnalysisEtTrackMatchCorrections.h AliAnalysisEtRecEffCorrection.cxx AliAnalysisEtRecEffCorrection.h AliAnalysisEt.cxx AliAnalysisEt.h AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarlo.h AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtMonteCarloPhos.h AliAnalysisEtMonteCarloEmcal.cxx AliAnalysisEtMonteCarloEmcal.h AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructed.h AliAnalysisEtReconstructedPhos.cxx AliAnalysisEtReconstructedPhos.h AliAnalysisEtReconstructedEmcal.cxx AliAnalysisEtReconstructedEmcal.h AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskTransverseEnergy.h AliAnalysisEmEtMonteCarlo.cxx AliAnalysisEmEtMonteCarlo.h AliAnalysisEmEtReconstructed.cxx AliAnalysisEmEtReconstructed.h AliAnalysisTaskTotEt.cxx AliAnalysisTaskTotEt.h badchannels.root corrections.root calocorrections.root ConfigEtMonteCarlo.C ConfigEtReconstructed.C");
  }
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
  plugin->SetSplitMaxInputFileNumber(50);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //plugin->SetMaxInitFailed(50);
  // Optionally resubmit threshold.
  //plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("TaskEt.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1); 
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");


  plugin->SetMergeViaJDL();

  return plugin;
} 
