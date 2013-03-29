AliAnalysisGrid* CreateAlienHandlerHadEt(Int_t dataset, Bool_t data, Bool_t test, Int_t material, Bool_t altV0Scale = kFALSE, bool runCompiledVersion = kFALSE, int simflag = 0)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output results from a previous session
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  if(test){
    plugin->SetRunMode("test");  // VERY IMPORTANT - DECRIBED BELOW
  }
  else{
    plugin->SetRunMode("full");  // VERY IMPORTANT - DECRIBED BELOW
  }
  //needed for local testing?
  //plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
  // Set versions of used packages 
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-34-02-1");
   plugin->SetAliROOTVersion("v5-04-34-AN");
  // Declare input data to be processed.

   plugin->AddIncludePath("-I$ALICE_ROOT/PWGUD/base");
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

  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("tag.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  if(dataset==20100){//PbPb 2.76 TeV
    if(data){//185 jobs
      cout<<"Running over data"<<endl;
      plugin->SetGridDataDir("/alice/data/2010/LHC10h");//PbPb data
      plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
      plugin->SetRunPrefix("000");   // real data
    }
    else{
      cout<<"Running over MC"<<endl;
      if(material==0){
	  plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	//plugin->SetGridDataDir("/alice/sim/LHC11a4_bis");//PbPb simulation
      }
      if(material==1){//7% more material
	plugin->SetGridDataDir("/alice/sim/LHC10h9");//PbPb simulation
      }
      if(material==-1){//7% less material
	plugin->SetGridDataDir("/alice/sim/LHC10h10");//PbPb simulation
      }
      if(material==10){//OK it's a cheat but this runs on DPM jet
	plugin->SetGridDataDir("/alice/sim/LHC10h2");//PbPb simulation
      }
      if(material==11){//OK it's a cheat but this runs on AMPT
	plugin->SetGridDataDir(" /alice/sim/LHC11a9a");//PbPb simulation
      }
      plugin->SetDataPattern("*ESDs.root");
      //plugin->SetGridWorkingDir("etPbPbSim");
    }
    plugin->AddRunNumber(139465);
//     plugin->AddRunNumber(137366);
//     plugin->AddRunNumber(137161);
  }
  else{
    if(dataset==2009){//pp 900 GeV
      if(data){//only 233 jobs!
	cout<<"Running over data"<<endl;
	plugin->SetGridDataDir("/alice/data/2010/LHC10c");//PbPb data
	//plugin->SetDataPattern("*ESDs/pass3/*ESDs.root");
	plugin->SetDataPattern("*ESDs/pass3/*ESDs.root");
	plugin->SetRunPrefix("000");   // real data
      }
      else{//sim over 1500 jobs, some get killed because it's above quota
	cout<<"Running over MC"<<endl;
	if(simflag==0){
	  if(material==0){
	    plugin->SetGridDataDir("/alice/sim/LHC11b1a");//PbPb simulation
	  }
	  if(material==1){//10% more material budget
	    plugin->SetGridDataDir("/alice/sim/LHC11b1b");//PbPb simulation
	  }
	  if(material==-1){//10% less material budget
	    plugin->SetGridDataDir("/alice/sim/LHC11b1c");//PbPb simulation
	  }
	}
	if(simflag==1)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h1a");//PYTHIA
	if(simflag==2)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h1b");//PHOJET
	if(simflag==3)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h1c");//PYTHIA Flat
	plugin->SetDataPattern("*ESDs.root");
      }
      plugin->AddRunNumber(118506);
//       plugin->AddRunNumber(121040);
//       plugin->AddRunNumber(121039);
//       plugin->AddRunNumber(118561);
//       plugin->AddRunNumber(118560);
//       plugin->AddRunNumber(118558);
//       plugin->AddRunNumber(118557);
//       plugin->AddRunNumber(118556);
//       plugin->AddRunNumber(118518);
//       plugin->AddRunNumber(118512);
//       plugin->AddRunNumber(118507);
//       plugin->AddRunNumber(118506);
//if(data){
// 	plugin->AddRunNumber(118504);
// 	plugin->AddRunNumber(118503);
//    }
    }
  
    if(dataset==20111){//pp 2.76 TeV 
      if(data){//257 jobs
	cout<<"Running over data"<<endl;
	plugin->SetGridDataDir("/alice/data/2011/LHC11a");//
	plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
	plugin->SetRunPrefix("000");   // real data
      }
      else{//sim - 332 jobs
	cout<<"Running over MC"<<endl;
	if(simflag==0)	plugin->SetGridDataDir("/alice/sim/LHC11b10a");//
	if(simflag==1)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h5a");//PYTHIA
	if(simflag==2)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h5b");//PHOJET
	if(simflag==3)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h5c");//PYTHIA Flat
	plugin->SetDataPattern("*ESDs.root");
      }
//       plugin->AddRunNumber(146860);
//       plugin->AddRunNumber(146859);
//       plugin->AddRunNumber(146856);
//       plugin->AddRunNumber(146824);
//       plugin->AddRunNumber(146817);
//       plugin->AddRunNumber(146806);
//       plugin->AddRunNumber(146805);
//       plugin->AddRunNumber(146804);
//       plugin->AddRunNumber(146803);
//       plugin->AddRunNumber(146802);
//       plugin->AddRunNumber(146801);
//       plugin->AddRunNumber(146748);
//       plugin->AddRunNumber(146747);
//       plugin->AddRunNumber(146746);

  //     plugin->AddRunNumber(146860);
//       plugin->AddRunNumber(146859);
//       plugin->AddRunNumber(146858);
//       plugin->AddRunNumber(146857);
//       plugin->AddRunNumber(146856);
      //        plugin->AddRunNumber(146824);
      //        if(data){//these productions are not yet done for MC
      // 	 plugin->AddRunNumber(146817);
      // 	 plugin->AddRunNumber(146807);
      // 	 plugin->AddRunNumber(146806);
      // 	 plugin->AddRunNumber(146805);
      // 	 plugin->AddRunNumber(146804);
      // 	 plugin->AddRunNumber(146803);
      // 	 plugin->AddRunNumber(146802);
      //        }
      plugin->AddRunNumber(146805);

    }
    if(dataset==2010){//pp 7 TeV
      if(data){//data - 569 jobs
	cout<<"Running over 7 TeV data"<<endl;
	plugin->SetGridDataDir("/alice/data/2010/LHC10e");//PbPb data
	plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
	cout<<"Setting run prefix to be 000"<<endl;
	plugin->SetRunPrefix("000");   // real data
      }
      else{//sim- 346 jobs
	if(simflag==0)	plugin->SetGridDataDir("/alice/sim/LHC10e20");//
	if(simflag==1)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h4a");//PYTHIA
	if(simflag==2)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h4b");//PHOJET
	if(simflag==3)	plugin->SetGridDataDir("/alice/sim/2011/LHC11h4c");//PYTHIA Flat
	plugin->SetDataPattern("*ESDs.root");
      }
      plugin->AddRunNumber("130795");
// 	plugin->AddRunNumber("130840");
// 	plugin->AddRunNumber("130834");
// 	plugin->AddRunNumber("130833");
// 	plugin->AddRunNumber("130831");
// 	plugin->AddRunNumber("130804");
// 	plugin->AddRunNumber("130803");
// 	plugin->AddRunNumber("130802");
// 	plugin->AddRunNumber("130799");
// 	plugin->AddRunNumber("130798");
// 	plugin->AddRunNumber("130795");
    }

    if(dataset==2012){//pp 8 TeV
      if(data){//data - 569 jobs
	cout<<"Running over 8 TeV data"<<endl;
	plugin->SetGridDataDir("/alice/data/2012/LHC12b");//PbPb data
	plugin->SetDataPattern("*ESDs/pass1/*ESDs.root");
	cout<<"Setting run prefix to be 000"<<endl;
	plugin->SetRunPrefix("000");   // real data
      }
      else{//sim- 346 jobs
	plugin->SetGridDataDir("/alice/sim/2012/LHC12c1b");//
      }
      plugin->AddRunNumber("178030");
    }
    if(dataset==2013){//pPb 
      if(data){//data - 569 jobs
	cout<<"Running over 8 TeV data"<<endl;
	plugin->SetGridDataDir("/alice/data/2013/LHC13b");//PbPb data
	plugin->SetDataPattern("*ESDs/pass2/*ESDs.root");
	cout<<"Setting run prefix to be 000"<<endl;
	plugin->SetRunPrefix("000");   // real data
      }
      else{//sim- 346 jobs
	plugin->SetGridDataDir(" /alice/sim/2013/LHC13b3");//
      }
      plugin->AddRunNumber("195483");
    }
  }


  if(dataset==20100){//PbPb 2.76 TeV
    if(data){
      plugin->SetGridWorkingDir("etPbPbData");
    }
    else{
      if(material==0){plugin->SetGridWorkingDir("etPbPbSim");}
      if(material==-1) plugin->SetGridWorkingDir("etPbPbSimMatBudLow");
      if(material==1) plugin->SetGridWorkingDir("etPbPbSimMatBudHigh");
      if(material==10)  plugin->SetGridWorkingDir("etPbPbSimDPMJET");
       if(material==11)  plugin->SetGridWorkingDir("etPbPbSimAMPT");
    }
  }
  else{
    if(dataset==2009){//pp 900 GeV
      if(data){
	plugin->SetGridWorkingDir("etpp900GeVData");
      }
      else{
	if(simflag==0){
	  if(material==0) plugin->SetGridWorkingDir("etpp900GeVSim");
	  if(material==-1) plugin->SetGridWorkingDir("etpp900GeVSimMatBudLow");
	  if(material==1) plugin->SetGridWorkingDir("etpp900GeVSimMatBudHigh");
	}
	if(simflag==1)	plugin->SetGridWorkingDir("etpp900GeVSimPYTHIA");//PYTHIA
	if(simflag==2)	plugin->SetGridWorkingDir("etpp900GeVSimPHOJET");//PHOJET
	if(simflag==3)	plugin->SetGridWorkingDir("etpp900GeVSimPYTHIAFLAT");//PYTHIA Flat
      }
    }
    if(dataset==20111){//pp 2.76 TeV
      if(data){
	plugin->SetGridWorkingDir("etpp276TeVData");
      }
      else{
	if(altV0Scale) plugin->SetGridWorkingDir("etpp276TeVSimAlt");
	else{
	  // plugin->SetGridWorkingDir("etpp276TeVSim");
	  
	  if(simflag==0)	plugin->SetGridWorkingDir("etpp276TeVSim");//
	  if(simflag==1)	plugin->SetGridWorkingDir("etpp276TeVSimPYTHIA");//PYTHIA
	  if(simflag==2)	plugin->SetGridWorkingDir("etpp276TeVSimPHOJET");//PHOJET
	  if(simflag==3)	plugin->SetGridWorkingDir("etpp276TeVSimPYTHIAFLAT");//PYTHIA Flat
	}
      }
    }
    if(dataset==2010){//pp 7 TeV
      if(data){
	plugin->SetGridWorkingDir("etpp7TeVData");
      }
      else{
	if(simflag==0)	plugin->SetGridWorkingDir("etpp7TeVSim");//
	if(simflag==1)	plugin->SetGridWorkingDir("etpp7TeVSimPYTHIA");//PYTHIA
	if(simflag==2)	plugin->SetGridWorkingDir("etpp7TeVSimPHOJET");//PHOJET
	if(simflag==3)	plugin->SetGridWorkingDir("etpp7TeVSimPYTHIAFLAT");//PYTHIA Flat
      }
    }
    if(dataset==2012){//pp 8 TeV
      if(data){
	plugin->SetGridWorkingDir("etpp8TeVData");
      }
      else{
	plugin->SetGridWorkingDir("etpp8TeVSim");
      }
    }
    if(dataset==2013){//pPb
      if(data){
	plugin->SetGridWorkingDir("etpPb5TeVData");
      }
      else{
	plugin->SetGridWorkingDir("etpPb5TeVSim");
      }
    }
  }


  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  //plugin->SetGridWorkingDir("et");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime IN THE SAME ORDER THEY ARE LISTED
  // using ACLiC on the worker nodes.
  //plugin->SetAnalysisSource("AliAnalysisTaskHadEt.cxx");
  //plugin->SetAnalysisSource("AliAnalysisEt.cxx AliAnalysisEtMonteCarlo.cxx AliAnalysisEtMonteCarloPhos.cxx AliAnalysisEtReconstructed.cxx AliAnalysisEtReconstructedPhos.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskHadEt.cxx AliAnalysisTaskTotEt.cxx");
  //TString sourcefiles = "AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskHadEt.cxx";
  //plugin->SetAnalysisSource(sourcefiles.Data());
  if(!runCompiledVersion){
    plugin->SetAnalysisSource("AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.cxx AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisTaskHadEt.cxx");
  }
   
  //cout<<"Setting source files "<<sourcefiles<<endl;
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  //TString additionallibs = "AliAnalysisEtCuts.h AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.h AliAnalysisHadEtCorrections.cxx  AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisEtCommon.h AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskHadEt.cxx AliAnalysisHadEt.h AliAnalysisHadEtMonteCarlo.h AliAnalysisHadEtReconstructed.h AliAnalysisTaskHadEt.h  AliAnalysisEtSelectionContainer.h AliAnalysisEtSelectionHandler.h AliAnalysisTaskTransverseEnergy.h corrections.root ConfigHadEtAnalysis.C ConfigHadEtMonteCarlo.C ConfigHadEtReconstructed.C physicsSelections.root";
  //TString additionallibs = "AliAnalysisEtCuts.h AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.h AliAnalysisHadEtCorrections.cxx  AliAnalysisEtSelectionContainer.cxx AliAnalysisEtSelectionHandler.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisEtCommon.h AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskHadEt.cxx AliAnalysisHadEt.h AliAnalysisHadEtMonteCarlo.h AliAnalysisHadEtReconstructed.h AliAnalysisTaskHadEt.h  AliAnalysisEtSelectionContainer.h AliAnalysisEtSelectionHandler.h AliAnalysisTaskTransverseEnergy.h physicsSelections.root ConfigHadEtMonteCarlo.C  ConfigHadEtReconstructed.C corrections.root";
  //plugin->SetAdditionalLibs(additionallibs.Data());
  if(!runCompiledVersion){
    plugin->SetAdditionalLibs( "AliAnalysisEtCuts.h AliAnalysisEtCuts.cxx AliAnalysisHadEtCorrections.h AliAnalysisHadEtCorrections.cxx AliAnalysisTaskTransverseEnergy.cxx AliAnalysisEtCommon.h AliAnalysisEtCommon.cxx AliAnalysisHadEt.cxx AliAnalysisHadEtMonteCarlo.cxx AliAnalysisHadEtReconstructed.cxx AliAnalysisTaskHadEt.cxx AliAnalysisHadEt.h AliAnalysisHadEtMonteCarlo.h AliAnalysisHadEtReconstructed.h AliAnalysisTaskHadEt.h AliAnalysisTaskTransverseEnergy.h ConfigHadEtMonteCarlo.C  ConfigHadEtReconstructed.C corrections.root libPWGUDbase.so");
  }
  else{
    plugin->SetAdditionalLibs( "ConfigHadEtMonteCarlo.C  ConfigHadEtReconstructed.C corrections.root libPWGUDbase.so libPWGLFtotEt.so");
  }
  // No need for output file names. Procedure is automatic. <-- not true
  //plugin->SetDefaultOutputs(kFALSE);
  //plugin->SetOutputFiles("Et.ESD.new.sim.root");
  // No need define the files to be archived. Note that this is handled automatically by the plugin.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
  plugin->SetAnalysisMacro("ChristinesEtAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore). The optimum for an analysis
  // is correlated with the run time - count few hours TTL per job, not minutes !
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //plugin->SetMaxInitFailed(5);
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
  plugin->SetTerminateFiles("event_stat.root") ;
  plugin->SetKeepLogs();
  return plugin;
} 
